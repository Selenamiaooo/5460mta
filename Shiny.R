library(shiny)
library(bslib)
library(dplyr)
library(ggplot2)
library(plotly)
library(scales)
library(lubridate)
library(broom)
library(leaflet)
library(readr)
library(tibble)
library(sf)

# -------------------------
# Data loading (hard fail if missing)
# -------------------------
stopifnot(file.exists("data/monthly_stats.rds"))
stopifnot(file.exists("data/station_effects.rds"))

monthly_stats <- readRDS("data/monthly_stats.rds")
station_effects <- readRDS("data/station_effects.rds")

# Ensure Date type
if (!inherits(monthly_stats$month_date, "Date")) {
  monthly_stats <- monthly_stats %>% mutate(month_date = as.Date(month_date))
}

default_start <- min(monthly_stats$month_date, na.rm = TRUE)
default_end   <- max(monthly_stats$month_date, na.rm = TRUE)

# -------------------------
# Helpers
# -------------------------
compute_beta_ci <- function(df) {
  df <- df %>% arrange(month_date)
  if (nrow(df) < 6) {
    return(list(beta = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_, r2 = NA_real_))
  }
  m <- lm(mean_ridership ~ time_index, data = df)
  tr <- tidy(m, conf.int = TRUE) %>% filter(term == "time_index")
  gl <- glance(m)
  list(beta = tr$estimate, lo = tr$conf.low, hi = tr$conf.high, p = tr$p.value, r2 = gl$r.squared)
}

compute_high_share <- function(df, threshold) {
  df %>%
    mutate(sd_safe = ifelse(is.na(sd_ridership) | sd_ridership <= 0, NA_real_, sd_ridership)) %>%
    mutate(share_high = ifelse(
      is.na(sd_safe),
      NA_real_,
      pmax(0, pmin(1, 1 - pnorm((threshold - mean_ridership) / sd_safe)))
    )) %>%
    select(-sd_safe)
}

fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("< 0.001")
  sprintf("%.3f", p)
}

safe_read_stations <- function() {
  candidates <- c("stationsmap.csv", "data/stationsmap.csv")
  path <- candidates[file.exists(candidates)][1]
  if (is.na(path)) return(NULL)
  
  df <- suppressWarnings(readr::read_csv(path, show_col_types = FALSE))
  names(df) <- trimws(names(df))
  
  if (!("Longitude" %in% names(df)) && ("longitude" %in% names(df))) df$Longitude <- df$longitude
  if (!("Latitude"  %in% names(df)) && ("latitude"  %in% names(df))) df$Latitude  <- df$latitude
  
  if (!("borough" %in% names(df)) && ("Borough" %in% names(df))) df$borough <- df$Borough
  
  if (!("station_complex_id" %in% names(df)) && ("STATION_COMPLEX_ID" %in% names(df))) {
    df$station_complex_id <- df$STATION_COMPLEX_ID
  }
  
  if (!("station_complex" %in% names(df))) {
    if ("station" %in% names(df)) df$station_complex <- df$station
    if ("station_name" %in% names(df)) df$station_complex <- df$station_name
    if ("STATION_COMPLEX" %in% names(df)) df$station_complex <- df$STATION_COMPLEX
  }
  
  if ("station_complex_id" %in% names(df)) df$station_complex_id <- as.character(df$station_complex_id)
  if ("station_complex" %in% names(df)) df$station_complex <- as.character(df$station_complex)
  
  if ("borough" %in% names(df)) {
    df$borough <- as.character(df$borough)
    df$borough <- ifelse(is.na(df$borough), NA_character_, tools::toTitleCase(tolower(df$borough)))
  }
  
  df
}

stations <- safe_read_stations()

# -------------------------
# Spatial objects (sf)
# -------------------------
stations_sf <- NULL
nn_distances <- NULL

if (!is.null(stations) && "Longitude" %in% names(stations) && "Latitude" %in% names(stations)) {
  stations_sf <- stations %>%
    filter(!is.na(Longitude), !is.na(Latitude)) %>%
    st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%
    st_transform(crs = 32618)  # UTM 18N for distance in meters
  
  # Join station_effects
  eff <- station_effects %>% mutate(station_complex_id = as.character(station_complex_id))
  if ("station_complex_id" %in% names(stations_sf)) {
    stations_sf <- stations_sf %>%
      mutate(station_complex_id = as.character(station_complex_id)) %>%
      left_join(eff, by = "station_complex_id", suffix = c("", ".eff"))
  } else if ("station_complex" %in% names(stations_sf)) {
    stations_sf <- stations_sf %>%
      left_join(eff, by = "station_complex", suffix = c("", ".eff"))
  }
  
  # Compute nearest-neighbor distance for each station
  dist_mat <- st_distance(stations_sf)
  diag(dist_mat) <- Inf
  stations_sf$nn_dist_m <- apply(dist_mat, 1, min) %>% as.numeric()
  stations_sf$nn_dist_km <- stations_sf$nn_dist_m / 1000
  
  # Find index of nearest neighbor
  stations_sf$nn_idx <- apply(dist_mat, 1, which.min)
  if ("station_complex" %in% names(stations_sf)) {
    stations_sf$nn_name <- stations_sf$station_complex[stations_sf$nn_idx]
  }
}

# Load study area polygon if available
study_area_sf <- NULL
area_candidates <- c("my_area.geojson", "data/my_area.geojson")
area_path <- area_candidates[file.exists(area_candidates)][1]
if (!is.na(area_path)) {
  study_area_sf <- tryCatch(st_read(area_path, quiet = TRUE) %>% st_transform(32618), error = function(e) NULL)
}

# -------------------------
# Forecast helpers (ITS + counterfactual)
# -------------------------
make_its_frame <- function(df, policy_date) {
  df <- df %>% arrange(month_date)
  
  df <- df %>%
    mutate(
      time_index = row_number(),
      post = as.integer(month_date >= policy_date)
    )
  
  idx0 <- which(df$month_date >= policy_date)[1]
  if (is.na(idx0)) {
    df <- df %>% mutate(time_after = 0L)
  } else {
    df <- df %>% mutate(time_after = pmax(0L, time_index - idx0 + 1L))
  }
  
  df
}

fit_its <- function(df_its) {
  lm(mean_ridership ~ time_index + post + time_after, data = df_its)
}

predict_with_se <- function(model, newdata) {
  pr <- predict(model, newdata = newdata, se.fit = TRUE)
  out <- newdata
  out$fit <- as.numeric(pr$fit)
  out$se  <- as.numeric(pr$se.fit)
  out$lo  <- out$fit - 1.96 * out$se
  out$hi  <- out$fit + 1.96 * out$se
  out
}

# -------------------------
# UI
# -------------------------
ui <- tagList(
  tags$div(
    style = "background: linear-gradient(90deg, #0ea5e9, #22c55e);
             color:white; text-align:center; padding:14px 10px;
             font-size:22px; font-weight:800; letter-spacing:0.2px;",
    "NYC Subway Early Warning Dashboard — Congestion Pricing (MVP)"
  ),
  
  page_navbar(
    title = "Early Warning System",
    theme = bs_theme(version = 5, bootswatch = "flatly", base_font = font_google("Inter")),
    fillable = FALSE,
    
    nav_panel(
      "Overview",
      layout_sidebar(
        sidebar = sidebar(
          dateRangeInput(
            "date_rng",
            "Select date range (monthly)",
            start = default_start,
            end   = default_end,
            min   = default_start,
            max   = default_end
          ),
          sliderInput(
            "thr",
            "High-demand threshold (riders/hour)",
            min = 20, max = 300, value = 100, step = 10
          ),
          selectInput(
            "view_mode",
            "Trend view",
            choices = c("Mean with 95% CI" = "mean_ci", "Mean with ±1 SD bands" = "mean_sd"),
            selected = "mean_ci"
          ),
          hr(),
          p(strong("How to read this page:")),
          p("Chart 1 shows the system-wide monthly mean of hourly subway entries with uncertainty. Chart 2 converts ridership into a simple early-warning pressure signal using your chosen threshold. Hover any point to see exact values.")
        ),
        
        layout_columns(
          col_widths = c(6, 6),
          
          value_box(
            title = "Temporal metric: monthly trend β (95% CI)",
            uiOutput("vb_beta"),
            style = "background: linear-gradient(90deg, #0ea5e9, #22c55e); color: white;"
          ),
          
          value_box(
            title = "Current high-demand share (last month)",
            uiOutput("vb_high"),
            style = "background: linear-gradient(90deg, #f97316, #ef4444); color: white;"
          )
        ),
        
        card(
          card_header("1) System trend over time (interactive)"),
          plotlyOutput("p_trend", height = "420px")
        ),
        
        card(
          card_header("2) Early-warning pressure over time (interactive)"),
          plotlyOutput("p_high", height = "360px")
        ),
        
        card(
          card_header("Plain-language interpretation"),
          uiOutput("txt_overview")
        )
      )
    ),
    
    nav_panel(
      "Map",
      layout_sidebar(
        sidebar = sidebar(
          selectInput(
            "map_metric",
            "Map metric",
            choices = c(
              "Station trend β (riders/hour/month)" = "beta",
              "Model fit R²" = "r2",
              "Volatility (sigma)" = "sigma"
            ),
            selected = "beta"
          ),
          checkboxInput(
            "map_sig",
            "Show only statistically significant trends (95% CI excludes 0)",
            value = FALSE
          ),
          uiOutput("map_borough_ui"),
          hr(),
          p(strong("How to read this map:")),
          p("Each marker represents a station complex from your station coordinate file. Color encodes the selected station-level metric computed from monthly station means (2021–2025). Hover a marker to see the station name and numeric details.")
        ),
        
        card(
          card_header("Station-level distributional risk (interactive map)"),
          leafletOutput("leaf_map", height = "650px")
        )
      )
    ),
    
    nav_panel(
      "Details",
      layout_sidebar(
        sidebar = sidebar(
          sliderInput(
            "thr2",
            "Threshold (sync with Overview)",
            min = 20, max = 300, value = 100, step = 10
          ),
          hr(),
          p(strong("Why these charts exist:")),
          p("They provide alternative views that are easier to read for non-technical audiences, including seasonality and the typical range of the pressure indicator.")
        ),
        
        card(
          card_header("A) Seasonal pattern (heatmap): high-demand share by month-of-year"),
          plotlyOutput("p_heat", height = "380px")
        ),
        
        card(
          card_header("B) Distribution: high-demand share across months"),
          plotlyOutput("p_dist", height = "360px")
        ),
        
        card(
          card_header("Explanation & limits"),
          uiOutput("txt_limits")
        )
      )
    ),
    
    nav_panel(
      "Forecast (12 months)",
      layout_sidebar(
        sidebar = sidebar(
          dateInput(
            "policy_date",
            "Congestion pricing start date",
            value = as.Date("2025-01-05"),
            min = default_start,
            max = default_end
          ),
          hr(),
          uiOutput("fc_group_ui"),
          uiOutput("fc_group_select_ui"),
          hr(),
          p(strong("How to read this page:")),
          p("Black points are observed monthly ridership. Blue line is the projected path if congestion pricing continues (Policy ON). Gray line is the counterfactual projection if the policy had not started (Policy OFF). The shaded ribbons show 95% uncertainty bands.")
        ),
        
        layout_columns(
          col_widths = c(6, 6),
          
          value_box(
            title = "Forecasted Policy Effect (12 mo avg gap)",
            uiOutput("vb_fc_gap"),
            style = "background: linear-gradient(90deg, #0ea5e9, #22c55e); color: white;"
          ),
          
          value_box(
            title = "Forecast horizon",
            HTML("<div style='font-size:26px; font-weight:950;'>12 months</div><div style='font-size:13px; opacity:0.95;'>From last observed month</div>"),
            style = "background: linear-gradient(90deg, #f97316, #ef4444); color: white;"
          )
        ),
        
        card(
          card_header("Observed + 12-month forecast (Policy ON vs counterfactual)"),
          plotlyOutput("p_forecast", height = "520px")
        ),
        
        card(
          card_header("Projected policy gap by month (Policy ON − Policy OFF)"),
          plotlyOutput("p_gap_monthly", height = "320px")
        ),
        
        card(
          card_header("Cumulative projected gap (sum over months)"),
          plotlyOutput("p_gap_cum", height = "320px"),
          uiOutput("txt_gap_explain")
        ),
        
        card(
          card_header("Plain-language summary"),
          uiOutput("txt_forecast")
        )
      )
    ),
    
    nav_panel(
      "Spatial Analysis",
      layout_sidebar(
        sidebar = sidebar(
          selectInput(
            "sp_station",
            "Select a station",
            choices = if (!is.null(stations_sf) && "station_complex" %in% names(stations_sf))
              sort(unique(stations_sf$station_complex)) else "No stations",
            selected = NULL
          ),
          sliderInput(
            "sp_radius_km",
            "Neighborhood radius (km)",
            min = 0.5, max = 5, value = 1.5, step = 0.5
          ),
          selectInput(
            "sp_borough",
            "Borough filter",
            choices = if (!is.null(stations_sf) && "borough" %in% names(stations_sf))
              c("All" = "ALL", sort(unique(na.omit(stations_sf$borough)))) else "All",
            selected = "ALL"
          ),
          hr(),
          p(strong("How to read this page:")),
          p("This page uses sf-powered spatial calculations to analyze station proximity and neighborhood-level ridership trends. Select a station and radius to see which stations fall within that service area and their average trend characteristics.")
        ),
        
        layout_columns(
          col_widths = c(6, 6),
          value_box(
            title = "Spatial metric: mean nearest-neighbor distance",
            uiOutput("vb_nn_dist"),
            style = "background: linear-gradient(90deg, #8b5cf6, #6366f1); color: white;"
          ),
          value_box(
            title = "Spatial average β within radius",
            uiOutput("vb_spatial_beta"),
            style = "background: linear-gradient(90deg, #0ea5e9, #22c55e); color: white;"
          )
        ),
        
        card(
          card_header("1) Station proximity map — nearest-neighbor distance & service area"),
          leafletOutput("sp_map", height = "520px")
        ),
        
        card(
          card_header("2) Nearest-neighbor distance vs. ridership trend (β)"),
          plotlyOutput("sp_scatter", height = "400px")
        ),
        
        card(
          card_header("Spatial neighborhood summary"),
          uiOutput("sp_txt_neighborhood")
        ),
        
        card(
          card_header("Station isolation & clustering interpretation"),
          uiOutput("sp_txt_clustering")
        )
      )
    ),
    
    nav_panel(
      "About this site",
      fluidPage(
        br(),
        card(
          card_header("About This Dashboard"),
          card_body(
            p("This site presents a minimum viable early-warning dashboard for policy blind spots related to NYC congestion pricing. It translates large-scale MTA subway ridership data into stakeholder-readable indicators of system change and distributional risk."),
            p("The intent is that a non-technical viewer can understand what is changing over time, whether the change is statistically meaningful, and where potential pressure is concentrating across the network.")
          )
        ),
        br(),
        card(
          card_header("Data Sources"),
          card_body(
            p("The dashboard uses New York State Open Data’s MTA Subway Hourly Ridership datasets and aggregates them to monthly summaries for stability and performance."),
            p("The interactive map uses a station coordinate reference file named stationsmap.csv (located either in the project root or in data/).")
          )
        ),
        br(),
        card(
          card_header("How to Run (Posit Cloud)"),
          card_body(
            p("Put this file as app.R in the project root. Ensure data/monthly_stats.rds and data/station_effects.rds exist. Then click Run App.")
          )
        ),
        br()
      )
    )
  )
)

# -------------------------
# Server
# -------------------------
server <- function(input, output, session) {
  
  observeEvent(input$thr, {
    updateSliderInput(session, "thr2", value = input$thr)
  }, ignoreInit = TRUE)
  
  filtered_monthly <- reactive({
    req(input$date_rng)
    monthly_stats %>%
      filter(month_date >= input$date_rng[1], month_date <= input$date_rng[2]) %>%
      arrange(month_date)
  })
  
  beta_obj <- reactive({
    df <- filtered_monthly() %>% mutate(time_index = row_number())
    compute_beta_ci(df)
  })
  
  high_df <- reactive({
    compute_high_share(filtered_monthly(), input$thr)
  })
  
  output$vb_beta <- renderUI({
    b <- beta_obj()
    if (is.na(b$beta)) {
      return(HTML("<div style='font-size:18px;'>Not enough months selected to estimate a stable trend.</div>"))
    }
    HTML(paste0(
      "<div style='font-size:22px; font-weight:900;'>β = ", sprintf("%.3f", b$beta),
      " <span style='font-size:14px; font-weight:700;'>riders/hour/month</span></div>",
      "<div style='font-size:15px;'>95% CI [", sprintf("%.3f", b$lo), ", ", sprintf("%.3f", b$hi), "]</div>",
      "<div style='font-size:13px; opacity:0.95;'>p = ", fmt_p(b$p), " · R² = ", sprintf("%.2f", b$r2), "</div>"
    ))
  })
  
  output$vb_high <- renderUI({
    df <- high_df()
    if (nrow(df) == 0) return("NA")
    last_share <- df$share_high[nrow(df)]
    if (is.na(last_share)) return("NA")
    HTML(paste0(
      "<div style='font-size:26px; font-weight:950;'>", percent(last_share, accuracy = 0.1), "</div>",
      "<div style='font-size:13px; opacity:0.95;'>Estimated share above threshold in the most recent month</div>"
    ))
  })
  
  output$p_trend <- renderPlotly({
    df <- filtered_monthly()
    if (nrow(df) == 0) return(plotly_empty())
    
    g_mean_ci <- ggplot(df, aes(x = month_date)) +
      geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = "#93c5fd", alpha = 0.35) +
      geom_line(aes(y = mean_ridership), color = "#0f172a", linewidth = 1.2) +
      geom_point(
        aes(
          y = mean_ridership,
          text = paste0(
            "<b>Month:</b> ", format(month_date, "%Y-%m"),
            "<br><b>Mean:</b> ", round(mean_ridership, 2),
            "<br><b>95% CI:</b> [", round(ci_lower, 2), ", ", round(ci_upper, 2), "]",
            "<br><b>SD:</b> ", round(sd_ridership, 2),
            "<br><b>N:</b> ", comma(n_obs)
          )
        ),
        color = "#2563eb", size = 2.3
      ) +
      scale_y_continuous(labels = comma) +
      labs(x = "Month", y = "Hourly subway entries (riders/hour)") +
      theme_minimal(base_size = 13) +
      theme(panel.grid.minor = element_blank())
    
    g_mean_sd <- ggplot(df, aes(x = month_date)) +
      geom_line(aes(y = mean_ridership), color = "#0f172a", linewidth = 1.2) +
      geom_point(
        aes(
          y = mean_ridership,
          text = paste0(
            "<b>Month:</b> ", format(month_date, "%Y-%m"),
            "<br><b>Mean:</b> ", round(mean_ridership, 2),
            "<br><b>SD:</b> ", round(sd_ridership, 2)
          )
        ),
        color = "#22c55e", size = 2.3
      ) +
      geom_line(aes(y = mean_ridership + sd_ridership), color = "#f97316", linewidth = 0.9, linetype = "dashed") +
      geom_line(aes(y = pmax(0, mean_ridership - sd_ridership)), color = "#f97316", linewidth = 0.9, linetype = "dashed") +
      scale_y_continuous(labels = comma) +
      labs(x = "Month", y = "Hourly subway entries (riders/hour)") +
      theme_minimal(base_size = 13) +
      theme(panel.grid.minor = element_blank())
    
    p <- if (input$view_mode == "mean_ci") g_mean_ci else g_mean_sd
    ggplotly(p, tooltip = "text") %>% layout(hoverlabel = list(bgcolor = "white"))
  })
  
  output$p_high <- renderPlotly({
    df <- high_df()
    if (nrow(df) == 0) return(plotly_empty())
    
    df <- df %>%
      mutate(
        text = paste0(
          "<b>Month:</b> ", format(month_date, "%Y-%m"),
          "<br><b>Threshold:</b> ", input$thr,
          "<br><b>High-demand share:</b> ", ifelse(is.na(share_high), "NA", percent(share_high, accuracy = 0.1))
        )
      )
    
    g <- ggplot(df, aes(x = month_date, y = share_high)) +
      geom_area(fill = "#fed7aa", alpha = 0.8, na.rm = TRUE) +
      geom_line(color = "#ea580c", linewidth = 1.1, na.rm = TRUE) +
      geom_point(aes(text = text), color = "#9a3412", size = 2.2, na.rm = TRUE) +
      scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
      labs(x = "Month", y = "Estimated share of high-demand station-hours") +
      theme_minimal(base_size = 13) +
      theme(panel.grid.minor = element_blank())
    
    ggplotly(g, tooltip = "text") %>% layout(hoverlabel = list(bgcolor = "white"))
  })
  
  output$txt_overview <- renderUI({
    df <- filtered_monthly()
    hd <- high_df()
    b  <- beta_obj()
    if (nrow(df) == 0) return("No data in selected window.")
    
    start_m <- format(min(df$month_date), "%Y-%m")
    end_m   <- format(max(df$month_date), "%Y-%m")
    mean_start <- df$mean_ridership[1]
    mean_end   <- df$mean_ridership[nrow(df)]
    
    hd_start <- hd$share_high[1]
    hd_end   <- hd$share_high[nrow(hd)]
    change_pp <- (hd_end - hd_start) * 100
    
    trend_sentence <- if (is.na(b$beta)) {
      "Trend β cannot be estimated reliably because the selected window is too short."
    } else {
      paste0(
        "Over ", start_m, " to ", end_m,
        ", the estimated monthly trend is β = ", sprintf("%.3f", b$beta),
        " riders/hour/month (95% CI [", sprintf("%.3f", b$lo), ", ", sprintf("%.3f", b$hi),
        "], p ", ifelse(b$p < 0.001, "< 0.001", paste0("= ", sprintf("%.3f", b$p))), ")."
      )
    }
    
    HTML(paste0(
      "<div style='font-size:15px; line-height:1.55;'>",
      "<p><b>System change:</b> Mean ridership moved from <b>", round(mean_start, 2), "</b> to <b>", round(mean_end, 2), "</b> riders/hour. ", trend_sentence, "</p>",
      "<p><b>Early-warning pressure:</b> With threshold <b>", input$thr, "</b>, estimated high-demand share changed from <b>", ifelse(is.na(hd_start), "NA", percent(hd_start, 0.1)),
      "</b> to <b>", ifelse(is.na(hd_end), "NA", percent(hd_end, 0.1)), "</b> (", sprintf("%+.2f", change_pp), " percentage points).</p>",
      "</div>"
    ))
  })
  
  high_df2 <- reactive({
    compute_high_share(filtered_monthly(), input$thr2)
  })
  
  output$p_heat <- renderPlotly({
    df <- high_df2()
    if (nrow(df) == 0) return(plotly_empty())
    
    heat <- df %>%
      mutate(
        mon_num = month(month_date),
        year = year(month_date)
      ) %>%
      group_by(year, mon_num) %>%
      summarize(share_high = mean(share_high, na.rm = TRUE), .groups = "drop") %>%
      mutate(
        mon = lubridate::month(ymd(paste0("2020-", sprintf("%02d", mon_num), "-01")), label = TRUE, abbr = TRUE),
        text = paste0(
          "<b>Year:</b> ", year,
          "<br><b>Month:</b> ", mon,
          "<br><b>High-demand share:</b> ", ifelse(is.na(share_high), "NA", percent(share_high, 0.1)),
          "<br><b>Threshold:</b> ", input$thr2
        )
      )
    
    g <- ggplot(heat, aes(x = mon, y = factor(year), fill = share_high, text = text)) +
      geom_tile(color = "white", linewidth = 0.5) +
      scale_fill_gradient(low = "#e0f2fe", high = "#ef4444", labels = percent_format(accuracy = 1), na.value = "grey85") +
      labs(x = "Month of year", y = "Year", fill = "High-demand share") +
      theme_minimal(base_size = 13) +
      theme(panel.grid = element_blank())
    
    ggplotly(g, tooltip = "text") %>% layout(hoverlabel = list(bgcolor = "white"))
  })
  
  output$p_dist <- renderPlotly({
    df <- high_df2()
    if (nrow(df) == 0) return(plotly_empty())
    
    df <- df %>%
      filter(!is.na(share_high)) %>%
      mutate(
        text = paste0(
          "<b>Month:</b> ", format(month_date, "%Y-%m"),
          "<br><b>High-demand share:</b> ", percent(share_high, 0.1)
        )
      )
    
    if (nrow(df) == 0) return(plotly_empty())
    
    g <- ggplot(df, aes(x = share_high, text = text)) +
      geom_histogram(fill = "#a78bfa", color = "white", bins = 18, alpha = 0.9) +
      geom_vline(xintercept = median(df$share_high, na.rm = TRUE), color = "#111827", linewidth = 1.0, linetype = "dashed") +
      scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
      labs(x = "High-demand share (across months)", y = "Number of months") +
      theme_minimal(base_size = 13) +
      theme(panel.grid.minor = element_blank())
    
    ggplotly(g, tooltip = "text") %>% layout(hoverlabel = list(bgcolor = "white"))
  })
  
  output$txt_limits <- renderUI({
    HTML(paste0(
      "<div style='font-size:15px; line-height:1.55;'>",
      "<p><b>Why multiple charts:</b> Different chart types help non-technical viewers understand direction (trend lines), seasonality (heatmap), and typical variability (distribution).</p>",
      "<p><b>Important limitation:</b> The high-demand share is an approximation based on monthly mean and SD because this MVP stores monthly aggregates for speed and reliability in Posit Cloud.</p>",
      "</div>"
    ))
  })
  
  output$map_borough_ui <- renderUI({
    if (is.null(stations)) {
      return(helpText("Map needs stationsmap.csv (with Longitude/Latitude). Place it in the project root or in data/."))
    }
    if (!("borough" %in% names(stations))) return(NULL)
    boroughs <- sort(unique(na.omit(stations$borough)))
    if (length(boroughs) == 0) return(NULL)
    selectInput("map_borough", "Borough filter", choices = c("All" = "ALL", boroughs), selected = "ALL")
  })
  
  map_df <- reactive({
    if (is.null(stations)) return(NULL)
    
    df <- stations
    if (!is.null(input$map_borough) && input$map_borough != "ALL" && "borough" %in% names(df)) {
      df <- df %>% filter(borough == input$map_borough)
    }
    
    eff <- station_effects %>% mutate(station_complex_id = as.character(station_complex_id))
    
    if (isTRUE(input$map_sig)) {
      eff <- eff %>% filter(sig_95)
    }
    
    if ("station_complex_id" %in% names(df)) {
      df <- df %>%
        mutate(station_complex_id = as.character(station_complex_id)) %>%
        left_join(eff, by = "station_complex_id")
    } else if ("station_complex" %in% names(df)) {
      df <- df %>% left_join(eff, by = "station_complex")
    } else {
      df <- df %>% mutate(beta = NA_real_, r2 = NA_real_, sigma = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p = NA_real_)
    }
    
    df
  })
  
  output$leaf_map <- renderLeaflet({
    base <- leaflet() %>%
      addProviderTiles(providers$CartoDB.Positron) %>%
      setView(lng = -73.94, lat = 40.70, zoom = 10)
    
    if (is.null(stations)) return(base)
    
    df <- map_df()
    if (is.null(df) || !("Longitude" %in% names(df)) || !("Latitude" %in% names(df))) return(base)
    
    metric <- input$map_metric
    
    if (!(metric %in% names(df))) {
      df$metric_value <- NA_real_
    } else {
      df$metric_value <- suppressWarnings(as.numeric(df[[metric]]))
    }
    
    pal <- colorNumeric(
      palette = c("#93c5fd", "#f97316", "#ef4444"),
      domain  = df$metric_value,
      na.color = "#9ca3af"
    )
    
    label_text <- paste0(
      if ("station_complex" %in% names(df)) df$station_complex else "Station",
      if ("borough" %in% names(df) && !all(is.na(df$borough))) paste0(" (", df$borough, ")") else "",
      "<br><b>β:</b> ", ifelse(is.na(df$beta), "NA", round(df$beta, 3)),
      "  <b>95% CI:</b> ", ifelse(is.na(df$ci_low), "NA", paste0("[", round(df$ci_low, 3), ", ", round(df$ci_high, 3), "]")),
      "<br><b>R²:</b> ", ifelse(is.na(df$r2), "NA", round(df$r2, 3)),
      "  <b>sigma:</b> ", ifelse(is.na(df$sigma), "NA", round(df$sigma, 2))
    )
    
    leaflet(df) %>%
      addProviderTiles(providers$CartoDB.Positron) %>%
      setView(lng = -73.94, lat = 40.70, zoom = 10) %>%
      addCircleMarkers(
        lng = ~Longitude,
        lat = ~Latitude,
        radius = 4,
        stroke = FALSE,
        fillOpacity = 0.85,
        color = ~pal(metric_value),
        label = lapply(label_text, HTML)
      ) %>%
      addLegend(
        position = "bottomright",
        pal = pal,
        values = df$metric_value,
        title = paste0("Map metric: ", metric),
        opacity = 1
      )
  })
  
  # -------------------------
  # Forecast page logic
  # -------------------------
  detect_group_var <- reactive({
    cand <- c("borough", "Borough", "group", "segment", "region", "line")
    cand <- cand[cand %in% names(monthly_stats)]
    if (length(cand) == 0) return(NULL)
    cand[1]
  })
  
  output$fc_group_ui <- renderUI({
    gv <- detect_group_var()
    if (is.null(gv)) {
      return(helpText("Forecast is system-wide because monthly_stats has no grouping column (e.g., borough)."))
    }
    
    lab_group <- paste0("By group (", gv, ")")
    choices_vec <- setNames(
      c("system", "group"),
      c("System-wide (all)", lab_group)
    )
    
    selectInput(
      "fc_mode",
      "Forecast level",
      choices = choices_vec,
      selected = "system"
    )
  })
  
  output$fc_group_select_ui <- renderUI({
    gv <- detect_group_var()
    if (is.null(gv)) return(NULL)
    
    req(input$fc_mode)
    if (input$fc_mode != "group") return(NULL)
    
    groups <- sort(unique(na.omit(as.character(monthly_stats[[gv]]))))
    if (length(groups) == 0) return(NULL)
    
    selectInput("fc_group", paste0("Select ", gv), choices = groups, selected = groups[1])
  })
  
  forecast_series <- reactive({
    gv <- detect_group_var()
    df <- monthly_stats %>% arrange(month_date)
    
    if (!is.null(gv) && !is.null(input$fc_mode) && input$fc_mode == "group") {
      req(input$fc_group)
      df <- df %>% filter(as.character(.data[[gv]]) == as.character(input$fc_group))
    }
    
    df
  })
  
  forecast_bundle <- reactive({
    df0 <- forecast_series()
    if (nrow(df0) < 8) return(NULL)
    req(input$policy_date)
    
    policy_date <- as.Date(input$policy_date)
    
    df_its <- make_its_frame(df0, policy_date)
    mod <- fit_its(df_its)
    
    last_date <- max(df_its$month_date, na.rm = TRUE)
    future_months <- seq(
      from = floor_date(last_date %m+% months(1), "month"),
      by = "month",
      length.out = 12
    )
    
    last_t <- max(df_its$time_index, na.rm = TRUE)
    future_base <- tibble(
      month_date = future_months,
      time_index = (last_t + 1):(last_t + 12)
    )
    
    idx0 <- which(df_its$month_date >= policy_date)[1]
    if (is.na(idx0)) idx0 <- Inf
    
    future_on <- future_base %>%
      mutate(
        post = as.integer(month_date >= policy_date),
        time_after = ifelse(is.finite(idx0), pmax(0L, time_index - idx0 + 1L), 0L)
      )
    
    future_off <- future_base %>%
      mutate(
        post = 0L,
        time_after = 0L
      )
    
    pred_on  <- predict_with_se(mod, future_on)  %>% mutate(scenario = "Policy ON")
    pred_off <- predict_with_se(mod, future_off) %>% mutate(scenario = "Policy OFF (counterfactual)")
    
    fit_in_sample <- predict_with_se(mod, df_its %>% select(month_date, time_index, post, time_after)) %>%
      mutate(scenario = "Fitted (Policy ON)")
    
    list(
      observed = df_its,
      fitted   = fit_in_sample,
      fc_on    = pred_on,
      fc_off   = pred_off,
      model    = mod
    )
  })
  
  output$vb_fc_gap <- renderUI({
    b <- forecast_bundle()
    if (is.null(b)) {
      return(HTML("<div style='font-size:18px;'>Not enough months to forecast. Select a longer range or ensure the dataset has at least 8 monthly points.</div>"))
    }
    on <- b$fc_on
    off <- b$fc_off
    gap <- mean(on$fit - off$fit, na.rm = TRUE)
    
    HTML(paste0(
      "<div style='font-size:26px; font-weight:950;'>",
      ifelse(is.na(gap), "NA", comma(round(gap, 2))), "</div>",
      "<div style='font-size:13px; opacity:0.95;'>Average Policy ON − Policy OFF over next 12 months</div>"
    ))
  })
  
  output$p_forecast <- renderPlotly({
    b <- forecast_bundle()
    if (is.null(b)) return(plotly_empty())
    
    obs <- b$observed
    fit <- b$fitted
    on  <- b$fc_on
    off <- b$fc_off
    
    obs <- obs %>%
      mutate(text = paste0(
        "<b>Month:</b> ", format(month_date, "%Y-%m"),
        "<br><b>Observed mean:</b> ", round(mean_ridership, 2)
      ))
    
    fit <- fit %>%
      mutate(text = paste0(
        "<b>Month:</b> ", format(month_date, "%Y-%m"),
        "<br><b>Fitted (Policy ON):</b> ", round(fit, 2),
        "<br><b>95%:</b> [", round(lo, 2), ", ", round(hi, 2), "]"
      ))
    
    on <- on %>%
      mutate(text = paste0(
        "<b>Month:</b> ", format(month_date, "%Y-%m"),
        "<br><b>Forecast Policy ON:</b> ", round(fit, 2),
        "<br><b>95%:</b> [", round(lo, 2), ", ", round(hi, 2), "]"
      ))
    
    off <- off %>%
      mutate(text = paste0(
        "<b>Month:</b> ", format(month_date, "%Y-%m"),
        "<br><b>Forecast Policy OFF:</b> ", round(fit, 2),
        "<br><b>95%:</b> [", round(lo, 2), ", ", round(hi, 2), "]"
      ))
    
    g <- ggplot() +
      geom_point(
        data = obs,
        aes(x = month_date, y = mean_ridership, text = text),
        color = "#0f172a", size = 2.2, alpha = 0.9
      ) +
      geom_ribbon(
        data = fit,
        aes(x = month_date, ymin = lo, ymax = hi),
        fill = "#93c5fd", alpha = 0.22
      ) +
      geom_line(
        data = fit,
        aes(x = month_date, y = fit),
        color = "#2563eb", linewidth = 1.0
      ) +
      geom_ribbon(
        data = off,
        aes(x = month_date, ymin = lo, ymax = hi),
        fill = "grey70", alpha = 0.20
      ) +
      geom_line(
        data = off,
        aes(x = month_date, y = fit),
        color = "grey40", linewidth = 1.1
      ) +
      geom_ribbon(
        data = on,
        aes(x = month_date, ymin = lo, ymax = hi),
        fill = "#93c5fd", alpha = 0.18
      ) +
      geom_line(
        data = on,
        aes(x = month_date, y = fit),
        color = "#2563eb", linewidth = 1.4
      ) +
      scale_y_continuous(labels = comma) +
      labs(x = "Month", y = "Hourly subway entries (riders/hour)") +
      theme_minimal(base_size = 13) +
      theme(panel.grid.minor = element_blank())
    
    ggplotly(g, tooltip = "text") %>% layout(hoverlabel = list(bgcolor = "white"))
  })
  
  # -------------------------
  # FIXED GAP PLOTS (no vjust referencing columns outside aes)
  # -------------------------
  
  output$p_gap_monthly <- renderPlotly({
    b <- forecast_bundle()
    if (is.null(b)) return(plotly_empty())
    
    on  <- b$fc_on
    off <- b$fc_off
    
    gap <- on %>%
      select(month_date, fit_on = fit, lo_on = lo, hi_on = hi) %>%
      left_join(off %>% select(month_date, fit_off = fit, lo_off = lo, hi_off = hi), by = "month_date") %>%
      mutate(
        gap = fit_on - fit_off,
        se_on  = (hi_on  - fit_on) / 1.96,
        se_off = (hi_off - fit_off) / 1.96,
        se_gap = sqrt(pmax(0, se_on^2 + se_off^2)),
        lo_gap = gap - 1.96 * se_gap,
        hi_gap = gap + 1.96 * se_gap,
        label = round(gap, 1),
        text = paste0(
          "<b>Month:</b> ", format(month_date, "%Y-%m"),
          "<br><b>Gap (ON − OFF):</b> ", comma(round(gap, 2)),
          "<br><b>95% (approx):</b> [", comma(round(lo_gap, 2)), ", ", comma(round(hi_gap, 2)), "]"
        )
      ) %>%
      filter(!is.na(gap))
    
    if (nrow(gap) == 0) return(plotly_empty())
    
    pad <- 0.15 * max(1, diff(range(gap$gap, na.rm = TRUE)))
    y_min <- min(gap$lo_gap, na.rm = TRUE) - pad
    y_max <- max(gap$hi_gap, na.rm = TRUE) + pad
    
    gap <- gap %>%
      mutate(
        label_y = ifelse(gap >= 0, hi_gap + 0.06 * (y_max - y_min), lo_gap - 0.06 * (y_max - y_min))
      )
    
    g <- ggplot(gap, aes(x = month_date, y = gap, text = text)) +
      geom_hline(yintercept = 0, linewidth = 0.8, color = "#111827", alpha = 0.50) +
      geom_col(fill = "#22c55e", alpha = 0.80, width = 25) +
      geom_errorbar(aes(ymin = lo_gap, ymax = hi_gap), width = 8, alpha = 0.60) +
      geom_text(aes(y = label_y, label = label), size = 3.6) +
      scale_y_continuous(labels = comma, limits = c(y_min, y_max)) +
      scale_x_date(date_breaks = "3 months", date_labels = "%b %Y") +
      labs(x = "Month", y = "Projected gap (riders/hour)") +
      theme_minimal(base_size = 13) +
      theme(panel.grid.minor = element_blank())
    
    ggplotly(g, tooltip = "text") %>% layout(hoverlabel = list(bgcolor = "white"))
  })
  
  output$p_gap_cum <- renderPlotly({
    b <- forecast_bundle()
    if (is.null(b)) return(plotly_empty())
    
    on  <- b$fc_on
    off <- b$fc_off
    
    gap <- on %>%
      select(month_date, fit_on = fit) %>%
      left_join(off %>% select(month_date, fit_off = fit), by = "month_date") %>%
      mutate(
        gap = fit_on - fit_off,
        gap0 = ifelse(is.na(gap), 0, gap),
        gap_cum = cumsum(gap0),
        label = round(gap_cum, 0),
        text = paste0(
          "<b>Month:</b> ", format(month_date, "%Y-%m"),
          "<br><b>Monthly gap:</b> ", comma(round(gap, 2)),
          "<br><b>Cumulative:</b> ", comma(round(gap_cum, 2))
        )
      ) %>%
      filter(!is.na(gap_cum))
    
    if (nrow(gap) == 0) return(plotly_empty())
    
    yr <- range(gap$gap_cum, na.rm = TRUE)
    span <- max(1, diff(yr))
    
    # Push labels upward so they don't sit on the points
    offset <- 0.10 * span
    gap <- gap %>%
      mutate(label_y = gap_cum + offset)
    
    g <- ggplot(gap, aes(x = month_date, y = gap_cum, text = text)) +
      geom_hline(yintercept = 0, linewidth = 0.8, color = "#111827", alpha = 0.50) +
      geom_line(color = "#2563eb", linewidth = 1.3) +
      geom_point(color = "#2563eb", size = 2.4, alpha = 0.95) +
      geom_text(aes(y = label_y, label = label), size = 3.6) +
      scale_y_continuous(
        labels = comma,
        expand = expansion(mult = c(0.06, 0.20))
      ) +
      scale_x_date(date_breaks = "3 months", date_labels = "%b %Y") +
      labs(x = "Month", y = "Cumulative gap (riders/hour·month)") +
      theme_minimal(base_size = 13) +
      theme(panel.grid.minor = element_blank())
    
    ggplotly(g, tooltip = "text") %>% layout(hoverlabel = list(bgcolor = "white"))
  })
  
  output$txt_gap_explain <- renderUI({
    HTML(paste0(
      "<div style='font-size:14px; line-height:1.55; padding-top:8px;'>",
      "<p><b>Interpretation for this dashboard:</b> In the current view, the Policy ON and Policy OFF forecasts appear very similar. ",
      "This reflects the model output under the selected time window and policy start date: the estimated post-policy level shift (<i>post</i>) and/or slope change (<i>time_after</i>) is small. ",
      "Therefore, the projected 12-month gap (Policy ON − Policy OFF) remains modest and the cumulative gap evolves smoothly.</p>",
      "<p><b>What would create clearer separation:</b> A larger visual gap should come from stronger signal in the underlying series rather than from the plotting style. ",
      "In practice, this often involves adjusting the policy start date, forecasting a more policy-sensitive subgroup (e.g., Manhattan vs. non-Manhattan), ",
      "or using an outcome that more directly reflects the policy mechanism. ",
      "With the current modeling setup, the <b>cumulative gap</b> is the most transparent summary of small month-to-month differences.</p>",
      "</div>"
    ))
  })
  
  output$txt_forecast <- renderUI({
    b <- forecast_bundle()
    if (is.null(b)) {
      return(HTML("<div style='font-size:15px; line-height:1.55;'><p><b>Forecast unavailable:</b> The selected series does not have enough months to support an interrupted time-series model. Choose a longer range or verify the dataset contains at least 8 monthly points.</p></div>"))
    }
    
    on <- b$fc_on
    off <- b$fc_off
    
    last_obs <- max(b$observed$month_date, na.rm = TRUE)
    end_fc   <- max(on$month_date, na.rm = TRUE)
    
    gap <- mean(on$fit - off$fit, na.rm = TRUE)
    gap_txt <- ifelse(is.na(gap), "NA", paste0(comma(round(gap, 2)), " riders/hour"))
    
    HTML(paste0(
      "<div style='font-size:15px; line-height:1.55;'>",
      "<p><b>What this forecast means:</b> We extend the historical trend forward by 12 months, comparing two scenarios. ",
      "<b>Policy ON</b> follows the post-policy level/slope structure estimated from the data, while <b>Policy OFF</b> continues the pre-policy trend only.</p>",
      "<p><b>Time window:</b> Forecast runs from <b>", format(last_obs %m+% months(1), "%Y-%m"),
      "</b> to <b>", format(end_fc, "%Y-%m"), "</b>. Across these 12 months, the average projected gap (Policy ON − Policy OFF) is <b>", gap_txt, "</b>.</p>",
      "</div>"
    ))
  })
  
  # -------------------------
  # Spatial Analysis page logic
  # -------------------------
  
  sp_filtered <- reactive({
    if (is.null(stations_sf)) return(NULL)
    df <- stations_sf
    if (!is.null(input$sp_borough) && input$sp_borough != "ALL" && "borough" %in% names(df)) {
      df <- df %>% filter(borough == input$sp_borough)
    }
    df
  })
  
  sp_neighbors <- reactive({
    df <- sp_filtered()
    if (is.null(df) || nrow(df) < 2) return(NULL)
    req(input$sp_station)
    
    sel <- df %>% filter(station_complex == input$sp_station)
    if (nrow(sel) == 0) return(NULL)
    
    radius_m <- input$sp_radius_km * 1000
    buf <- st_buffer(sel[1, ], dist = radius_m)
    within <- st_within(df, buf, sparse = FALSE)[, 1]
    neighbors <- df[within, ]
    
    list(
      selected = sel[1, ],
      buffer = buf,
      neighbors = neighbors
    )
  })
  
  output$vb_nn_dist <- renderUI({
    df <- sp_filtered()
    if (is.null(df) || nrow(df) == 0) {
      return(HTML("<div style='font-size:18px;'>No spatial data available.</div>"))
    }
    
    mean_nn <- mean(df$nn_dist_km, na.rm = TRUE)
    med_nn <- median(df$nn_dist_km, na.rm = TRUE)
    
    HTML(paste0(
      "<div style='font-size:22px; font-weight:900;'>", sprintf("%.2f", mean_nn),
      " <span style='font-size:14px; font-weight:700;'>km (mean)</span></div>",
      "<div style='font-size:15px;'>Median: ", sprintf("%.2f", med_nn), " km</div>",
      "<div style='font-size:13px; opacity:0.95;'>Computed via sf::st_distance across ", nrow(df), " stations</div>"
    ))
  })
  
  output$vb_spatial_beta <- renderUI({
    nb <- sp_neighbors()
    if (is.null(nb)) {
      return(HTML("<div style='font-size:18px;'>Select a station to compute spatial average.</div>"))
    }
    
    n <- nb$neighbors
    if (!"beta" %in% names(n)) {
      return(HTML("<div style='font-size:18px;'>No β data for these stations.</div>"))
    }
    
    avg_beta <- mean(n$beta, na.rm = TRUE)
    n_in <- nrow(n)
    
    HTML(paste0(
      "<div style='font-size:22px; font-weight:900;'>", sprintf("%.3f", avg_beta),
      " <span style='font-size:14px; font-weight:700;'>riders/hour/month</span></div>",
      "<div style='font-size:15px;'>Spatial avg across ", n_in, " stations within ", input$sp_radius_km, " km</div>",
      "<div style='font-size:13px; opacity:0.95;'>Service area computed via sf::st_buffer + sf::st_within</div>"
    ))
  })
  
  output$sp_map <- renderLeaflet({
    df <- sp_filtered()
    if (is.null(df) || nrow(df) == 0) {
      return(leaflet() %>% addProviderTiles(providers$CartoDB.Positron) %>% setView(-73.94, 40.70, 10))
    }
    
    df_wgs <- st_transform(df, 4326)
    coords <- st_coordinates(df_wgs)
    df_plot <- df_wgs %>% mutate(lng = coords[, 1], lat = coords[, 2])
    
    pal <- colorNumeric(
      palette = c("#22c55e", "#eab308", "#ef4444"),
      domain = df_plot$nn_dist_km,
      na.color = "#9ca3af"
    )
    
    label_text <- paste0(
      "<b>", df_plot$station_complex, "</b>",
      if ("borough" %in% names(df_plot)) paste0(" (", df_plot$borough, ")") else "",
      "<br>NN distance: ", round(df_plot$nn_dist_km, 2), " km",
      "<br>Nearest: ", df_plot$nn_name,
      "<br>β: ", ifelse(is.na(df_plot$beta), "NA", round(df_plot$beta, 3))
    )
    
    m <- leaflet(df_plot) %>%
      addProviderTiles(providers$CartoDB.Positron) %>%
      addCircleMarkers(
        lng = ~lng, lat = ~lat,
        radius = 5, stroke = FALSE, fillOpacity = 0.85,
        color = ~pal(nn_dist_km),
        label = lapply(label_text, HTML)
      ) %>%
      addLegend(
        position = "bottomright", pal = pal,
        values = df_plot$nn_dist_km,
        title = "NN distance (km)", opacity = 1
      )
    
    # Add buffer circle if station selected
    nb <- sp_neighbors()
    if (!is.null(nb)) {
      buf_wgs <- st_transform(nb$buffer, 4326)
      sel_wgs <- st_transform(nb$selected, 4326)
      sel_coords <- st_coordinates(sel_wgs)
      
      m <- m %>%
        addPolygons(data = buf_wgs, fillColor = "#6366f1", fillOpacity = 0.12,
                    color = "#6366f1", weight = 2, dashArray = "5,5") %>%
        addCircleMarkers(lng = sel_coords[1, 1], lat = sel_coords[1, 2],
                         radius = 9, color = "#6366f1", fillColor = "#6366f1",
                         fillOpacity = 1, stroke = TRUE, weight = 2,
                         label = HTML(paste0("<b>Selected: ", nb$selected$station_complex, "</b>")))
    }
    
    m
  })
  
  output$sp_scatter <- renderPlotly({
    df <- sp_filtered()
    if (is.null(df) || nrow(df) == 0 || !"beta" %in% names(df)) return(plotly_empty())
    
    df_plot <- df %>%
      filter(!is.na(beta), !is.na(nn_dist_km)) %>%
      mutate(
        text = paste0(
          "<b>", station_complex, "</b>",
          "<br>NN dist: ", round(nn_dist_km, 2), " km",
          "<br>β: ", round(beta, 3),
          if ("borough" %in% names(df)) paste0("<br>Borough: ", borough) else ""
        )
      )
    
    if (nrow(df_plot) == 0) return(plotly_empty())
    
    g <- ggplot(df_plot, aes(x = nn_dist_km, y = beta, text = text)) +
      geom_point(aes(color = borough),
                 size = 2.5, alpha = 0.8) +
      geom_smooth(method = "lm", se = TRUE, color = "#6366f1", linewidth = 1, alpha = 0.15) +
      labs(x = "Nearest-neighbor distance (km)", y = "Station trend β (riders/hour/month)",
           color = "Borough") +
      theme_minimal(base_size = 13) +
      theme(panel.grid.minor = element_blank())
    
    ggplotly(g, tooltip = "text") %>% layout(hoverlabel = list(bgcolor = "white"))
  })
  
  output$sp_txt_neighborhood <- renderUI({
    nb <- sp_neighbors()
    if (is.null(nb)) {
      return(HTML("<div style='font-size:15px;'>Select a station to see its spatial neighborhood summary.</div>"))
    }
    
    sel <- nb$selected
    n <- nb$neighbors
    n_count <- nrow(n)
    avg_beta <- mean(n$beta, na.rm = TRUE)
    sel_beta <- if (!is.na(sel$beta)) round(sel$beta, 3) else "NA"
    sel_nn <- round(sel$nn_dist_km, 2)
    
    boroughs_in <- if ("borough" %in% names(n)) paste(sort(unique(na.omit(n$borough))), collapse = ", ") else "N/A"
    max_beta_station <- if (any(!is.na(n$beta))) n$station_complex[which.max(n$beta)] else "N/A"
    min_beta_station <- if (any(!is.na(n$beta))) n$station_complex[which.min(n$beta)] else "N/A"
    
    # Compute area of the buffer in km²
    area_km2 <- as.numeric(st_area(nb$buffer)) / 1e6
    density <- n_count / area_km2
    
    HTML(paste0(
      "<div style='font-size:15px; line-height:1.55;'>",
      "<p><b>Selected station:</b> ", sel$station_complex,
      " — own β = ", sel_beta,
      ", nearest neighbor at ", sel_nn, " km (", sel$nn_name, ").</p>",
      "<p><b>Service area (", input$sp_radius_km, " km radius):</b> ",
      n_count, " stations within a ", round(area_km2, 2), " km² buffer ",
      "(station density: ", round(density, 1), " per km²). ",
      "Boroughs represented: ", boroughs_in, ".</p>",
      "<p><b>Spatial average β:</b> ", round(avg_beta, 3), " riders/hour/month across neighbors. ",
      "Strongest growth: <b>", max_beta_station, "</b>; ",
      "weakest: <b>", min_beta_station, "</b>.</p>",
      "</div>"
    ))
  })
  
  output$sp_txt_clustering <- renderUI({
    df <- sp_filtered()
    if (is.null(df) || nrow(df) == 0) {
      return(HTML("<div style='font-size:15px;'>No spatial data available.</div>"))
    }
    
    mean_nn <- mean(df$nn_dist_km, na.rm = TRUE)
    sd_nn <- sd(df$nn_dist_km, na.rm = TRUE)
    
    # Identify most isolated and most clustered stations
    most_isolated <- df$station_complex[which.max(df$nn_dist_km)]
    max_dist <- round(max(df$nn_dist_km, na.rm = TRUE), 2)
    most_clustered <- df$station_complex[which.min(df$nn_dist_km)]
    min_dist <- round(min(df$nn_dist_km, na.rm = TRUE), 2)
    
    # Correlation between nn_dist and beta
    cor_val <- if ("beta" %in% names(df)) {
      ct <- cor.test(df$nn_dist_km, df$beta, use = "pairwise.complete.obs")
      list(r = ct$estimate, p = ct$p.value)
    } else {
      list(r = NA, p = NA)
    }
    
    HTML(paste0(
      "<div style='font-size:15px; line-height:1.55;'>",
      "<p><b>Station spacing:</b> Across the ", nrow(df), " stations shown, ",
      "average nearest-neighbor distance is <b>", round(mean_nn, 2), " km</b> (SD = ", round(sd_nn, 2), " km). ",
      "Most isolated: <b>", most_isolated, "</b> (", max_dist, " km to nearest). ",
      "Most clustered: <b>", most_clustered, "</b> (", min_dist, " km to nearest).</p>",
      if (!is.na(cor_val$r)) paste0(
        "<p><b>Spatial-trend relationship:</b> Correlation between station isolation (NN distance) and ridership trend β is ",
        "r = ", sprintf("%.3f", cor_val$r), " (p ", ifelse(cor_val$p < 0.001, "< 0.001", paste0("= ", sprintf("%.3f", cor_val$p))),
        "). ",
        ifelse(abs(cor_val$r) < 0.1, "There is essentially no linear association between station spacing and growth trend.",
               ifelse(cor_val$r > 0, "More isolated stations tend to show slightly stronger ridership growth.",
                      "More clustered (central) stations tend to show slightly stronger ridership growth.")),
        "</p>"
      ) else "",
      "</div>"
    ))
  })
}

shinyApp(ui = ui, server = server)