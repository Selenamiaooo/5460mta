library(httr)
library(jsonlite)
library(dplyr)
library(lubridate)
library(broom)

dir.create("data", showWarnings = FALSE)

base_2020_2024 <- "https://data.ny.gov/resource/wujg-7c2s.json"
base_2025_on   <- "https://data.ny.gov/resource/5wq4-mkjj.json"

soql_get <- function(base, query_list) {
  url <- modify_url(base, query = query_list)
  res <- GET(url)
  if (http_error(res)) {
    message("Request failed: ", url)
    message(content(res, as = "text", encoding = "UTF-8"))
    stop_for_status(res)
  }
  fromJSON(content(res, as = "text", encoding = "UTF-8"))
}

start_date <- "2021-01-01T00:00:00.000"
end_date   <- "2025-12-31T23:59:59.999"
mode_keep  <- "subway"

where_clause <- paste0(
  "transit_timestamp >= '", start_date, "' ",
  "AND transit_timestamp <= '", end_date, "' ",
  "AND transit_mode = '", mode_keep, "'"
)

get_monthly_stats <- function(base) {
  soql_get(
    base,
    list(
      `$select` = paste0(
        "date_trunc_ym(transit_timestamp) as month, ",
        "avg(ridership) as mean_ridership, ",
        "stddev_samp(ridership) as sd_ridership, ",
        "count(*) as n_obs"
      ),
      `$where` = where_clause,
      `$group` = "month",
      `$order` = "month"
    )
  )
}

monthly_stats <- bind_rows(
  get_monthly_stats(base_2020_2024),
  get_monthly_stats(base_2025_on)
) %>%
  mutate(
    mean_ridership = as.numeric(mean_ridership),
    sd_ridership   = as.numeric(sd_ridership),
    n_obs          = as.numeric(n_obs),
    month_date     = as.Date(substr(month, 1, 10))
  ) %>%
  arrange(month_date) %>%
  mutate(
    se       = sd_ridership / sqrt(n_obs),
    ci_lower = mean_ridership - 1.96 * se,
    ci_upper = mean_ridership + 1.96 * se,
    time_index = (year(month_date) - year(min(month_date))) * 12 +
      (month(month_date) - month(min(month_date))),
    year = year(month_date),
    mon  = month(month_date, label = TRUE, abbr = TRUE)
  )

saveRDS(monthly_stats, "data/monthly_stats.rds")
message("Saved: data/monthly_stats.rds (rows = ", nrow(monthly_stats), ")")

get_station_month <- function(base) {
  soql_get(
    base,
    list(
      `$select` = paste0(
        "date_trunc_ym(transit_timestamp) as month, ",
        "station_complex_id, ",
        "station_complex, ",
        "avg(ridership) as mean_ridership"
      ),
      `$where` = where_clause,
      `$group` = "month, station_complex_id, station_complex",
      `$order` = "station_complex_id, month",
      `$limit` = "50000"
    )
  )
}

station_month <- bind_rows(
  get_station_month(base_2020_2024),
  get_station_month(base_2025_on)
) %>%
  mutate(
    station_complex_id = as.character(station_complex_id),
    station_complex    = as.character(station_complex),
    mean_ridership     = as.numeric(mean_ridership),
    month_date         = as.Date(substr(month, 1, 10))
  ) %>%
  arrange(station_complex_id, month_date) %>%
  mutate(
    time_index = (year(month_date) - year(min(month_date))) * 12 +
      (month(month_date) - month(min(month_date)))
  )

station_effects <- station_month %>%
  group_by(station_complex_id, station_complex) %>%
  filter(sum(!is.na(mean_ridership)) >= 24) %>%
  group_modify(~{
    m <- lm(mean_ridership ~ time_index, data = .x)
    tr <- tidy(m, conf.int = TRUE) %>% filter(term == "time_index")
    gl <- glance(m)
    tibble(
      n_months = nrow(.x),
      beta     = tr$estimate,
      se       = tr$std.error,
      p        = tr$p.value,
      ci_low   = tr$conf.low,
      ci_high  = tr$conf.high,
      r2       = gl$r.squared,
      sigma    = gl$sigma
    )
  }) %>%
  ungroup() %>%
  mutate(sig_95 = (ci_low > 0) | (ci_high < 0))

saveRDS(station_effects, "data/station_effects.rds")
message("Saved: data/station_effects.rds (stations = ", nrow(station_effects), ")")

message("Scope: ", format(min(monthly_stats$month_date)), " to ", format(max(monthly_stats$month_date)))