**NYC Subway Early Warning Dashboard — Congestion Pricing (MVP)**

**About This Dashboard:**
This dashboard presents an early-warning decision-support system for understanding changes in NYC subway ridership under congestion pricing conditions. It integrates system-wide trends, station-level spatial patterns, and short-term forecasts to help users identify where demand pressure is increasing and how changes are distributed across the subway network.
Temporal indicators on the Overview and Forecast pages summarize long-term trends and projected differences between Policy ON and Policy OFF scenarios over the next 12 months. The Map and Details pages provide station-level and seasonal views that help interpret how demand variability evolves across locations and time.
The Spatial Analysis page introduces sf-based spatial calculations, including nearest-neighbor station distance, service-area buffers, and spatial averages of local ridership trends. These spatial metrics help evaluate station spacing, accessibility, and neighborhood-level clustering patterns, supporting interpretation of how local network structure relates to observed ridership change across the system.

**Data Sources:**
The dashboard uses New York State Open Data’s MTA Subway Hourly Ridership datasets and aggregates them to monthly summaries for stability and performance.
The interactive map uses a station coordinate reference file named stationsmap.csv (located either in the project root or in data/).

**How to run:**
1. Open this project in RStudio: Shiny.R
2. Ensure the /data folder contains:
   - monthly_stats.rds
   - station_effects.rds
   - stationsmap.csv
   - my_area.geojson
3. Run:
   shiny::runApp()
