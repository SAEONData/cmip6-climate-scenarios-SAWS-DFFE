# Climate Indices for Rainfall

                    | Short Name  | Description                                              | Units    |
                    | ----------- | -------------------------------------------------------- | -------- |
                    | **CDD**     | Max consecutive dry days (PR < 1 mm/day)                 | days     |
                    | **CWD**     | Max consecutive wet days (PR ≥ 1 mm/day)                 | days     |
                    | **R10mm**   | Days with ≥ 10 mm of rainfall                            | days     |
                    | **R20mm**   | Days with ≥ 20 mm of rainfall                            | days     |
                    | **PRCPTOT** | Annual rainfall on wet days (PR ≥ 1 mm/day)              | mm       |
                    | **R95pTOT** | % of rainfall from very wet days (>95th percentile)      | %        |
                    | **R99pTOT** | % of rainfall from extremely wet days (>99th percentile) | %        |
                    | **R95p**    | Total rainfall from very wet days (>95th percentile)     | mm       |
                    | **R99p**    | Total rainfall from extreme days (>99th percentile)      | mm       |
                    | **SPI**     | Standardised Precipitation Index (3, 6, 12 months)       | unitless |
                    | **SDII**    | Rainfall intensity on wet days                           | mm/day   |
                    | **Rx1day**  | Wettest day of the year                                  | mm       |
                    | **Rx5day**  | Wettest 5-day period of the year                         | mm       |
                    

  # Climate Indices for Temperature
                    | Short Name  | Description                                                       | Units    |
                    | ----------- | ----------------------------------------------------------------- | -------- |
                    | **FD**      | Days with minimum temperature below 0°C (frost days)              | days     |
                    | **TNlt2**   | Days with minimum temperature below 2°C                           | days     |
                    | **TXx**     | Warmest daily maximum temperature                                 | °C       |
                    | **TNn**     | Coldest daily minimum temperature                                 | °C       |
                    | **WSDI**    | Warm spell duration: TX > 90th percentile for ≥ 6 consecutive days| days     |
                    | **CSDI**    | Cold spell duration: TN < 10th percentile for ≥ 6 consecutive days| days     |
                    | **TXgt50p** | % of days with TX above the 50th percentile                       | %        |
                    | **TXge30**  | Days with TX ≥ 30°C                                               | days     |
                    | **TXdTNd**  | Consecutive days where both TX & TN > 95th percentile             | events   |
                    | **TNx**     | Warmest daily minimum temperature (hottest night)                 | °C       |
                    | **TXn**     | Coldest daily maximum temperature (coldest day)                   | °C       |
                    | **TX10p**   | % of days with TX < 10th percentile (cool days)                   | %        |
                    | **TX90p**   | % of days with TX > 90th percentile (hot days)                    | %        |
                    | **TN10p**   | % of days with TN < 10th percentile (cold nights)                 | %        |
                    | **TN90p**   | % of days with TN > 90th percentile (warm nights)                 | %        |


# CMIP Variable Descriptions
            | Variable  | Full Name                               | Unit       | Description                                                                 |
            | --------- | ---------------------------------------- | ---------- | --------------------------------------------------------------------------- |
            | **hurs**  | Relative Humidity                        | %          | Near-surface relative humidity (usually at 2 meters above ground).          |
            | **huss**  | Specific Humidity                        | kg/kg      | Near-surface specific humidity, the ratio of water vapor mass to air mass.  |
            | **pr**    | Precipitation Rate                       | kg/m²/s    | Total precipitation rate. Multiply by 86400 to get mm/day.                 |
            | **rssds** | Surface Downwelling Shortwave Radiation  | W/m²       | Solar radiation received at the surface (clear-sky and all-sky).            |
            | **sfcWind**| Near-surface Wind Speed                 | m/s        | Magnitude of wind speed at 10 meters above the surface.                     |
            | **tas**   | Near-surface Air Temperature             | K          | Daily average temperature at 2 meters above ground.                         |
            | **tasmax**| Daily Maximum Near-surface Temperature   | K          | Maximum daily temperature at 2 meters above the surface.                    |
            | **tasmin**| Daily Minimum Near-surface Temperature   | K          | Minimum daily temperature at 2 meters above the surface.    
