# cmip6-climate-scenarios-SAWS-DFFE
This repository provides tools and scripts for calculating climate indices developed by the South African Weather Service using CMIP6 climate projections across South Africa. The analysis combines data from bio-climatic regions and urban centers to asses extreme climate indices.

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

## Getting started
git clone https://github.com/MfopaC cmip6-climate-scenarios-SAWS-DFFE.git

cd cmip6-climate-scenarios-SAWS-DFFE

cd scripts

cd pr

## Install libraries
pip install xarray pandas geopandas regionmask matplotlib numpy cftime 

## Run script

