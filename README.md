# cmip6-climate-scenarios-SAWS-DFFE
This repository provides tools and scripts to calculate climate indices using CMIP6 climate projections. It utilises data from NASA’s Earth Exchange Global Daily Downscaled Projections (NEX-GDDP-CMIP6) to analyse climate extremes across South Africa. The analysis encompasses multiple Shared Socioeconomic Pathways (SSP) climate change scenarios spanning various future timeframes for South Africa. It also integrates information from bio-climatic regions and urban centres to support regional climate impact assessments.

## List of supported models 
Configure to download data from a subset of *25 CMIP6 models* selected for data consistency and availability.

    ACCESS-CM2
    ACCESS-ESM1-5
    BCC-CSM2-MR
    CMCC-CM2-SR5
    CMCC-ESM2
    CanESM5
    EC-Earth3
    EC-Earth3-Veg-LR
    GFDL-CM4
    GFDL-CM4_gr2
    GFDL-ESM4
    IITM-ESM
    INM-CM4-8
    INM-CM5-0
    IPSL-CM6A-LR
    KACE-1-0-G
    KIOST-ESM
    MIROC6
    MPI-ESM1-2-HR
    MPI-ESM1-2-LR
    MRI-ESM2-0
    NESM3
    NorESM2-LM
    NorESM2-MM
    TaiESM1

## Ensemble member
All model downloads use the same ensemble member
    *r1i1p1f1*   

## Grid label consideration
Please note the grid label varies between models, therefor adjust accordingly for data download.
    "gn"
    "gr" 
    "gr1" and "gr2" 

## Getting started
      git clone https://github.com/MfopaC cmip6-climate-scenarios-SAWS-DFFE.git
      cd cmip6-climate-scenarios-SAWS-DFFE
      cd scripts
      cd pr

## Install libraries
    pip install xarray pandas geopandas regionmask matplotlib numpy cftime 

## Run script

