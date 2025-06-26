import os
import glob
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# ------------------ Configuration ------------------ #
data_path = "/home/caroline/nex-sa-tools/data/pr" #  Replace with your actual path
lat_bounds = [-35, -22]
lon_bounds = [16, 33]
threshold = 1.0  # mm/day (wet day threshold)

# ------------------ Find Historical Files ------------------ #
all_hist_files = sorted(glob.glob(os.path.join(data_path, "**", "historical", "*.nc"), recursive=True))
ensemble_prcptot = []

# ------------------ Process Each File ------------------ #
for file in all_hist_files:
    try:
        ds = xr.open_dataset(file)
        pr = ds['pr'] * 86400  # Convert from kg/m²/s to mm/day

        # Subset to region
        pr = pr.sel(lat=slice(*lat_bounds), lon=slice(*lon_bounds))

        # Optional: restrict to specific time period
        # pr = pr.sel(time=slice("1981-01-01", "2010-12-31"))

        # Mask out dry days (< 1.0 mm/day)
        wet_pr = pr.where(pr >= threshold)

        # Compute annual wet-day total precipitation
        annual_wet_total = wet_pr.resample(time='YE').sum(dim='time')

        # Compute long-term mean over all years
        prcptot = annual_wet_total.mean(dim='time')

        ensemble_prcptot.append(prcptot)
        print(f"✅ Processed: {os.path.basename(file)}")

    except Exception as e:
        print(f"❌ Error in {file}: {e}")

if ensemble_prcptot:
    prcptot_stack = xr.concat(ensemble_prcptot, dim="model")
    ensemble_mean_prcptot = prcptot_stack.mean(dim="model")
    ensemble_mean_prcptot.name = "PRCPTOT"
    ensemble_mean_prcptot.to_netcdf("ensemble_prcptot_mean.nc")
    print("✅ Saved to 'ensemble_prcptot_mean.nc'")

    # ------------------ Plot ------------------ #
    fig = plt.figure(figsize=(10, 6))
    ax = plt.axes(projection=ccrs.PlateCarree())

    ensemble_mean_prcptot.plot(
        ax=ax,
        transform=ccrs.PlateCarree(),
        cmap='GnBu',
        cbar_kwargs={'label': 'Total Wet-Day Precipitation (mm/year)'},
    )

    ax.coastlines()
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.set_title("CMIP6 PRCPTOT (Historical)", fontsize=14)
    ax.set_extent([16, 33, -35, -22])
    plt.tight_layout()
    plt.show()
