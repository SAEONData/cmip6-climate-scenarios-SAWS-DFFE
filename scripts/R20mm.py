import os
import glob
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# ------------------ Configuration ------------------ #
data_path = "/home/caroline/nex-gddp-sa-tools/pr"  #  Replace with your actual CMIP6 root path
lat_bounds = [-35, -22]
lon_bounds = [16, 33]
threshold = 20.0  # mm/day for R20mm

# ------------------ Find Only Historical Files ------------------ #
all_hist_files = sorted(glob.glob(os.path.join(data_path, "**", "historical", "*.nc"), recursive=True))
ensemble_r20 = []

# ------------------ Process Each File ------------------ #
for file in all_hist_files:
    try:
        ds = xr.open_dataset(file)

        # Convert PR from kg/m²/s to mm/day
        pr = ds['pr'] * 86400

        # Subset to region of interest
        pr = pr.sel(lat=slice(*lat_bounds), lon=slice(*lon_bounds))

        # (Optional) Time slice e.g., 1981–2010
        # pr = pr.sel(time=slice("1981-01-01", "2010-12-31"))

        # Count days where PR ≥ 20 mm
        r20 = (pr >= threshold).sum(dim="time")

        ensemble_r20.append(r20)
        print(f"✅ Processed: {os.path.basename(file)}")

    except Exception as e:
        print(f"❌ Error in {file}: {e}")

# ------------------ Ensemble Mean ------------------ #
if ensemble_r20:
    r20_stack = xr.concat(ensemble_r20, dim="model")
    ensemble_mean_r20 = r20_stack.mean(dim="model")
    ensemble_mean_r20.name = "R20mm"
    ensemble_mean_r20.to_netcdf("ensemble_r20mm_mean.nc")
    print("✅ Saved to 'ensemble_r20mm_mean.nc'")

 # ------------------ Plot ------------------ #
    fig = plt.figure(figsize=(10, 6))
    ax = plt.axes(projection=ccrs.PlateCarree())

    ensemble_mean_r20.plot(
        ax=ax,
        transform=ccrs.PlateCarree(),
        cmap='PuBuGn',
        cbar_kwargs={'label': 'Days with PR ≥ 20 mm'},
    )

    ax.coastlines()
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.set_title("Ensemble Mean: Rain Days ≥ 20 mm (Historical)", fontsize=14)
    ax.set_extent([16, 33, -35, -22])
    plt.tight_layout()
    plt.show()