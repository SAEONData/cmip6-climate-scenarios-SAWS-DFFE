import os
import glob
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# ------------------ Configuration ------------------ #
data_path = "/home//caroline/nex-gddp-sa-tools/data/pr"  #  Replace with your actual path
lat_bounds = [-35, -22]
lon_bounds = [16, 33]
threshold = 10.0  # mm/day threshold

# ------------------ Find Only Historical Files ------------------ #
all_hist_files = sorted(glob.glob(os.path.join(data_path, "**", "historical", "*.nc"), recursive=True))
ensemble_r10 = []

# ------------------ Process Each File ------------------ #
for file in all_hist_files:
    try:
        ds = xr.open_dataset(file)

        # Convert PR to mm/day
        pr = ds['pr'] * 86400  # kg/m²/s → mm/day

        # Subset to region of interest
        pr = pr.sel(lat=slice(*lat_bounds), lon=slice(*lon_bounds))

        # (Optional) Time slice e.g., 1981–2010
        # pr = pr.sel(time=slice("1981-01-01", "2010-12-31"))

        # Count days where PR >= 10 mm
        r10 = (pr >= threshold).sum(dim="time")

        ensemble_r10.append(r10)
        print(f"✅ Processed: {os.path.basename(file)}")

    except Exception as e:
        print(f"❌ Error in {file}: {e}")

# ------------------ Ensemble Mean ------------------ #
if ensemble_r10:
    r10_stack = xr.concat(ensemble_r10, dim="model")
    ensemble_mean_r10 = r10_stack.mean(dim="model")
    ensemble_mean_r10.name = "R10mm"
    ensemble_mean_r10.to_netcdf("ensemble_r10mm_mean.nc")
    print("✅ Saved to 'ensemble_r10mm_mean.nc'")

# ------------------ Plot ------------------ #
    fig = plt.figure(figsize=(10, 6))
    ax = plt.axes(projection=ccrs.PlateCarree())

    ensemble_mean_r10.plot(
        ax=ax,
        transform=ccrs.PlateCarree(),
        cmap='Blues',
        cbar_kwargs={'label': 'Days with PR ≥ 10 mm'},
    )

    ax.coastlines()
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.set_title("CMIP6 Mean: Rain Days ≥ 10 mm (Historical)", fontsize=14)
    ax.set_extent([16, 33, -35, -22])
    plt.tight_layout()
    plt.show()
