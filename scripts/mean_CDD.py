import os
import glob
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# ------------------ Configuration ------------------ #
data_path = "/home/caroline/nex-gddp-sa-tools/data/pr"  # Replace with your path
lat_bounds = [-35, -22]
lon_bounds = [16, 33]
threshold = 1.0  # mm/day

# ------------------ Helper Function ------------------ #
def max_consecutive_dry_days(precip, threshold=1.0):
    dry = precip < threshold
    dry = dry.astype(int)

    def _count(arr):
        # Count max consecutive ones
        padded = np.concatenate(([0], arr, [0]))
        diff = np.diff(padded)
        run_starts = np.where(diff == 1)[0]
        run_ends = np.where(diff == -1)[0]
        return (run_ends - run_starts).max() if run_starts.size else 0

    return xr.apply_ufunc(
        np.vectorize(_count, signature="(t)->()"),
        dry,
        input_core_dims=[["time"]],
        vectorize=True,
        dask="parallelized",
        output_dtypes=[int]
    )
# ------------------ Load and Process Files ------------------ #
all_hist_files = sorted(glob.glob(os.path.join(data_path, "**", "historical", "*.nc"), recursive=True))
ensemble_cdd = []

for file in all_hist_files:
    try:
        ds = xr.open_dataset(file)
        pr = ds['pr'] * 86400  # kg/m2/s to mm/day

        # Subset region
        pr = pr.sel(lat=slice(*lat_bounds), lon=slice(*lon_bounds))

        # Compute CDD per grid cell
        cdd = max_consecutive_dry_days(pr, threshold=threshold)

        ensemble_cdd.append(cdd)
        print(f"✅ Processed: {os.path.basename(file)}")

    except Exception as e:
        print(f"❌ Error processing {file}: {e}")

# ------------------ Compute Ensemble Mean ------------------ #
if ensemble_cdd:
    cdd_stack = xr.concat(ensemble_cdd, dim="model")
    ensemble_mean_cdd = cdd_stack.mean(dim="model")
    ensemble_mean_cdd.name = "CDD"
    ensemble_mean_cdd.to_netcdf("ensemble_cdd_mean.nc")
    print("✅ Saved to ensemble_cdd_mean.nc")

 # ------------------ Plot Map ------------------ #
    fig = plt.figure(figsize=(10, 6))
    ax = plt.axes(projection=ccrs.PlateCarree())

    ensemble_mean_cdd.plot(
        ax=ax,
        transform=ccrs.PlateCarree(),
        cmap='YlOrBr',
        cbar_kwargs={'label': 'Max Consecutive Dry Days'},
    )

    ax.coastlines()
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.set_title("CMIP6 Mean Annual Consecutive Dry Days CDD (Historical)", fontsize=14)
    ax.set_extent([16, 33, -35, -22])
    plt.tight_layout()
    plt.show()
else:
    print("⚠️ No historical files processed. Check your file paths.")
