import os
import glob
import numpy as np
import xarray as xr
import geopandas as gpd
import regionmask
from tqdm import tqdm
import matplotlib.pyplot as plt
import warnings

warnings.filterwarnings("ignore")

# ------------------ Configuration ------------------ #
input_dir = "/content/drive/MyDrive/Climate Data TTT/tasmin" # Folder with tasmax & tasmin NetCDFs
shapefile_path = "/content/drive/MyDrive/CLIMREG OVERLAPS REMOVED/cleaned_veg_biome_clim_reg.shp"
base_period = slice("1950-01-01", "2014-12-31")
min_duration = 3  # User-defined consecutive day threshold
sa_bounds = dict(lat=slice(-35, -22), lon=slice(16, 33))

# ------------------ Load Bioregions and Region Mask ------------------ #
bioregions = gpd.read_file(shapefile_path).to_crs("EPSG:4326")
region_mask = regionmask.Regions(
    outlines=bioregions.geometry,
    names=bioregions['Veg_Biome'],
    abbrevs=bioregions['Veg_Biome'],
    name="Bioregions"
)

# ------------------ Helper Functions ------------------ #
def compute_daily_percentile_thresholds(data, base_period, q):
    base = data.sel(time=base_period)
    if base.time.size == 0:
        raise ValueError("⚠️ Base period is empty.")

    time_vals = base.indexes["time"].to_list()
    doy = xr.DataArray(
        [t.timetuple().tm_yday for t in time_vals],
        coords={"time": base.time},
        dims="time"
    )
    base.coords["doy"] = doy

    def percentile_func(x):
        result = np.percentile(x, q, axis=0)
        coords = {k: v for k, v in x.coords.items() if k not in ['time', 'doy']}
        return xr.DataArray(result, dims=('lat', 'lon'), coords=coords)

    return base.groupby("doy").apply(percentile_func).rename({"doy": "dayofyear"})

def identify_warm_spells(tx, tn, tx95p, tn95p, min_duration):
    time_vals = tx.indexes["time"].to_list()
    doy = xr.DataArray(
        [t.timetuple().tm_yday for t in time_vals],
        coords={"time": tx.time},
        dims="time"
    )
    tx.coords["doy"] = doy
    tn.coords["doy"] = doy

    tx_thresh = tx95p.sel(dayofyear=doy).drop_vars("dayofyear")
    tn_thresh = tn95p.sel(dayofyear=doy).drop_vars("dayofyear")
    warm_condition = (tx > tx_thresh) & (tn > tn_thresh)

    # ✅ Fixed: no duplicate 'dim' argument
    spells = warm_condition.rolling(time=min_duration).reduce(np.all)
    return spells.fillna(False)

# ------------------ Main Processing ------------------ #
model_dirs = sorted([d for d in os.listdir(input_dir) if os.path.isdir(os.path.join(input_dir, d))])
print(f"📂 Found {len(model_dirs)} model folders.")
model_bioregion_results = []

for model in tqdm(model_dirs, desc="🔁 Processing models"):
    model_path = os.path.join(input_dir, model)
    nc_files = sorted(glob.glob(os.path.join(model_path, "*.nc")))
    if not nc_files:
        print(f"⚠️ No NetCDF files found in {model_path}")
        continue

    try:
        tasmax_files = [f for f in nc_files if "tasmax" in f]
        tasmin_files = [f for f in nc_files if "tasmin" in f]
        if not tasmax_files or not tasmin_files:
            print(f"⚠️ Missing tasmax or tasmin in {model_path}")
            continue

        ds_tasmax = xr.open_mfdataset(tasmax_files, combine="by_coords", decode_times=True, use_cftime=True)
        ds_tasmin = xr.open_mfdataset(tasmin_files, combine="by_coords", decode_times=True, use_cftime=True)

        tasmax = ds_tasmax["tasmax"].sel(**sa_bounds) - 273.15
        tasmin = ds_tasmin["tasmin"].sel(**sa_bounds) - 273.15

        if tasmax.time.size < 365 * 30 or tasmin.time.size < 365 * 30:
            print(f"⚠️ Skipping {model} due to insufficient time coverage.")
            continue

        tx95p = compute_daily_percentile_thresholds(tasmax, base_period, 95)
        tn95p = compute_daily_percentile_thresholds(tasmin, base_period, 95)

        warm_spells = identify_warm_spells(tasmax, tasmin, tx95p, tn95p, min_duration)

        # Assign year per timestep (for cftime)
        time_vals = tasmax.indexes["time"].to_list()
        years = xr.DataArray(
            [t.year for t in time_vals],
            coords={"time": tasmax.time},
            dims="time"
        )
        warm_spells.coords["year"] = years
        spells_per_year = warm_spells.groupby("year").sum(dim="time")

        # Compute bioregion means
        region_ids = region_mask.mask(tasmax)
        bio_means = [
            spells_per_year.where(region_ids == r).mean(dim=["lat", "lon"], skipna=True)
            for r in range(len(region_mask))
        ]
        model_bioregion_results.append(xr.concat(bio_means, dim="region"))

        ds_tasmax.close()
        ds_tasmin.close()

    except Exception as e:
        print(f"❌ Failed to process {model}: {e}")

# ------------------ Aggregate & Plot ------------------ #
if model_bioregion_results:
    ensemble_bio_spells = xr.concat(model_bioregion_results, dim="model").mean(dim=["model", "year"])
    bioregions["HotSpells"] = ensemble_bio_spells.values

    fig, ax = plt.subplots(figsize=(10, 8))
    vmin = int(bioregions["HotSpells"].min())
    vmax = int(bioregions["HotSpells"].max())
    ticks = range(vmin, vmax + 1, 1)

    bioregions.plot(
        column="HotSpells",
        cmap="OrRd",
        linewidth=0.8,
        edgecolor="black",
        legend=True,
        legend_kwds={
            "label": f"Mean Hot Spells (TX & TN > 95th, ≥{min_duration} days)",
            "orientation": "vertical",
            "ticks": ticks
        },
        ax=ax
    )

    ax.set_title(f"CMIP6: Mean Hot Spells per Year by Bioregion\n(TX & TN > 95th Percentile, ≥{min_duration} Consecutive Days)", fontsize=14)
    ax.set_axis_off()
    plt.tight_layout()
    plt.show()