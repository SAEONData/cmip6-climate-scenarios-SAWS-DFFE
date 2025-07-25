import os
import glob
import numpy as np
import xarray as xr
import geopandas as gpd
import regionmask
from tqdm import tqdm
import matplotlib.pyplot as plt
import warnings
import cftime
warnings.filterwarnings("ignore")

# ------------------ Configuration ------------------ #
input_dir = "/content/drive/MyDrive/Climate Data TTT/models"
shapefile_path = "/content/drive/MyDrive/CLIMREG OVERLAPS REMOVED/cleaned_veg_biome_clim_reg.shp"
base_period = slice("1950-01-01", "2014-12-31")
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
def compute_tx90p(tasmax, base_period):
    base = tasmax.sel(time=base_period)
    if base.time.size == 0:
        raise ValueError("⚠️ Base period is empty.")

    # Robust dayofyear extraction for all calendars
    time_vals = base.indexes["time"].to_list()
    doy = xr.DataArray(
        [t.timetuple().tm_yday for t in time_vals],
        coords={"time": base.time},
        dims="time"
    )

    base.coords["doy"] = doy
    tx90p = base.groupby("doy").reduce(np.percentile, q=90, dim="time")
    return tx90p.rename({"doy": "dayofyear"})

def compute_hot_day_percentage(tasmax, tx90p):
    # Extract dayofyear for threshold lookup
    time_vals = tasmax.indexes["time"].to_list()
    doy = xr.DataArray(
        [t.timetuple().tm_yday for t in time_vals],
        coords={"time": tasmax.time},
        dims="time"
    )
    tasmax.coords["doy"] = doy

    threshold = tx90p.sel(dayofyear=doy)
    hot_days = tasmax > threshold

    # Extract year manually
    years = xr.DataArray(
        [t.year for t in time_vals],
        coords={"time": tasmax.time},
        dims="time"
    )
    hot_days.coords["year"] = years

    # Group by extracted year
    annual_hot_pct = 100 * hot_days.groupby("year").mean(dim="time")
    return annual_hot_pct

# ------------------ Main Processing ------------------ #
model_pct_bioregion = []
model_pct_national = []
model_dirs = sorted([d for d in os.listdir(input_dir) if os.path.isdir(os.path.join(input_dir, d))])
print(f"📂 Found {len(model_dirs)} model folders.")

for model in tqdm(model_dirs, desc="🔁 Processing models"):
    model_path = os.path.join(input_dir, model)
    nc_files = sorted(glob.glob(os.path.join(model_path, "*.nc")))
    if not nc_files:
        print(f"⚠️ No NetCDF files found in {model_path}")
        continue

    try:
        ds_list = [xr.open_dataset(f, decode_times=True, use_cftime=True) for f in nc_files]
        tasmax = xr.concat([ds["tasmax"] for ds in ds_list], dim="time") - 273.15
        tasmax = tasmax.sel(**sa_bounds)

        if tasmax.time.size < 365 * 30:
            print(f"⚠️ Skipping {model} due to insufficient time coverage.")
            continue

        tx90p = compute_tx90p(tasmax, base_period)
        annual_hot_pct = compute_hot_day_percentage(tasmax, tx90p)

        # National average
        model_pct_national.append(annual_hot_pct.mean(dim=["lat", "lon"]))

        # Bioregional average
        region_ids = region_mask.mask(tasmax)
        bio_means = [annual_hot_pct.where(region_ids == r).mean(dim=["lat", "lon"], skipna=True)
                     for r in range(len(region_mask))]
        model_pct_bioregion.append(xr.concat(bio_means, dim="region"))

        for ds in ds_list:
            ds.close()

    except Exception as e:
        print(f"❌ Failed to process {model}: {e}")

# ------------------ Save and Visualize ------------------ #
if model_pct_national:
    ensemble_pct = xr.concat(model_pct_national, dim="model").mean(dim="model")
    ensemble_pct.name = "HotDaysPct"
    ensemble_pct.attrs.update({
        "description": "Percentage of Hot Days (TX > 90th percentile)",
        "units": "%"
    })
    ensemble_pct.to_netcdf(os.path.join(input_dir, "ensemble_hot_days_pct_SA.nc"))
    print("✅ Saved national ensemble hot day percentage.")
else:
    print("⚠️ No valid national data.")

if model_pct_bioregion:
    ensemble_bio_pct = xr.concat(model_pct_bioregion, dim="model").mean(dim=["model", "year"])
    bioregions["HotDaysPct"] = ensemble_bio_pct.values

    fig, ax = plt.subplots(figsize=(10, 8))
    vmin = int(bioregions["HotDaysPct"].min())
    vmax = int(bioregions["HotDaysPct"].max())
    ticks = range(vmin, vmax + 5, 5)

    bioregions.plot(
        column="HotDaysPct",
        cmap="Reds",
        linewidth=0.8,
        edgecolor="black",
        legend=True,
        legend_kwds={"label": "% of Hot Days", "orientation": "vertical", "ticks": ticks},
        ax=ax
    )

    ax.set_title("CMIP6: % of Days with TX > 90th Percentile by Bioregion (1950–2014)", fontsize=14)
    ax.set_axis_off()
    plt.tight_layout()
    plt.show()