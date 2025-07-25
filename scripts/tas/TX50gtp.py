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
def compute_tx50p(tasmax, base_period):
    base = tasmax.sel(time=base_period)
    if base.time.size == 0:
        raise ValueError("⚠️ Base period is empty.")
    return base.groupby("time.dayofyear").reduce(np.percentile, q=50, dim="time")

def compute_above_median_percentage(tasmax, tx50p):
    doy = tasmax["time"].dt.dayofyear
    threshold = tx50p.sel(dayofyear=doy)
    hot_days = tasmax > threshold
    annual_pct = 100 * hot_days.groupby("time.year").mean(dim="time")
    return annual_pct

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
        ds_list = [xr.open_dataset(f) for f in nc_files]
        tasmax = xr.concat([ds["tasmax"] for ds in ds_list], dim="time") - 273.15
        tasmax = tasmax.sel(**sa_bounds)

        if tasmax.time.size < 365 * 30:
            print(f"⚠️ Skipping {model} due to insufficient time coverage.")
            continue

        tx50p = compute_tx50p(tasmax, base_period)
        annual_pct = compute_above_median_percentage(tasmax, tx50p)

        # National average
        model_pct_national.append(annual_pct.mean(dim=["lat", "lon"]))

        # Bioregional average
        region_ids = region_mask.mask(tasmax)
        bio_means = [annual_pct.where(region_ids == r).mean(dim=["lat", "lon"], skipna=True)
                     for r in range(len(region_mask))]
        model_pct_bioregion.append(xr.concat(bio_means, dim="region"))

        for ds in ds_list:
            ds.close()

    except Exception as e:
        print(f"❌ Failed to process {model}: {e}")

# ------------------ Save and Visualize ------------------ #
if model_pct_national:
    ensemble_pct = xr.concat(model_pct_national, dim="model").mean(dim="model")
    ensemble_pct.name = "Above50Pct"
    ensemble_pct.attrs.update({
        "description": "Percentage of Days with TX > 50th percentile",
        "units": "%"
    })
    ensemble_pct.to_netcdf(os.path.join(input_dir, "ensemble_above50pct_tx_SA.nc"))
    print("✅ Saved national ensemble percentage for TX > 50th percentile.")
else:
    print("⚠️ No valid national data.")

if model_pct_bioregion:
    ensemble_bio_pct = xr.concat(model_pct_bioregion, dim="model").mean(dim=["model", "year"])
    bioregions["Above50Pct"] = ensemble_bio_pct.values

    fig, ax = plt.subplots(figsize=(10, 8))
    vmin = int(bioregions["Above50Pct"].min())
    vmax = int(bioregions["Above50Pct"].max())
    ticks = range(vmin, vmax + 5, 5)

    bioregions.plot(
        column="Above50Pct",
        cmap="OrRd",
        linewidth=0.8,
        edgecolor="black",
        legend=True,
        legend_kwds={"label": "% Days TX > 50th percentile", "orientation": "vertical", "ticks": ticks},
        ax=ax
    )

    ax.set_title("CMIP6: % of Days with TX > 50th Percentile by Bioregion (1950–2014)", fontsize=14)
    ax.set_axis_off()
    plt.tight_layout()
    plt.show()

