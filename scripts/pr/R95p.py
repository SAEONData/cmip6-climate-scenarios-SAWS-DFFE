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
def compute_pr95p(pr, base_period):
    base = pr.sel(time=base_period)
    wet_base = base.where(base >= 1.0)
    return wet_base.reduce(np.nanpercentile, q=95, dim="time")

def compute_sum_above_95p(pr, pr95p):
    heavy = pr.where(pr > pr95p)
    return heavy.groupby("time.year").sum(dim="time", skipna=True)

# ------------------ Main Processing ------------------ #
model_sum_bioregion = []
model_sum_national = []
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
        pr = xr.concat([ds["pr"] for ds in ds_list], dim="time") * 86400  # kg/m²/s to mm/day
        pr = pr.sel(**sa_bounds)

        if pr.time.size < 365 * 30:
            print(f"⚠️ Skipping {model} due to insufficient time coverage.")
            continue

        pr95p = compute_pr95p(pr, base_period)
        annual_sum = compute_sum_above_95p(pr, pr95p)

        model_sum_national.append(annual_sum.mean(dim=["lat", "lon"]))

        region_ids = region_mask.mask(pr)
        bio_means = [annual_sum.where(region_ids == r).mean(dim=["lat", "lon"], skipna=True)
                     for r in range(len(region_mask))]
        model_sum_bioregion.append(xr.concat(bio_means, dim="region"))

        for ds in ds_list:
            ds.close()

    except Exception as e:
        print(f"❌ Failed to process {model}: {e}")

# ------------------ Save and Visualize ------------------ #
if model_sum_national:
    ensemble_sum = xr.concat(model_sum_national, dim="model").mean(dim="model")
    ensemble_sum.name = "SumPR95p"
    ensemble_sum.attrs.update({
        "description": "Sum of Daily Precipitation > 95th Percentile",
        "units": "mm/year"
    })
    ensemble_sum.to_netcdf(os.path.join(input_dir, "ensemble_sum_pr95p_SA.nc"))
    print("✅ Saved national ensemble sum PR > 95th percentile.")
else:
    print("⚠️ No valid national data.")

if model_sum_bioregion:
    ensemble_bio_sum = xr.concat(model_sum_bioregion, dim="model").mean(dim=["model", "year"])
    bioregions["SumPR95p"] = ensemble_bio_sum.values

    fig, ax = plt.subplots(figsize=(10, 8))
    vmin = int(bioregions["SumPR95p"].min())
    vmax = int(bioregions["SumPR95p"].max())
    step = max((vmax - vmin) // 5, 10)
    ticks = range(vmin, vmax + step, step)

    bioregions.plot(
        column="SumPR95p",
        cmap="YlGnBu",
        linewidth=0.8,
        edgecolor="black",
        legend=True,
        legend_kwds={"label": "Annual PR Sum > 95th (mm)", "orientation": "vertical", "ticks": ticks},
        ax=ax
    )

    ax.set_title("CMIP6: Annual Sum of PR > 95th Percentile by Bioregion (1950–2014)", fontsize=14)
    ax.set_axis_off()
    plt.tight_layout()
    plt.show()
