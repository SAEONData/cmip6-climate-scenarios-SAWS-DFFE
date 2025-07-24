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
input_dir = "/home/caroline/nex-gddp-sa-tools/data/pr"
shapefile_path = "/home/caroline/nex-gddp-sa-tools/climate_regions/cleaned_veg_biome_clim_reg.shp"
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
def compute_pr99p(pr, base_period):
    try:
        base = pr.sortby("time").sel(time=base_period)
        wet_base = base.where(base >= 1.0)
        return wet_base.reduce(np.nanpercentile, q=99, dim="time")
    except Exception as e:
        raise ValueError(f"Error computing 99th percentile during base period: {e}")

def compute_sum_above_99p(pr, pr99p):
    extreme = pr.where(pr > pr99p)
    return extreme.groupby("time.year").sum(dim="time", skipna=True)

# ------------------ Main Processing ------------------ #
model_sum_bioregion = []
model_sum_national = []
model_dirs = sorted([d for d in os.listdir(input_dir) if os.path.isdir(os.path.join(input_dir, d))])
print(f"📂 Found {len(model_dirs)} model folders.")

for model in tqdm(model_dirs, desc="🔁 Processing models"):
    model_path = os.path.join(input_dir, model)
    nc_files = sorted(glob.glob(os.path.join(model_path, "**", "*.nc"), recursive=True))

    if not nc_files:
        print(f"⚠️ No NetCDF files found in {model_path}")
        continue

    try:
        ds_list = [xr.open_dataset(f) for f in nc_files]
        pr_list = []
        for ds in ds_list:
            if "pr" not in ds:
                raise ValueError("Missing 'pr' variable in dataset")
            pr_list.append(ds["pr"])

        pr = xr.concat(pr_list, dim="time") * 86400  # Convert from kg/m²/s to mm/day
        pr = pr.sortby("time").sel(**sa_bounds)

        if pr.time.size < 365 * 30:
            print(f"⚠️ Skipping {model} due to insufficient time coverage.")
            continue

        pr99p = compute_pr99p(pr, base_period)
        annual_sum = compute_sum_above_99p(pr, pr99p)

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
    ensemble_sum.name = "SumPR99p"
    ensemble_sum.attrs.update({
        "description": "Sum of Daily Precipitation > 99th Percentile",
        "units": "mm/year"
    })
    ensemble_sum.to_netcdf(os.path.join(input_dir, "ensemble_sum_pr99p_SA.nc"))
    print("✅ Saved national ensemble sum PR > 99th percentile.")
else:
    print("⚠️ No valid national data.")

if model_sum_bioregion:
    ensemble_bio_sum = xr.concat(model_sum_bioregion, dim="model").mean(dim=["model", "year"])
    bioregions["SumPR99p"] = ensemble_bio_sum.values

    fig, ax = plt.subplots(figsize=(10, 8))
    vmin = int(np.floor(bioregions["SumPR99p"].min()))
    vmax = int(np.ceil(bioregions["SumPR99p"].max()))
    step = max((vmax - vmin) // 5, 10)
    ticks = range(vmin, vmax + step, step)

    bioregions.plot(
        column="SumPR99p",
        cmap="YlGnBu",
        linewidth=0.8,
        edgecolor="black",
        legend=True,
        legend_kwds={"label": "Annual PR Sum > 99th (mm)", "orientation": "vertical", "ticks": ticks},
        ax=ax
    )

    ax.set_title("CMIP6: Annual Sum of PR > 99th Percentile by Bioregion (1950–2014)", fontsize=14)
    ax.set_axis_off()
    plt.tight_layout()
    plt.show()
