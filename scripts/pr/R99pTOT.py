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
def compute_pr99p(pr, base_period):
    base = pr.sel(time=base_period)
    wet_base = base.where(base >= 1.0)
    return wet_base.reduce(np.nanpercentile, q=99, dim="time")

def compute_fraction_above_99p(pr, pr99p):
    wet = pr.where(pr >= 1.0)
    heavy = pr.where(pr > pr99p)
    annual_wet_total = wet.groupby("time.year").sum(dim="time", skipna=True)
    annual_heavy_total = heavy.groupby("time.year").sum(dim="time", skipna=True)
    return annual_heavy_total / annual_wet_total

# ------------------ Main Processing ------------------ #
model_frac_bioregion = []
model_frac_national = []
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

        pr99p = compute_pr99p(pr, base_period)
        annual_frac = compute_fraction_above_99p(pr, pr99p)

        model_frac_national.append(annual_frac.mean(dim=["lat", "lon"]))

        region_ids = region_mask.mask(pr)
        bio_means = [annual_frac.where(region_ids == r).mean(dim=["lat", "lon"], skipna=True)
                     for r in range(len(region_mask))]
        model_frac_bioregion.append(xr.concat(bio_means, dim="region"))

        for ds in ds_list:
            ds.close()

    except Exception as e:
        print(f"❌ Failed to process {model}: {e}")

# ------------------ Save and Visualize ------------------ #
if model_frac_national:
    ensemble_frac = xr.concat(model_frac_national, dim="model").mean(dim="model")
    ensemble_frac.name = "FracWetPR99p"
    ensemble_frac.attrs.update({
        "description": "Fraction of Wet-Day Rainfall from Days with PR > 99th percentile",
        "units": "fraction"
    })
    ensemble_frac.to_netcdf(os.path.join(input_dir, "ensemble_frac_wet_pr99p_SA.nc"))
    print("✅ Saved national ensemble PR fraction > 99th percentile.")
else:
    print("⚠️ No valid national data.")

if model_frac_bioregion:
    ensemble_bio_frac = xr.concat(model_frac_bioregion, dim="model").mean(dim=["model", "year"])
    bioregions["FracWetPR99p"] = ensemble_bio_frac.values

    fig, ax = plt.subplots(figsize=(10, 8))
    vmin = 0
    vmax = round(bioregions["FracWetPR99p"].max(), 2)
    ticks = np.linspace(vmin, vmax, 6)

    bioregions.plot(
        column="FracWetPR99p",
        cmap="Purples",
        linewidth=0.8,
        edgecolor="black",
        legend=True,
        legend_kwds={"label": "Fraction PR > 99th on Wet Days", "orientation": "vertical", "ticks": ticks},
        ax=ax
    )

    ax.set_title("CMIP6: Fraction of Wet-Day PR from >99th Percentile by Bioregion (1950–2014)", fontsize=14)
    ax.set_axis_off()
    plt.tight_layout()
    plt.show()
