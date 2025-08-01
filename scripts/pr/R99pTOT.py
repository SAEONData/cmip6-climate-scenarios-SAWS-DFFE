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
input_dir = "/content/drive/MyDrive/Climate Data TTT/pr"
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
def compute_pr99(pr, base_period):
    base = pr.sel(time=base_period)
    if base.time.size == 0:
        raise ValueError("⚠️ Base period is empty.")
    return base.groupby("time.dayofyear").reduce(np.percentile, q=99, dim="time")

def compute_fraction_exceeding_99th(pr, pr99):
    doy = pr["time"].dt.dayofyear
    threshold = pr99.sel(dayofyear=doy)
    wet_days = pr > 1.0  # wet-day threshold
    extreme_days = pr > threshold

    total_wet_day_rain = pr.where(wet_days).groupby("time.year").sum(dim="time")
    extreme_wet_day_rain = pr.where(extreme_days).groupby("time.year").sum(dim="time")
    
    fraction = 100 * extreme_wet_day_rain / total_wet_day_rain
    return fraction

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
        pr = xr.concat([ds["pr"] for ds in ds_list], dim="time") * 86400  # Convert from kg/m2/s to mm/day
        pr = pr.sel(**sa_bounds)

        if pr.time.size < 365 * 30:
            print(f"⚠️ Skipping {model} due to insufficient time coverage.")
            continue

        pr99 = compute_pr99(pr, base_period)
        frac_99 = compute_fraction_exceeding_99th(pr, pr99)

        # National average
        model_pct_national.append(frac_99.mean(dim=["lat", "lon"]))

        # Bioregional average
        region_ids = region_mask.mask(pr)
        bio_means = [frac_99.where(region_ids == r).mean(dim=["lat", "lon"], skipna=True)
                     for r in range(len(region_mask))]
        model_pct_bioregion.append(xr.concat(bio_means, dim="region"))

        for ds in ds_list:
            ds.close()

    except Exception as e:
        print(f"❌ Failed to process {model}: {e}")

# ------------------ Save and Visualize ------------------ #
if model_pct_national:
    ensemble_pct = xr.concat(model_pct_national, dim="model").mean(dim="model")
    ensemble_pct.name = "WetDayPR_Above99thPct"
    ensemble_pct.attrs.update({
        "description": "Fraction of Wet-Day Rainfall from PR > 99th percentile",
        "units": "%"
    })
    ensemble_pct.to_netcdf(os.path.join(input_dir, "ensemble_wetday_99th_pct_SA.nc"))
    print("✅ Saved national ensemble fraction.")

if model_pct_bioregion:
    ensemble_bio_pct = xr.concat(model_pct_bioregion, dim="model").mean(dim=["model", "year"])
    bioregions["Wet99Pct"] = ensemble_bio_pct.values

    fig, ax = plt.subplots(figsize=(10, 8))
    vmin = int(bioregions["Wet99Pct"].min())
    vmax = int(bioregions["Wet99Pct"].max())
    ticks = range(vmin, vmax + 5, 5)

    bioregions.plot(
        column="Wet99Pct",
        cmap="Blues",
        linewidth=0.8,
        edgecolor="black",
        legend=True,
        legend_kwds={"label": "% of Rainfall (PR > 99th percentile)", "orientation": "vertical", "ticks": ticks},
        ax=ax
    )

    ax.set_title("CMIP6: Fraction of Wet-Day Rainfall from PR > 99th Percentile by Bioregion (1950–2014)", fontsize=14)
    ax.set_axis_off()
    plt.tight_layout()
    plt.show()
