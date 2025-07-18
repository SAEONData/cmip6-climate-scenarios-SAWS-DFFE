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
input_dir = "/content/drive/MyDrive/Climate Data TTT/tasmin"
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
def compute_tn90p(tasmin, base_period):
    base = tasmin.sel(time=base_period)
    if base.time.size == 0:
        raise ValueError("⚠️ Base period is empty.")

    time_vals = base.indexes["time"].to_list()
    doy = xr.DataArray(
        [t.timetuple().tm_yday for t in time_vals],
        coords={"time": base.time},
        dims="time"
    )

    base.coords["doy"] = doy
    tn90p = base.groupby("doy").reduce(np.percentile, q=90, dim="time")
    return tn90p.rename({"doy": "dayofyear"})

def compute_warm_night_percentage(tasmin, tn90p):
    time_vals = tasmin.indexes["time"].to_list()
    doy = xr.DataArray(
        [t.timetuple().tm_yday for t in time_vals],
        coords={"time": tasmin.time},
        dims="time"
    )
    tasmin.coords["doy"] = doy

    threshold = tn90p.sel(dayofyear=doy)
    warm_nights = tasmin > threshold

    years = xr.DataArray(
        [t.year for t in time_vals],
        coords={"time": tasmin.time},
        dims="time"
    )
    warm_nights.coords["year"] = years

    annual_pct = 100 * warm_nights.groupby("year").mean(dim="time")
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
        ds_list = [xr.open_dataset(f, decode_times=True, use_cftime=True) for f in nc_files]
        tasmin = xr.concat([ds["tasmin"] for ds in ds_list], dim="time") - 273.15
        tasmin = tasmin.sel(**sa_bounds)

        if tasmin.time.size < 365 * 30:
            print(f"⚠️ Skipping {model} due to insufficient time coverage.")
            continue

        tn90p = compute_tn90p(tasmin, base_period)
        annual_pct = compute_warm_night_percentage(tasmin, tn90p)

        # National average
        model_pct_national.append(annual_pct.mean(dim=["lat", "lon"]))

        # Bioregional average
        region_ids = region_mask.mask(tasmin)
        bio_means = [annual_pct.where(region_ids == r).mean(dim=["lat", "lon"], skipna=True)
                     for r in range(len(region_mask))]
        model_pct_bioregion.append(xr.concat(bio_means, dim="region"))

        for ds in ds_list:
            ds.close()

    except Exception as e:
        print(f" Failed to process {model}: {e}")

# ------------------ Save and Visualize ------------------ #
if model_pct_national:
    ensemble_pct = xr.concat(model_pct_national, dim="model").mean(dim="model")
    ensemble_pct.name = "WarmNightsPct"
    ensemble_pct.attrs.update({
        "description": "Percentage of Warm Nights (TN > 90th percentile)",
        "units": "%"
    })
    ensemble_pct.to_netcdf(os.path.join(input_dir, "ensemble_warm_nights_pct_SA.nc"))
    print("✅ Saved national ensemble warm night percentage.")
else:
    print("⚠️ No valid national data.")

if model_pct_bioregion:
    ensemble_bio_pct = xr.concat(model_pct_bioregion, dim="model").mean(dim=["model", "year"])
    bioregions["WarmNightsPct"] = ensemble_bio_pct.values

    fig, ax = plt.subplots(figsize=(10, 8))
    vmin = int(bioregions["WarmNightsPct"].min())
    vmax = int(bioregions["WarmNightsPct"].max())
    ticks = range(vmin, vmax + 5, 5)

    bioregions.plot(
        column="WarmNightsPct",
        cmap="Oranges",
        linewidth=0.8,
        edgecolor="black",
        legend=True,
        legend_kwds={"label": "% of Warm Nights", "orientation": "vertical", "ticks": ticks},
        ax=ax
    )

    ax.set_title("CMIP6: % of Days with TN > 90th Percentile by Bioregion (1950–2014)", fontsize=14)
    ax.set_axis_off()
    plt.tight_layout()
    plt.show()
