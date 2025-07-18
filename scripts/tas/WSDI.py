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
data_path = "/content/drive/MyDrive/Climate Data TTT/models"
shapefile_path = "/content/drive/MyDrive/CLIMREG OVERLAPS REMOVED/cleaned_veg_biome_clim_reg.shp"
base_period = slice("1950-01-01", "2014-12-31")
min_spell_length = 6
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
def compute_daily_tx90p(tasmax, base_period):
    base = tasmax.sel(time=base_period)
    if base.time.size == 0:
        raise ValueError("⚠️ Base period is empty for this model.")
    return base.groupby("time.dayofyear").reduce(np.percentile, q=90, dim="time")

def get_warm_spell_mask(tasmax, tx90p, min_length):
    doy = tasmax['time'].dt.dayofyear
    threshold = tx90p.sel(dayofyear=doy)
    warm = tasmax > threshold
    warm_np = warm.values
    spell_mask = np.zeros_like(warm_np)

    for lat in range(warm.shape[1]):
        for lon in range(warm.shape[2]):
            series = warm_np[:, lat, lon]
            run = 0
            for t in range(len(series)):
                if series[t]:
                    run += 1
                else:
                    if run >= min_length:
                        spell_mask[t - run:t, lat, lon] = 1
                    run = 0
            if run >= min_length:
                spell_mask[len(series) - run:, lat, lon] = 1

    return xr.DataArray(spell_mask, coords=tasmax.coords, dims=tasmax.dims)

def compute_annual_wsdi(mask):
    return mask.groupby("time.year").sum(dim="time")

# ------------------ Main Processing ------------------ #
model_wsdi_bioregion = []
model_annual_wsdi = []
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

        tx90p = compute_daily_tx90p(tasmax, base_period)
        spell_mask = get_warm_spell_mask(tasmax, tx90p, min_length=min_spell_length)
        annual_wsdi = compute_annual_wsdi(spell_mask)

        # National WSDI
        model_annual_wsdi.append(annual_wsdi.mean(dim=["lat", "lon"]))

        # Bioregional WSDI
        region_ids = region_mask.mask(tasmax)
        bio_means = [annual_wsdi.where(region_ids == r).mean(dim=["lat", "lon"], skipna=True)
                     for r in range(len(region_mask))]
        model_wsdi_bioregion.append(xr.concat(bio_means, dim="region"))

        for ds in ds_list:
            ds.close()

    except Exception as e:
        print(f"❌ Failed to process {model}: {e}")

# ------------------ Save and Visualize ------------------ #
if model_annual_wsdi:
    ensemble_wsdi = xr.concat(model_annual_wsdi, dim="model").mean(dim="model")
    ensemble_wsdi.name = "WSDI"
    ensemble_wsdi.attrs.update({
        "description": "Annual Warm Spell Duration Indicator (TX > 90th percentile for ≥6 days)",
        "units": "days"
    })
    ensemble_wsdi.to_netcdf(os.path.join(input_dir, "ensemble_WSDI_annual_avg_SA.nc"))
    print("✅ Saved national ensemble WSDI.")
else:
    print("⚠️ No valid national WSDI data.")

if model_wsdi_bioregion:
    ensemble_wsdi_bioregion = xr.concat(model_wsdi_bioregion, dim="model").mean(dim=["model", "year"])
    bioregions["WSDI"] = ensemble_wsdi_bioregion.values

    fig, ax = plt.subplots(figsize=(10, 8))
    vmin = int(bioregions["WSDI"].min())
    vmax = int(bioregions["WSDI"].max())
    ticks = range(vmin, vmax + 5, 5)

    bioregions.plot(
        column="WSDI",
        cmap="YlOrRd",
        linewidth=0.8,
        edgecolor="black",
        legend=True,
        legend_kwds={"label": "Mean WSDI (Days)", "orientation": "vertical", "ticks": ticks},
        ax=ax
    )

    ax.set_title("CMIP6: WSDI, TX > 90th percentile for ≥ 6 days by Bioregion (1950–2014)", fontsize=14)
    ax.set_axis_off()
    plt.tight_layout()
    plt.show()

