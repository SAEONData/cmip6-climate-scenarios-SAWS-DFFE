import os
import glob
import pandas as pd
import numpy as np
import xarray as xr
import geopandas as gpd
import regionmask
import matplotlib.pyplot as plt

# ------------------ Config ------------------ #
data_path = "/home/caroline/nex-gddp-sa-tools/data/pr"
shapefile_path = "/home/caroline/nex-gddp-sa-tools/climate_regions/cleaned_veg_biome_clim_reg.shp"
towns_csv_path = "/home/caroline/nex-gddp-sa-tools/cities/cities.csv"

lat_bounds = [-35, -22]
lon_bounds = [16, 33]

# ------------------ Load Bioregions ------------------ #
bioregions = gpd.read_file(shapefile_path).to_crs("EPSG:4326")

region_mask = regionmask.Regions(
    outlines=bioregions.geometry,
    names=bioregions['Veg_Biome'],
    abbrevs=bioregions['Veg_Biome'],
    name="Bioregions"
)
# ------------------ Helper: Daily Precipitation Intensity ------------------ #
def daily_precip_intensity(precip, threshold=1.0):
    wet_days = precip.where(precip >= threshold)
    wet_day_sum = wet_days.sum(dim='time')
    wet_day_count = wet_days.count(dim='time')
    sdii = wet_day_sum / wet_day_count
    return sdii
# ------------------ Find files ------------------ #
all_hist_files = sorted(glob.glob(os.path.join(data_path, "**", "historical", "*.nc"), recursive=True))
print(f"📂 Found {len(all_hist_files)} historical NetCDF files.")

# ------------------ Process files ------------------ #
bioregion_sdii_data = []
model_names = []

for file in all_hist_files:
    try:
        ds = xr.open_dataset(file)

        # Subset to South Africa
        ds = ds.sel(lat=slice(*lat_bounds), lon=slice(*lon_bounds))

        # Convert precipitation to mm/day
        pr = ds['pr'] * 86400.0

        # Resample to annual
        sdii_annual = pr.resample(time='YE').map(daily_precip_intensity)

        # Long-term mean across years
        sdii_mean = sdii_annual.mean(dim='time')

        # Apply region mask
        region_mask_da = region_mask.mask(sdii_mean)

        # Aggregate by region
        regional_sdii = sdii_mean.groupby(region_mask_da).mean()

        # Append to list
        bioregion_sdii_data.append(regional_sdii)
        model_name = file.split(os.sep)[-3]  # e.g. CESM2
        model_names.append(model_name)

        print(f"✅ Processed: {os.path.basename(file)}")

    except Exception as e:
        print(f"❌ Error processing {os.path.basename(file)}: {e}")
        continue

# ------------------------ Ensemble Mean and Output ------------------------ #
ensemble_stack = xr.concat(bioregion_sdii_data, dim='model')
ensemble_stack['model'] = model_names

# Compute ensemble mean per region
ensemble_mean_by_region = ensemble_stack.mean(dim='model')

# Convert to DataFrame
df = pd.DataFrame({
    "Bioregion": region_mask.names,
    "SDII_mm_per_day": ensemble_mean_by_region.values
}).sort_values(by="SDII_mm_per_day", ascending=False)

# Print results
print("\n📊 Daily Precipitation Intensity (SDII, mm/day) (Ensemble Average):")
print(df)

# Merge with shapefile
bioregion_mean_df = pd.DataFrame({
    "Veg_Biome": region_mask.names,
    "SDII_mm_per_day": ensemble_mean_by_region.values
})

bioregions_merged = bioregions.merge(bioregion_mean_df, on="Veg_Biome")

# ------------------- Plot ------------------- #
towns_df = pd.read_csv(towns_csv_path, sep=';')
towns_df.columns = towns_df.columns.str.strip()

towns_gdf = gpd.GeoDataFrame(
    towns_df,
    geometry=gpd.points_from_xy(towns_df['lng'], towns_df['lat']),
    crs="EPSG:4326"
)

fig, ax = plt.subplots(figsize=(10, 8))

vmin = 0
vmax = bioregions_merged["SDII_mm_per_day"].max()
step = 2
ticks = np.arange(vmin, vmax + step, step)

bioregions_merged.plot(
    column="SDII_mm_per_day",
    cmap="Blues",
    linewidth=0.8,
    edgecolor="black",
    legend=True,
    legend_kwds={
        "label": "SDII (mm/day) – Mean Daily Precip on Wet Days",
        "orientation": "vertical",
        "ticks": ticks
    },
    ax=ax
)

towns_gdf.plot(ax=ax, color='red', markersize=40, zorder=5)

for x, y, label in zip(towns_gdf.geometry.x, towns_gdf.geometry.y, towns_gdf['city']):
    ax.text(x, y, label, fontsize=9, ha='left', va='bottom')

ax.set_title("CMIP6: Daily Precipitation Intensity (SDII, 1950–2014)", fontsize=14)
ax.set_axis_off()
plt.tight_layout()
plt.show()