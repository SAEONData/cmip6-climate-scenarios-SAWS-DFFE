import os
import glob
import pandas as pd
import numpy as np
import xarray as xr
import geopandas as gpd
import regionmask
import matplotlib.pyplot as plt

# ------------------ Config ------------------ #
data_path = "/content/drive/MyDrive/Climate Data TTT"
shapefile_path = "/content/drive/MyDrive/CLIMREG OVERLAPS REMOVED/cleaned_veg_biome_clim_reg.shp"
towns_csv_path = "/content/drive/MyDrive/CLIMREG OVERLAPS REMOVED/towns.csv"

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

# ------------------ Helper: Calculate R95pTOT (%) ------------------ #
def calculate_r95ptot(precip, wet_day_threshold=1.0):
    # Mask wet days (≥ 1.0 mm/day)
    wet_days = precip.where(precip >= wet_day_threshold, drop=True)

    # Compute 95th percentile for each grid point
    pr95 = wet_days.quantile(0.95, dim='time')

    # Sum of PR on days > 95th percentile (very wet days)
    pr_very_wet_sum = wet_days.where(wet_days > pr95).sum(dim='time')

    # Total wet-day rainfall sum
    pr_total_wet_sum = wet_days.sum(dim='time')

    # R95pTOT as %
    r95ptot_percent = (pr_very_wet_sum / pr_total_wet_sum) * 100.0

    return r95ptot_percent

# ------------------ Find files ------------------ #
all_hist_files = sorted(glob.glob(os.path.join(data_path, "**", "historical", "*.nc"), recursive=True))
print(f"📂 Found {len(all_hist_files)} historical NetCDF files.")

# ------------------ Process files ------------------ #
bioregion_r95ptot_data = []
model_names = []

for file in all_hist_files:
    try:
        ds = xr.open_dataset(file)

        # Subset to South Africa
        ds = ds.sel(lat=slice(*lat_bounds), lon=slice(*lon_bounds))

        # Convert precipitation to mm/day
        pr = ds['pr'] * 86400.0

        # Resample to annual
        pr_r95ptot_annual = pr.resample(time='YE').map(calculate_r95ptot)

        # Long-term mean across years
        pr_r95ptot_mean = pr_r95ptot_annual.mean(dim='time')

        # Apply region mask
        region_mask_da = region_mask.mask(pr_r95ptot_mean)

        # Aggregate by region
        regional_r95ptot = pr_r95ptot_mean.groupby(region_mask_da).mean()

        # Append to list
        bioregion_r95ptot_data.append(regional_r95ptot)
        model_name = file.split(os.sep)[-3]  # e.g. CESM2
        model_names.append(model_name)

        print(f"✅ Processed: {os.path.basename(file)}")

    except Exception as e:
        print(f"❌ Error processing {os.path.basename(file)}: {e}")
        continue

# ------------------------ Ensemble Mean and Output ------------------------ #
ensemble_stack = xr.concat(bioregion_r95ptot_data, dim='model')
ensemble_stack['model'] = model_names

# Compute ensemble mean per region
ensemble_mean_by_region = ensemble_stack.mean(dim='model')

# Convert to DataFrame
df = pd.DataFrame({
    "Bioregion": region_mask.names,
    "R95pTOT_percent": ensemble_mean_by_region.values
}).sort_values(by="R95pTOT_percent", ascending=False)

# Print results
print("\n📊 R95pTOT: % of wet-day rainfall from PR > 95th percentile (Ensemble Average):")
print(df)

# Merge mean values with bioregions GeoDataFrame
bioregion_mean_df = pd.DataFrame({
    "Veg_Biome": region_mask.names,
    "R95pTOT_percent": ensemble_mean_by_region.values
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
vmax = bioregions_merged["R95pTOT_percent"].max()
step = 5
ticks = np.arange(vmin, vmax + step, step)

bioregions_merged.plot(
    column="R95pTOT_percent",
    cmap="Blues",
    linewidth=0.8,
    edgecolor="black",
    legend=True,
    legend_kwds={
        "label": "R95pTOT (%) - Share of rainfall from very wet days",
        "orientation": "vertical",
        "ticks": ticks
    },
    ax=ax
)

towns_gdf.plot(ax=ax, color='red', markersize=40, zorder=5)

for x, y, label in zip(towns_gdf.geometry.x, towns_gdf.geometry.y, towns_gdf['city']):
    ax.text(x, y, label, fontsize=9, ha='left', va='bottom')

ax.set_title("CMIP6: % of wet-day P PR > 95th percentile (1950-2014)", fontsize=14)
ax.set_axis_off()
plt.tight_layout()
plt.show()
