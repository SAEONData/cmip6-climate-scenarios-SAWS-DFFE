import os
import glob
import pandas as pd
import xarray as xr
import geopandas as gpd
import regionmask
import matplotlib.pyplot as plt

# Define paths
data_path = "/home/caroline/nex-gddp-sa-tools/data/pr"
shapefile_path = "/home/caroline/nex-gddp-sa-tools/climate_regions/cleaned_veg_biome_clim_reg.shp"
towns_csv_path = "/home/caroline/nex-gddp-sa-tools/cities/cities.csv"

# South Africa bounding box (lat, lon)
lat_bounds = [-35, -22]
lon_bounds = [16, 33]

# Load shapefile and convert to WGS84
bioregions = gpd.read_file(shapefile_path).to_crs("EPSG:4326")

# Create region mask
region_mask = regionmask.Regions(
    outlines=bioregions.geometry,
    names=bioregions['Veg_Biome'],
    abbrevs=bioregions['Veg_Biome'],
    name="Bioregions"
)

# Find all historical NetCDF files recursively
all_hist_files = sorted(glob.glob(os.path.join(data_path, "**", "historical", "*.nc"), recursive=True))
print(f"📂 Found {len(all_hist_files)} historical NetCDF files.")

# Prepare containers for results
bioregion_annual_data = []
model_names = []

# Loop through each file
for file in all_hist_files:
    try:
        ds = xr.open_dataset(file)

        # Subset to South Africa
        ds = ds.sel(lat=slice(*lat_bounds), lon=slice(*lon_bounds))

        # Convert precipitation to mm/day
        pr = ds['pr'] * 86400

        # Resample to annual total precipitation (end of year)
        pr_annual = pr.resample(time='YE').sum(dim='time')

        # Calculate long-term mean annual rainfall
        pr_annual_mean = pr_annual.mean(dim='time')

        # Apply region mask as DataArray
        region_mask_da = region_mask.mask(pr_annual_mean)

        # Aggregate precipitation by region
        regional_means = pr_annual_mean.groupby(region_mask_da).mean()

        # Append to list
        bioregion_annual_data.append(regional_means)
        model_names.append(file.split(os.sep)[-3])  # Extract model name from path

        print(f"✅ Processed: {os.path.basename(file)}")

    except Exception as e:
        print(f"❌ Error processing {os.path.basename(file)}: {e}")
        continue

# ------------------------ Ensemble Mean and Output ------------------------ #

# Combine all models into one DataArray
ensemble_stack = xr.concat(bioregion_annual_data, dim='model')
ensemble_stack['model'] = model_names

# Compute ensemble mean per region
ensemble_mean_by_region = ensemble_stack.mean(dim='model')

# Convert to DataFrame
df = pd.DataFrame({
    "Bioregion": region_mask.names,
    "Mean_Annual_Rainfall_mm": ensemble_mean_by_region.values
}).sort_values(by="Mean_Annual_Rainfall_mm", ascending=False)

# Print results
print("\n📊 Mean Annual Rainfall by Bioregion (Ensemble Average):")
print(df)

# Merge mean rainfall values with bioregions GeoDataFrame
bioregion_mean_df = pd.DataFrame({
    "Veg_Biome": region_mask.names,
    "Mean_Annual_Rainfall_mm": ensemble_mean_by_region.values
})

bioregions_merged = bioregions.merge(bioregion_mean_df, on="Veg_Biome")

# ------------------- Plot ------------------- #
# Read CSV
towns_df = pd.read_csv(towns_csv_path, sep=';')  # or ',' if your file is comma-separated

# Strip spaces from column names
towns_df.columns = towns_df.columns.str.strip()

# Convert to GeoDataFrame
towns_gdf = gpd.GeoDataFrame(
    towns_df,
    geometry=gpd.points_from_xy(towns_df['lng'], towns_df['lat']),
    crs="EPSG:4326"
)

fig, ax = plt.subplots(figsize=(10, 8))

# Plot rainfall with steps
vmin = 0
vmax = bioregions_merged["Mean_Annual_Rainfall_mm"].max()
step = 100
ticks = range(int(vmin), int(vmax) + int(step), int(step))

bioregions_merged.plot(
    column="Mean_Annual_Rainfall_mm",
    cmap="Blues",
    linewidth=0.8,
    edgecolor="black",
    legend=True,
    legend_kwds={
        "label": "Mean Annual Rainfall (mm)",
        "orientation": "vertical",
        "ticks": ticks
    },
    ax=ax
)

# Plot towns
towns_gdf.plot(ax=ax, color='red', markersize=40, zorder=5)

# Add labels
for x, y, label in zip(towns_gdf.geometry.x, towns_gdf.geometry.y, towns_gdf['city']):
    ax.text(x, y, label, fontsize=9, ha='left', va='bottom')


# Finalize
ax.set_title("CMIP6 Mean Annual Rainfall by Vegetation Biome(1950-2014)", fontsize=14)
ax.set_axis_off()
plt.tight_layout()
plt.show()

