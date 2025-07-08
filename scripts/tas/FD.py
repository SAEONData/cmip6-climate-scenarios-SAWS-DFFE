import os
import glob
import pandas as pd
import xarray as xr
import geopandas as gpd
import regionmask
import matplotlib.pyplot as plt
import xclim
from xclim.indices import frost_days

# ------------------ Config ------------------ #
data_path = "/content/drive/MyDrive/Climate Data TTT"
shapefile_path = "/content/drive/MyDrive/CLIMREG OVERLAPS REMOVED/cleaned_veg_biome_clim_reg.shp"
towns_csv_path = "/content/drive/MyDrive/CLIMREG OVERLAPS REMOVED/cities.csv"

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

# ------------------ Find Files ------------------ #
all_hist_files = sorted(glob.glob(os.path.join(data_path, "*", "tasmin_day_*historical*.nc")))
print(f"📂 Found {len(all_hist_files)} historical NetCDF files.")

# ------------------ Process Files ------------------ #
bioregion_annual_data = []
model_names = []

for file in all_hist_files:
    try:
        ds = xr.open_dataset(file)
        ds = ds.sel(lat=slice(*lat_bounds), lon=slice(*lon_bounds))

        if 'tasmin' not in ds:
            print(f"⚠️ Skipping {file}, 'tasmin' not found.")
            continue

        tasmin = ds['tasmin'] - 273.15  # Convert to °C
        tasmin.attrs['units'] = 'degC'
        tasmin.attrs['standard_name'] = 'air_temperature'

        frost_per_year = frost_days(tasmin, freq='YS')
        frost_mean = frost_per_year.mean(dim='time')

        # Apply region mask
        mask = region_mask.mask(frost_mean)
        regional_means = frost_mean.groupby(mask).mean()

        bioregion_annual_data.append(regional_means)
        model_names.append(file.split(os.sep)[-2])

        print(f"✅ Processed: {os.path.basename(file)}")

    except Exception as e:
        print(f"❌ Error processing {os.path.basename(file)}: {e}")
        continue

# ------------------ Ensemble Mean ------------------ #
ensemble_stack = xr.concat(bioregion_annual_data, dim='model')
ensemble_stack['model'] = model_names
ensemble_mean_by_region = ensemble_stack.mean(dim='model')

# ------------------ Tabular Output ------------------ #
df = pd.DataFrame({
    "Bioregion": region_mask.names,
    "Frost_Days": ensemble_mean_by_region.values
}).sort_values(by="Frost_Days", ascending=False)

print("\n📊 Frost Days (TN < 0°C) by Bioregion (Ensemble Average):")
print(df)

# ------------------ Spatial Output ------------------ #
bioregion_mean_df = pd.DataFrame({
    "Veg_Biome": region_mask.names,
    "Frost_Days": ensemble_mean_by_region.values
})
bioregions_merged = bioregions.merge(bioregion_mean_df, on="Veg_Biome")

towns_df = pd.read_csv(towns_csv_path, sep=';')
towns_df.columns = towns_df.columns.str.strip()

towns_gdf = gpd.GeoDataFrame(
    towns_df,
    geometry=gpd.points_from_xy(towns_df['lng'], towns_df['lat']),
    crs="EPSG:4326"
)

# ------------------ Plot Map ------------------ #
fig, ax = plt.subplots(figsize=(10, 8))

vmin = 0
vmax = bioregions_merged["Frost_Days"].max()
step = 10
ticks = range(int(vmin), int(vmax) + step, step)

bioregions_merged.plot(
    column="Frost_Days",
    cmap="Blues",
    linewidth=0.8,
    edgecolor="black",
    legend=True,
    legend_kwds={
        "label": "Mean Annual Frost Days (TN < 0°C)",
        "orientation": "vertical",
        "ticks": ticks
    },
    ax=ax
)

towns_gdf.plot(ax=ax, color='red', markersize=40, zorder=5)
for x, y, label in zip(towns_gdf.geometry.x, towns_gdf.geometry.y, towns_gdf['city']):
    ax.text(x, y, label, fontsize=9, ha='left', va='bottom')

ax.set_title("CMIP6: FD, Number of days when TN < 0°C by Vegetation Biomes (1950–2014)", fontsize=14)
ax.set_axis_off()
plt.tight_layout()
plt.show()
