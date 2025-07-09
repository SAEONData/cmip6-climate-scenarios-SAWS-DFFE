import os
import glob
import pandas as pd
import xarray as xr
import geopandas as gpd
import regionmask
import matplotlib.pyplot as plt
from xclim.indices import tx_min  # ✅ For Coldest Day (TXn)

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
all_hist_files = sorted(glob.glob(os.path.join(data_path, "*", "tasmax_day_*historical*.nc")))
print(f"📂 Found {len(all_hist_files)} historical NetCDF files.")

# ------------------ Process Files ------------------ #
bioregion_txn_data = []
model_names = []

for file in all_hist_files:
    try:
        ds = xr.open_dataset(file)
        ds = ds.sel(lat=slice(*lat_bounds), lon=slice(*lon_bounds))

        if 'tasmax' not in ds:
            print(f"⚠️ Skipping {file}, 'tasmax' not found.")
            continue

        tasmax = ds['tasmax'] - 273.15  # Convert from K to °C
        tasmax.attrs['units'] = 'degC'
        tasmax.attrs['standard_name'] = 'air_temperature'

        # ✅ Coldest Day (TXn)
        txn = tx_min(tasmax, freq='YS')
        txn_mean = txn.mean(dim='time')

        # Apply region mask
        mask = region_mask.mask(txn_mean)
        regional_means = txn_mean.groupby(mask).mean()

        bioregion_txn_data.append(regional_means)
        model_names.append(file.split(os.sep)[-2])

        print(f"✅ Processed: {os.path.basename(file)}")

    except Exception as e:
        print(f"❌ Error processing {os.path.basename(file)}: {e}")
        continue

# ------------------ Ensemble Mean ------------------ #
ensemble_stack = xr.concat(bioregion_txn_data, dim='model')
ensemble_stack['model'] = model_names
ensemble_mean_by_region = ensemble_stack.mean(dim='model')

# ------------------ Tabular Output ------------------ #
df = pd.DataFrame({
    "Bioregion": region_mask.names,
    "Min_TX": ensemble_mean_by_region.values
}).sort_values(by="Min_TX")

print("\n📊 Coldest Daily Maximum Temperature (TXn, °C) by Bioregion (Ensemble Average):")
print(df)

# ------------------ Spatial Output ------------------ #
bioregion_mean_df = pd.DataFrame({
    "Veg_Biome": region_mask.names,
    "Min_TX": ensemble_mean_by_region.values
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

vmin = bioregions_merged["Min_TX"].min()
vmax = bioregions_merged["Min_TX"].max()
step = 2
ticks = range(int(vmin), int(vmax) + step, step)

bioregions_merged.plot(
    column="Min_TX",
    cmap="coolwarm",
    linewidth=0.8,
    edgecolor="black",
    legend=True,
    legend_kwds={
        "label": "Coldest Daily Maximum Temperature (TXn, °C)",
        "orientation": "vertical",
        "ticks": ticks
    },
    ax=ax
)

towns_gdf.plot(ax=ax, color='black', markersize=40, zorder=5)
for x, y, label in zip(towns_gdf.geometry.x, towns_gdf.geometry.y, towns_gdf['city']):
    ax.text(x, y, label, fontsize=9, ha='left', va='bottom')

ax.set_title("CMIP6: Coldest Day of the Year (TXn, °C) by Bioregion (1950–2014)", fontsize=14)
ax.set_axis_off()
plt.tight_layout()
plt.show()
