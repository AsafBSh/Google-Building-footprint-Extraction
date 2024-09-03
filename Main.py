import pandas as pd
import geopandas as gpd
import argparse
import re
from rtree import index
import multiprocessing as mp
from shapely.geometry import Polygon, MultiPolygon, box
from shapely import wkt
from tqdm import tqdm
import json
import os
import shutil

def create_spatial_index(gdf):
    """Create a spatial index for faster geometric operations."""
    idx = index.Index()
    for id, geometry in enumerate(gdf['geometry']):
        if geometry.is_empty:
            continue
        idx.insert(id, geometry.bounds)
    return idx

def process_chunk(args):
    """Process a chunk of the DataFrame, filtering polygons that intersect with the bounding box."""
    chunk, bb, spatial_index = args
    # Use spatial index to quickly find potential matches
    potential_matches_idxs = list(spatial_index.intersection(bb.bounds))
    filtered_chunk = chunk.iloc[potential_matches_idxs]
    # Return only geometries that actually intersect with the bounding box
    return filtered_chunk[filtered_chunk.geometry.intersects(bb)]

def parse_polygon(polygon_str):
    """Parse a POLYGON or MULTIPOLYGON string into a Shapely geometry object."""
    if isinstance(polygon_str, str):
        try:
            return wkt.loads(polygon_str)
        except:
            coords_str = re.findall(r'\(\((.*?)\)\)', polygon_str)
            if coords_str:
                all_coords = [
                    [tuple(map(float, pair.split())) for pair in coord.split(',')]
                    for coord in coords_str
                ]
                if len(all_coords) == 1:
                    return Polygon(all_coords[0])
                else:
                    return MultiPolygon([Polygon(coords) for coords in all_coords])
    return polygon_str

def load_and_filter_polygons(input_file, output_file, top_left, bottom_right, chunksize=100000, use_chunks=False):
    """
    Load polygons from chunked GeoJSON files or CSV, filter them based on a bounding box,
    and save the filtered polygons to a new GeoJSON file.
    """
    bb = box(min(top_left[1], bottom_right[1]),
             min(top_left[0], bottom_right[0]),
             max(top_left[1], bottom_right[1]),
             max(top_left[0], bottom_right[0]))

    filtered_polygons = []
    total_processed = 0

    if use_chunks:
        # Working with geographical chunks (GeoJSON files)
        chunk_boundaries = gpd.read_file(os.path.join(input_file, "chunk_boundaries.geojson"))
        relevant_chunks = chunk_boundaries[chunk_boundaries.intersects(bb)]

        with tqdm(total=len(relevant_chunks), desc="Processing chunks") as pbar:
            for _, chunk in relevant_chunks.iterrows():
                chunk_file = os.path.join(input_file, f"chunk_{chunk['chunk_id']}.geojson")
                gdf_chunk = gpd.read_file(chunk_file)

                # Filter polygons within the chunk
                filtered_chunk = gdf_chunk[gdf_chunk.intersects(bb)]
                filtered_polygons.append(filtered_chunk)

                pbar.update(1)

    else:
        # Working with original CSV file
        chunks = pd.read_csv(input_file, chunksize=chunksize, dtype={'geometry': str})
        with tqdm(total=None, desc="Processing CSV chunks") as pbar:
            for chunk in chunks:

                # Convert string representation of polygons to Shapely geometries
                chunk['geometry'] = chunk['geometry'].apply(parse_polygon)
                gdf_chunk = gpd.GeoDataFrame(chunk, geometry='geometry')
                # Drop rows with invalid geometries
                gdf_chunk = gdf_chunk.dropna(subset=['geometry'])

                # Verify CRS
                if gdf_chunk.crs is None:
                    gdf_chunk.set_crs(epsg=4326, inplace=True)  # Assuming WGS84, adjust if different

                # Filter polygons within the chunk
                filtered_chunk = gdf_chunk[gdf_chunk.intersects(bb)]
                filtered_polygons.append(filtered_chunk)
                pbar.update(len(chunk))

    # Combine all filtered chunks

    if filtered_polygons:
        final_result = gpd.GeoDataFrame(pd.concat(filtered_polygons, ignore_index=True))

        # Simplify geometries (optional, may reduce precision)
        final_result['geometry'] = final_result['geometry'].simplify(tolerance=0.0001)

        # Write to file (CSV or GeoJSON based on output_file extension)
        if output_file:
            if output_file.lower().endswith('.geojson'):
                final_result.to_file(output_file, driver='GeoJSON')
            else:
                final_result.to_csv(output_file, index=False)

        return final_result
    else:
        print("No polygons found within the specified bounding box.")
        return gpd.GeoDataFrame()

def parse_coordinates(coord_str):
    """Parse a string of coordinates into a tuple of floats."""
    try:
        lat, lon = map(float, coord_str.split(','))
        return (lat, lon)
    except ValueError:
        raise argparse.ArgumentTypeError("Coordinates must be in the format 'latitude,longitude'")

def load_tiles_geojson(file_path='tiles.geojson'):
    """ Load Google's tiles file to understand the regions of the relevant CSV data"""
    with open(file_path, 'r') as f:
        data = json.load(f)
    features = data['features']
    tiles = []
    for feature in features:
        tile = {
            'tile_id': feature['properties']['tile_id'],
            'tile_url': feature['properties']['tile_url'],
            'size_mb': feature['properties']['size_mb'],
            'geometry': Polygon(feature['geometry']['coordinates'][0])
        }
        tiles.append(tile)
    return tiles

def create_geographic_chunks(tile_geometry, num_chunks=1000):
    """Create geographic chunks within the given tile geometry."""
    minx, miny, maxx, maxy = tile_geometry.bounds
    dx = (maxx - minx) / int(num_chunks ** 0.5)
    dy = (maxy - miny) / int(num_chunks ** 0.5)

    chunks = []
    for i in range(int(num_chunks ** 0.5)):
        for j in range(int(num_chunks ** 0.5)):
            chunk_box = box(minx + i * dx, miny + j * dy, minx + (i + 1) * dx, miny + (j + 1) * dy)
            if chunk_box.intersects(tile_geometry):
                chunks.append(chunk_box)

    return chunks

def divide_database(tile_id, override=False, num_chunks=1000):
    """ Divide the CSV buildings footprint database into smaller pieces"""
    tiles = load_tiles_geojson()
    tile = next((t for t in tiles if t['tile_id'] == tile_id), None)

    if not tile:
        print(f"Error: Tile {tile_id} not found in tiles.geojson.")
        return

    file_name = f"{tile_id}_buildings.csv"
    output_folder = f"{tile_id}_chunks"

    if not os.path.exists(file_name):
        print(f"Error: File {file_name} not found.")
        return

    if os.path.exists(output_folder):
        if not override:
            user_input = input(f"Folder {output_folder} already exists. Do you want to override it? (y/n): ")
            if user_input.lower() != 'y':
                print("Operation cancelled.")
                return
        shutil.rmtree(output_folder)

    os.makedirs(output_folder)

    print(f"Dividing {file_name} into geographic chunks. This may take a while for large files...")

    # Create geographic chunks
    chunks = create_geographic_chunks(tile['geometry'], num_chunks)

    # Read the CSV file into a GeoDataFrame
    df = pd.read_csv(file_name)
    gdf = gpd.GeoDataFrame(df, geometry=gpd.GeoSeries.from_wkt(df['geometry']))

    # Set the CRS if it's not already set
    if gdf.crs is None:
        gdf.set_crs(epsg=4326, inplace=True)  # Assuming WGS84, adjust if different

    # Process each chunk
    for i, chunk in enumerate(tqdm(chunks, desc="Processing chunks")):
        chunk_gdf = gdf[gdf.intersects(chunk)]
        if not chunk_gdf.empty:
            chunk_gdf.to_file(f"{output_folder}/chunk_{i}.geojson", driver='GeoJSON')

    # Save chunk boundaries as a single GeoJSON file
    chunk_boundaries = gpd.GeoDataFrame(geometry=chunks)
    chunk_boundaries['chunk_id'] = range(len(chunks))
    chunk_boundaries.to_file(f"{output_folder}/chunk_boundaries.geojson", driver="GeoJSON")

    print(f"Database divided into chunks in folder: {output_folder}")

def main():
    parser = argparse.ArgumentParser(description="Filter polygons from a CSV file or chunked GeoJSON files based on a bounding box.")
    parser.add_argument('-i', '--input', help="Input CSV file or chunk folder")
    parser.add_argument('-o', '--output', default='cropped_buildings.geojson', help="Output file (default: cropped_buildings.geojson)")
    parser.add_argument('--top-left', type=parse_coordinates, help="Top-left coordinates of the bounding box (latitude,longitude)")
    parser.add_argument('--bottom-right', type=parse_coordinates, help="Bottom-right coordinates of the bounding box (latitude,longitude)")
    parser.add_argument('--fromdb', action='store_true', help="Use the tiled database approach")
    parser.add_argument('--divide', type=str, help="Divide the specified tile_id into chunks")
    parser.add_argument('--override', action='store_true', help="Override existing chunk folders")

    args = parser.parse_args()

    if args.divide:
        divide_database(args.divide, override=args.override)
        return

    # Check if required arguments are provided for filtering operation
    if not all([args.input, args.output, args.top_left, args.bottom_right]):
        parser.error("For filtering operations, the following arguments are required: -i/--input, -o/--output, --top-left, --bottom-right")

    use_chunks = args.fromdb

    if use_chunks:
        if not os.path.isdir(args.input):
            print(f"Error: {args.input} is not a directory. When using --fromdb, input should be a chunk folder.")
            return
    else:
        if not os.path.isfile(args.input):
            print(f"Error: {args.input} is not a file. When not using --fromdb, input should be a CSV file.")
            return

    filtered_polygons = load_and_filter_polygons(
        args.input, args.output, args.top_left, args.bottom_right, use_chunks=use_chunks
    )

    num_filtered = len(filtered_polygons)

    if num_filtered > 0:
        print(f"Operation completed successfully. {num_filtered} polygons were cropped and exported to {args.output}")
    else:
        print(f"Operation completed, but no polygons were found within the specified bounding box. The output file {args.output} is empty.")

if __name__ == "__main__":
    main()
