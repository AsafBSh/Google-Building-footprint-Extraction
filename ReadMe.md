# Google Buildings Footprint Processor

## Overview

This tool processes Google's building footprint data, allowing users to filter and extract building polygons from large CSV files or pre-divided geographical chunks. It supports operations such as dividing large datasets into manageable chunks and filtering buildings within specified geographical boundaries. 

## Features

Divide large CSV datasets into geographical chunks - Process building footprints from original CSV files or pre-divided chunks - Filter buildings within a specified bounding box - Support for both single-threaded and multi-threaded processing - Progress tracking for long-running operations - Output results in CSV or GeoJSON format.

## Acknowledgements

This tool processes building footprint data provided by Google's Open Buildings dataset.
[Google Open Buildings Footprint](https://sites.research.google/open-buildings/)
[Tiles geomemetry & URLs](https://openbuildings-public-dot-gweb-research.uw.r.appspot.com/public/tiles.geojson)

## Usage

Dividing a Large Dataset into Chunks To divide a large CSV file into geographical chunks:

```
python Main.py --divide <tile_id> [--override]
```

- `<tile_id>`: The ID of the tile to process (e.g., "145")
- `--override`: Optional flag to override existing chunk folders

Example:

```
python Main.py --divide 145 --override
```

This will create a folder named `<tile_id>_chunks` containing the divided GeoJSON files and a `chunk_boundaries.geojson` file.

### Filtering Buildings

To filter buildings within a specified bounding box in CRS:84 (EPSG:4326):

```
python Main.py -i <input> -o <output> --top-left <lat,lon> --bottom-right <lat,lon> [--fromdb]
```

- `-i, --input`: Input CSV file or chunk folder
- `-o, --output`: Output file (default: cropped_buildings.geojson)
- `--top-left`: Top-left coordinates of the bounding box (latitude,longitude)
- `--bottom-right`: Bottom-right coordinates of the bounding box (latitude,longitude)
- `--fromdb`: Optional flag to use the tiled database approach (pre-divided chunks)

Examples:

Processing from original CSV:

```
python Main.py -i 145_buildings.csv -o output.geojson --top-left 31.199039,27.621791 --bottom-right 31.165775,27.675779
```

* The bonding box can be extracted from any map source such as Google Maps/Earth or QGIS.



Processing from pre-divided chunks:

```
python Main.py -i 145_chunks --fromdb -o output.geojson --top-left 31.199039,27.621791 --bottom-right 31.165775,27.675779
```

## File Structure

- `Main.py`: The main script containing all the processing logic
- `tiles.geojson`: A GeoJSON file containing information about the tiles (required for dividing datasets)
- `<tile_id>_buildings.csv`: The original CSV file containing building footprints for a specific tile
- `<tile_id>_chunks/`: Folder containing divided chunks for a specific tile
  - `chunk_<n>.geojson`: Individual chunk files
  - `chunk_boundaries.geojson`: GeoJSON file containing the boundaries of all chunks

## How It Works

1. **Dividing Datasets**: 
   
   - Reads the original CSV file
   - Creates a grid of geographical chunks based on the tile's geometry
   - Filters buildings for each chunk and saves them as separate GeoJSON files

2. **Filtering Buildings**:
   
   - For CSV processing:
     - Reads the CSV file in chunks
     - Converts geometry strings to Shapely objects
     - Filters buildings that intersect with the specified bounding box
   - For pre-divided chunks:
     - Identifies relevant chunks based on the bounding box
     - Loads and processes only the relevant chunks
   - Combines filtered results and saves to the output file

## Performance Considerations

- The tool uses chunked reading for large CSV files to manage memory usage
- When processing from pre-divided chunks, only relevant chunks are loaded, improving performance for large datasets
- Multi-threading is used for CSV processing to speed up operations on multi-core systems

## Troubleshooting

- Ensure that the input CSV file has a 'geometry' column with valid WKT (Well-Known Text) strings
- Verify that the coordinates in the CSV file match the expected coordinate system (WGS84, EPSG:4326)
- For issues with chunk processing, check that the chunk files and `chunk_boundaries.geojson` are present in the chunk folder



## License

Copyright (c) <2024> <CC BY 4.0>

https://creativecommons.org/licenses/by/4.0/

## 
