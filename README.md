# dem2mesh

Fast generation of 2.5D meshes from elevation models.

## Dependencies

GDAL is the only requirement. To install it run:

```
sudo apt-get install -y libgdal-dev 
```

## Building

```bash
mkdir build && cd build
cmake ..
make
``` 

## Running

```bash
./dem2mesh -rtc -verbose -inputFile dem.tif -outputFile mesh.ply
```

## Options

| Flag | Description | Required |
| --- | --- | --- |
| -inputFile | Path to DEM raster datasource (GDAL compatible) | ✔ |
| -outputFile | Path to PLY mesh output | ✔ |
| -maxVertexCount | Target number of vertices for the output mesh. This number is not always guaranteed to be reached. Defaults to `100000` | |
| -maxTileLength | Max length of a tile. Smaller values take longer to process but reduce memory usage by splitting the meshing process into tiles. Defaults to `1000`. | |
| -rtc | Use Relative To Center (RTC) X/Y coordinates in the output PLY mesh. This can be useful since most 3D visualization software use floating coordinate precision to represent vertices and using absolute coordinates might lead to jittering artifacts. | |
| -verbose | Print verbose output. | |
