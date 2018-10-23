# dem2mesh

Fast generation of 2.5D meshes from elevation models.

![image](https://user-images.githubusercontent.com/1951843/47350997-15d3da00-d685-11e8-8d9f-e394fc17859e.png)

![image](https://user-images.githubusercontent.com/1951843/47351205-7e22bb80-d685-11e8-87c5-33b21ae05b75.png)

## Dependencies

GDAL is the only dependency. To install it run:

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

:warning: The DEM should use a coordinate reference system (CRS) that has matching horizontal and vertical units (for example, UTM). If you are getting stretched or deformed results, double check your CRS and use `gdalwarp` to reproject the raster DEM. https://www.gdal.org/gdalwarp.html

## Options

| Flag | Description | Required |
| --- | --- | --- |
| -inputFile | Path to DEM raster datasource (GDAL compatible) | ✔ |
| -outputFile | Path to PLY mesh output | ✔ |
| -maxVertexCount | Target number of vertices for the output mesh. This number is not always guaranteed to be reached. Defaults to `100000` | |
| -maxTileLength | Max length of a tile. Smaller values take longer to process but reduce memory usage by splitting the meshing process into tiles. Defaults to `1000`. | |
| -bandNum | Raster band # to use. Defaults to `1`. | |
| -rtc | Use Relative To Center (RTC) X/Y coordinates in the output PLY mesh. This can be useful since most 3D visualization software use floating coordinate precision to represent vertices and using absolute coordinates might lead to jittering artifacts. | |
| -verbose | Print verbose output. | |
