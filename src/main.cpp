/*
Quickly generate 2.5D meshes from elevation models.
Copyright (C) 2018 Piero Toffanin

dem2mesh is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

dem2mesh is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with dem2mesh.  If not, see <https://www.gnu.org/licenses/>.
*/
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <unordered_map>
#include <omp.h>
#include "CmdLineParser.h"
#include "Logger.h"
#include "Simplify.h"
#include "tin.h"

#include "gdal_priv.h"
#include "cpl_conv.h" // for CPLMalloc()

Logger logWriter;

typedef struct Point2D{
   double x;
   double y;

   Point2D(): x(0.0), y(0.0){}
   Point2D(double x, double y): x(x), y(y){}
} Point2D;

typedef struct BoundingBox{
    Point2D max;
    Point2D min;

    BoundingBox(): max(Point2D()), min(Point2D()){}
    BoundingBox(Point2D min, Point2D max): max(max), min(min){}
} BoundingBox;

typedef struct PlyPoint{
    float x;
    float y;
    float z;
} PlyPoint;
size_t psize = sizeof(float) * 3;

typedef struct PlyFace{
    uint32_t p1;
    uint32_t p2;
    uint32_t p3;
} PlyFace;
size_t fsize = sizeof(uint32_t) * 3;

int arr_width, arr_height;
int subX, subY;
int blockSizeX, blockSizeY;

#define IS_BIG_ENDIAN (*(uint16_t *)"\0\xff" < 0x100)

cmdLineParameter< char* >
    InputFile( "inputFile" ) ,
    OutputFile( "outputFile" );
cmdLineParameter< int >
    MaxVertexCount( "maxVertexCount" ) ,
    MaxTileLength( "maxTileLength" ) ,
    Aggressiveness ( "aggressiveness" ) ,
    BandNum ( "bandNum" ) ,
    MaxConcurrency ( "maxConcurrency" );
cmdLineParameter< float >
    EdgeSwapThreshold ( "edgeSwapThreshold" );
cmdLineReadable
    Rtc ( "rtc" ),
    Verbose( "verbose" );

cmdLineReadable* params[] = {
    &InputFile , &OutputFile , 
    &MaxVertexCount , &MaxTileLength, &Aggressiveness, &BandNum, &MaxConcurrency,
    &EdgeSwapThreshold,
    &Rtc, &Verbose ,
    NULL
};

void help(char *ex){
    std::cout << "Usage: " << ex << std::endl
              << "\t -" << InputFile.name << " <input DSM raster>" << std::endl
              << "\t -" << OutputFile.name << " <output PLY mesh>" << std::endl
              << "\t [-" << MaxVertexCount.name << " <target number vertices> (Default: 100000)]" << std::endl
              << "\t [-" << MaxTileLength.name << " <max length of a tile. Smaller values take longer to process but reduce memory usage by splitting the meshing process into tiles.> (Default: 1000)]" << std::endl
              << "\t [-" << Aggressiveness.name << " <simplification aggressiveness factor. Higher values simplify the mesh more aggressively but can decrease the fidelity of the mesh. ([1-10] Default: 5)]" << std::endl
              << "\t [-" << BandNum.name << " <Band number> (Default: 1)]" << std::endl
              << "\t [-" << MaxConcurrency.name << " <threads> (Default: all cpus)]" << std::endl
              << "\t [-" << EdgeSwapThreshold.name << " <dot product threshold (between 0-1) to perform edge collapses. Performing edge collapses can violate the 2.5D constraints, but leads to better triangles for vertical structures> (Default: disabled)]" << std::endl
              << "\t [-" << Rtc.name << "]" << std::endl
              << "\t [-" << Verbose.name << "]" << std::endl;
    exit(EXIT_FAILURE);
}


void logArgs(cmdLineReadable* params[], Logger& logWriter){
    logWriter("Running with parameters:\n");
    char str[1024];
    for( int i=0 ; params[i] ; i++ ){
        if( params[i]->set ){
            params[i]->writeValue( str );
            if( strlen( str ) ) logWriter( "\t-%s %s\n" , params[i]->name , str );
            else                logWriter( "\t-%s\n" , params[i]->name );
        }
    }
}

Point2D geoLoc(float xloc, float yloc, double *affine){
    return Point2D(affine[0] + xloc*affine[1] + yloc*affine[2], affine[3] + xloc*affine[4] + yloc*affine[5]);
}

BoundingBox getExtent(GDALDataset *dataset){
    double affine[6];
    dataset->GetGeoTransform(affine);
    return BoundingBox(geoLoc(0, dataset->GetRasterYSize(), affine), geoLoc(dataset->GetRasterXSize(), 0, affine));
}

void writePly(const std::string &filename, int thread){
    // Start writing ply file
    std::ofstream f (filename, std::ios::binary);
    f << "ply" << std::endl;

    if (IS_BIG_ENDIAN){
      f << "format binary_big_endian 1.0" << std::endl;
    }else{
      f << "format binary_little_endian 1.0" << std::endl;
    }

    f   << "element vertex " << Simplify::vertices[thread]->size() << std::endl
        << "property float x" << std::endl
        << "property float y" << std::endl
        << "property float z" << std::endl
        << "element face " << Simplify::triangles[thread]->size() << std::endl
        << "property list uchar int vertex_indices" << std::endl
        << "end_header" << std::endl;

    PlyPoint p;
    for(Simplify::Vertex &v : *Simplify::vertices[thread]){
        p.x = static_cast<float>(v.p.x);
        p.y = static_cast<float>(v.p.y);
        p.z = static_cast<float>(v.p.z);
        f.write(reinterpret_cast<char *>(&p), psize);
    }

    PlyFace face;
    uint8_t three = 3;
    for(Simplify::Triangle &t : *Simplify::triangles[thread]){
        face.p1 = static_cast<uint32_t>(t.v[0]);
        face.p2 = static_cast<uint32_t>(t.v[1]);
        face.p3 = static_cast<uint32_t>(t.v[2]);

        f.write(reinterpret_cast<char *>(&three), sizeof(three));
        f.write(reinterpret_cast<char *>(&face), fsize);
    }

    f.close();
}

void writeBin(const std::string &filename, int blockWidth, int blockHeight, int thread){
    std::ofstream f (filename, std::ios::binary);

    unsigned long vsize = Simplify::vertices[thread]->size();
    unsigned long tsize = Simplify::triangles[thread]->size();

    f.write(reinterpret_cast<char *>(&blockWidth), sizeof(blockWidth));
    f.write(reinterpret_cast<char *>(&blockHeight), sizeof(blockHeight));
    f.write(reinterpret_cast<char *>(&vsize), sizeof(vsize));
    f.write(reinterpret_cast<char *>(&tsize), sizeof(tsize));

    PlyPoint p;
    for(Simplify::Vertex &v : *Simplify::vertices[thread]){
        p.x = static_cast<float>(v.p.x);
        p.y = static_cast<float>(v.p.y);
        p.z = static_cast<float>(v.p.z);
        f.write(reinterpret_cast<char *>(&p), psize);
    }

    PlyFace face;
    for(Simplify::Triangle &t : *Simplify::triangles[thread]){
        face.p1 = static_cast<uint32_t>(t.v[0]);
        face.p2 = static_cast<uint32_t>(t.v[1]);
        face.p3 = static_cast<uint32_t>(t.v[2]);
        f.write(reinterpret_cast<char *>(&face), fsize);
    }

    f.close();
}


void saveFinal(const std::string &filename, int thread, float edgeSwapThreshold){
    if (edgeSwapThreshold == -1.0){
        logWriter("Writing to file... ");
        writePly(filename, thread);
        logWriter(" done!\n");
    }else{
        logWriter("Performing edge collapses...\n");

        try{
            JMesh::init();
            Triangulation tin;

            unsigned long nv = Simplify::vertices[thread]->size();

            for(Simplify::Vertex &v : *Simplify::vertices[thread]){
                tin.V.appendTail(new Vertex(v.p.x,v.p.y,v.p.z));
            }

            Vertex *v;
            Node *n;
            unsigned long int i = 0;
            ExtVertex **var = (ExtVertex **)malloc(sizeof(ExtVertex *)*nv);

            TIN_FOREACHVERTEX(v, n) var[i++] = new ExtVertex(v);

            for(Simplify::Triangle &t : *Simplify::triangles[thread]){
                tin.CreateIndexedTriangle(var, t.v[0], t.v[1], t.v[2]);
            }

            for (i=0; i<nv; i++) delete(var[i]);
            free(var);

            tin.removeVertices();
            tin.cutAndStitch();
            tin.forceNormalConsistence();
            tin.duplicateNonManifoldVertices();
            tin.removeDuplicatedTriangles();

            tin.eulerUpdate();
            tin.verticalEdgeTagging(edgeSwapThreshold);
            tin.iterativeEdgeSwaps();
            tin.removeUnlinkedElements();
            tin.mergeCoincidentEdges();
            tin.removeDegenerateTriangles();

            logWriter("Writing to file... ");
            tin.savePLY(filename.c_str(), false);
            logWriter(" done!\n");
        }catch(const std::runtime_error &e){
            logWriter(e.what());
            logWriter("\nWriting to file... ");
            writePly(filename, thread);
            logWriter(" done!\n");
        }
    }
}

// Keep track of edge points's vertex IDs
// During the merge step we need to replace
// duplicate points with existing point's vertex IDs
// Half the array maps to X coords
// the other half maps to Y coords
//   1  2
// __|__|__ 3
// __|__|__ 4
//   |  |
//
// 1---2---3---4
long *pointToVertexIdMap = nullptr;

// (current vertex index) --> (existing vertex index) for duplicate points during merge
std::unordered_map<int, int> vertexToVertexMap;

void readBin(const std::string &filename, int blockX, int blockY, int thread){
    std::ifstream f (filename, std::ios::binary);

    int blockWidth, blockHeight;
    unsigned long vcount, tcount;
    unsigned long voffset = Simplify::vertices[thread]->size();
    float vertices[3];
    uint32_t triangles[3];

    f.read(reinterpret_cast<char *>(&blockWidth), sizeof(blockWidth));
    f.read(reinterpret_cast<char *>(&blockHeight), sizeof(blockHeight));
    f.read(reinterpret_cast<char *>(&vcount), sizeof(vcount));
    f.read(reinterpret_cast<char *>(&tcount), sizeof(tcount));

    int blockXPad = blockWidth - blockSizeX;
    int blockYPad = blockHeight - blockSizeY;
    int xOffset = blockX * blockSizeX - blockXPad;
    int yOffset = blockY * blockSizeY - blockYPad;

    int ptvidMapSize = arr_height * (subY - 1) + arr_width * (subX - 1) + 1;
    int ptvidYOffset = arr_height * (subX - 1) + 1;
    if (pointToVertexIdMap == nullptr){
        pointToVertexIdMap = new long[ptvidMapSize];
        memset(pointToVertexIdMap, -1, ptvidMapSize*sizeof(*pointToVertexIdMap));
    }

    int x, y;
    for (unsigned long i = 0; i < vcount; i++){
        f.read(reinterpret_cast<char *>(&vertices), sizeof(float) * 3);

        x = static_cast<int>(vertices[0]);
        y = static_cast<int>(vertices[1]);

        // Detect edge points
        if ((blockX > 0 && x == xOffset) || (blockX < subX - 1 && x == xOffset + blockWidth - 1)){
            if (pointToVertexIdMap[y + (x / blockSizeX) * arr_height] == -1){
                pointToVertexIdMap[y + (x / blockSizeX) * arr_height] = static_cast<long>(i + voffset);
            }else{
                // Already added from another block
                // Keep a reference to the point from the other block
                vertexToVertexMap[i + voffset] = pointToVertexIdMap[y + (x / blockSizeX) * arr_height];
            }
        }

        else if ((blockY > 0 && y == yOffset) || (blockY < subY - 1 && y == yOffset + blockHeight - 1)){
            if (pointToVertexIdMap[ptvidYOffset + x + (y / blockSizeY) * arr_width] == -1){
                pointToVertexIdMap[ptvidYOffset + x + (y / blockSizeY) * arr_width] = i + voffset;
            }else{
                // Already added from another block
                // Keep a reference to the point from the other block
                vertexToVertexMap[i + voffset] = pointToVertexIdMap[ptvidYOffset + x + (y / blockSizeY) * arr_width];
            }
        }

        Simplify::Vertex v;
        v.p.x = vertices[0];
        v.p.y = vertices[1];
        v.p.z = vertices[2];

        Simplify::vertices[thread]->push_back(v);
    }

    for (unsigned long i = 0; i < tcount; i++){
        f.read(reinterpret_cast<char *>(&triangles), sizeof(uint32_t) * 3);
        Simplify::Triangle t;
        loopk(0, 3) t.v[k] = triangles[k] + voffset;

        // Check vertices for substitutions
        loopk(0, 3) if (vertexToVertexMap.find(t.v[k]) != vertexToVertexMap.end()) t.v[k] = vertexToVertexMap[t.v[k]];

        // Skip degenerate triangles
        if (t.v[0] != t.v[1] && t.v[0] != t.v[2] && t.v[1] != t.v[2]){
            t.deleted = 0;
            Simplify::triangles[thread]->push_back(t);
        }
    }
}

void simplify(int target_count, int thread){
    unsigned long start_size = Simplify::triangles[thread]->size();
    if (target_count >= static_cast<int>(start_size)){
        logWriter("No simplification needed\n");
        Simplify::compact_mesh(thread);
        return;
    }

    Simplify::simplify_mesh(target_count, static_cast<double>(Aggressiveness.value), Verbose.set, thread);
    if ( Simplify::triangles[thread]->size() >= start_size) {
        logWriter("Unable to reduce mesh. We tried!\n");
    }
}

// From grid X/Y to world X/Y
void transform(const BoundingBox &extent, int thread){
    double ext_width = extent.max.x - extent.min.x;
    double ext_height = extent.max.y - extent.min.y;

    double center_x = 0.0;
    double center_y = 0.0;

    if (Rtc.set){
        center_x = (extent.min.x + extent.max.x) / 2.0;
        center_y = (extent.min.y + extent.max.y) / 2.0;
        logWriter("Relative to center (RTC) coordinate set to (%f, %f)\n", center_x, center_y);
    }

    for(Simplify::Vertex &v : *Simplify::vertices[thread]){
        v.p.x = extent.min.x + (static_cast<double>(v.p.x) / static_cast<double>(arr_width)) * ext_width - center_x;
        v.p.y = extent.max.y - (static_cast<double>(v.p.y) / static_cast<double>(arr_height)) * ext_height - center_y;
    }
}

BoundingBox extent;

int main(int argc, char **argv) {
    cmdLineParse( argc-1 , &argv[1] , params );
    if( !InputFile.set || !OutputFile.set ) help(argv[0]);
    if ( !MaxVertexCount.set ) MaxVertexCount.value = 100000;
    if ( !MaxTileLength.set ) MaxTileLength.value = 1000;
    if ( !Aggressiveness.set ) Aggressiveness.value = 5;
    if ( !BandNum.set ) BandNum.value = 1;
    if ( !MaxConcurrency.set ) MaxConcurrency.value = omp_get_max_threads();
    if ( MaxConcurrency.value == 0 ) MaxConcurrency.value = 1;
    if ( !EdgeSwapThreshold.set ) EdgeSwapThreshold.value = -1.0;

    Aggressiveness.value = std::min(10, std::max(1, Aggressiveness.value));
    logWriter.verbose = Verbose.set;
    logWriter.outputFile = NULL;
    logArgs(params, logWriter);

    GDALDataset  *dataset;
    GDALAllRegister();
    dataset = (GDALDataset *) GDALOpen( InputFile.value, GA_ReadOnly );
    if( dataset != NULL )
    {
        arr_width = dataset->GetRasterXSize();
        arr_height = dataset->GetRasterYSize();

        logWriter("Raster Size is %dx%d\n", arr_width, arr_height);
        extent = getExtent(dataset);

        logWriter("Extent is (%f, %f), (%f, %f)\n", extent.min.x, extent.max.x, extent.min.y, extent.max.y);

        GDALRasterBand *band = dataset->GetRasterBand(BandNum.value);

        int hasNoData = FALSE;
        double nodata = band->GetNoDataValue(&hasNoData);

        if (hasNoData){
            logWriter("NoData value: %.18g\n", nodata);
        }
        logWriter("Description: %s\n", band->GetDescription());

        unsigned long long int vertex_count = static_cast<unsigned long long int>(arr_height) *
                                              static_cast<unsigned long long int>(arr_width);

        logWriter("Reading raster...\n");
        logWriter("Total vertices before simplification: %llu\n", vertex_count);

        int qtreeLevelsX = 0;
        while(true){
            subX = (int)pow(2, qtreeLevelsX);
            blockSizeX = arr_width / subX;
            if (blockSizeX > MaxTileLength.value){
                qtreeLevelsX++;
            }else{
                break;
            }
        }
        int qtreeLevelsY = 0;
        while(true){
            subY = (int)pow(2, qtreeLevelsY);
            blockSizeY = arr_height / subY;
            if (blockSizeY > MaxTileLength.value){
                qtreeLevelsY++;
            }else{
                break;
            }
        }

        int numBlocks = subX * subY;

        logWriter("Subdivisions: %dx%d\n", subX, subY);
        logWriter("Splitting area in %d\n", numBlocks);
        logWriter("Block size is %d, %d\n", blockSizeX, blockSizeY);

        logWriter("Concurrency set to %d\n", MaxConcurrency.value);
        Simplify::allocate(MaxConcurrency.value);
        omp_set_num_threads(MaxConcurrency.value);

        int rasterDataBlocks = std::min(MaxConcurrency.value, numBlocks);
        logWriter("Allocating %d raster data blocks of %d bytes\n", rasterDataBlocks, sizeof(float) * (blockSizeX + 1));
        float *rasterData = new float[rasterDataBlocks * (blockSizeX + 1)];

        omp_lock_t readLock;
        omp_init_lock(&readLock);

        // If there's a nodata value, we need to perform a first pass
        // on the raster by checking for nodata values, as some blocks
        // might be completely empty. If we have lots of empty blocks,
        // our estimate for the triangle count will be terribly off
        int emptyBlocks = 0;
        bool *emptyTable = new bool[subX * subY];
        memset(emptyTable, false, subX * subY);

        if (hasNoData){
            #pragma omp parallel for collapse(2)
            for (int blockX = 0; blockX < subX; blockX++){
                for (int blockY = 0; blockY < subY; blockY++){
                    int t = omp_get_thread_num();
                    int blockXPad = blockX == 0 ? 0 : 1; // Blocks > 0 need to re-add a column for seamless meshing
                    int blockYPad = blockY == 0 ? 0 : 1; // Blocks > 0 need to re-add a row for seamless meshing
                    int xOffset = blockX * blockSizeX - blockXPad;
                    int yOffset = blockY * blockSizeY - blockYPad;

                    logWriter("Scanning block (%d,%d)\n", blockX, blockY);

                    bool empty = true;

                    omp_set_lock(&readLock);
                    for (int y = 0; y < blockSizeY + blockYPad; y++){
                        if (band->RasterIO( GF_Read, xOffset, yOffset + y, blockSizeX + blockXPad, 1,
                                            rasterData + t * (blockSizeX + 1), blockSizeX + blockXPad, 1, GDT_Float32, 0, 0 ) == CE_Failure){
                            std::cerr << "Cannot access raster data" << std::endl;
                            exit(EXIT_FAILURE);
                        }
                        for (int x = 0; x < blockSizeX + blockXPad; x++){
                            float z = (rasterData + t * (blockSizeX + 1))[x];
                            if (z != nodata){
                                empty = false;
                                y = blockSizeY + blockYPad; // Exit outer loop
                                break;
                            }
                        }
                    }
                    omp_unset_lock(&readLock);

                    if (empty){
                        #pragma omp atomic
                        emptyBlocks++;
                        emptyTable[blockY * subX + blockX] = true;
                    }
                }
            }
        }

        logWriter("Empty blocks: %d\n", emptyBlocks);
        numBlocks -= emptyBlocks;

        #pragma omp parallel for collapse(2)
        for (int blockX = 0; blockX < subX; blockX++){
            for (int blockY = 0; blockY < subY; blockY++){
                if (emptyTable[blockY * subX + blockX]) continue;

                int t = omp_get_thread_num();
                int blockXPad = blockX == 0 ? 0 : 1; // Blocks > 0 need to re-add a column for seamless meshing
                int blockYPad = blockY == 0 ? 0 : 1; // Blocks > 0 need to re-add a row for seamless meshing
                int xOffset = blockX * blockSizeX - blockXPad;
                int yOffset = blockY * blockSizeY - blockYPad;

                Simplify::vertices[t]->clear();
                Simplify::triangles[t]->clear();

                logWriter("Processing block (%d,%d)\n", blockX, blockY);

                omp_set_lock(&readLock);
                for (int y = 0; y < blockSizeY + blockYPad; y++){
                    if (band->RasterIO( GF_Read, xOffset, yOffset + y, blockSizeX + blockXPad, 1,
                                        rasterData + t * (blockSizeX + 1), blockSizeX + blockXPad, 1, GDT_Float32, 0, 0 ) == CE_Failure){
                        std::cerr << "Cannot access raster data" << std::endl;
                        exit(EXIT_FAILURE);
                    }
                    for (int x = 0; x < blockSizeX + blockXPad; x++){
                        Simplify::Vertex v;
                        v.p.x = xOffset + x;
                        v.p.y = yOffset + y;
                        v.p.z = (rasterData + t * (blockSizeX + 1))[x];

                        Simplify::vertices[t]->push_back(v);
                    }
                }
                omp_unset_lock(&readLock);

                unsigned int cols = blockSizeX + blockXPad;
                unsigned int rows = blockSizeY + blockYPad;

                for (unsigned int y = 0; y < rows - 1; y++){
                    for (unsigned int x = 0; x < cols - 1; x++){
                        Simplify::Triangle t1;
                        t1.v[0] = cols * (y + 1) + x;
                        t1.v[1] = cols * y + x + 1;
                        t1.v[2] = cols * y + x;
                        if (y == 0 || x == 0 || y == rows - 2 || x == cols - 2) t1.deleted = -1; // freeze
                        else t1.deleted = 0;

                        if (!hasNoData ||
                           ((*Simplify::vertices[t])[t1.v[0]].p.z != nodata &&
                            (*Simplify::vertices[t])[t1.v[1]].p.z != nodata &&
                            (*Simplify::vertices[t])[t1.v[2]].p.z != nodata)){
                            Simplify::triangles[t]->push_back(t1);
                        }

                        Simplify::Triangle t2;
                        t2.v[0] = cols * (y + 1) + x;
                        t2.v[1] = cols * (y + 1) + x + 1;
                        t2.v[2] = cols * y + x + 1;
                        if (y == 0 || x == 0 || y == rows - 2 || x == cols - 2) t2.deleted = -1; // freeze
                        else t2.deleted = 0;

                        if (!hasNoData ||
                           ((*Simplify::vertices[t])[t2.v[0]].p.z != nodata &&
                            (*Simplify::vertices[t])[t2.v[1]].p.z != nodata &&
                            (*Simplify::vertices[t])[t2.v[2]].p.z != nodata)){
                            Simplify::triangles[t]->push_back(t2);
                        }

                    }
                }

                int trianglesPerBlock = (MaxVertexCount.value * 2) / numBlocks;

                // If we have a merge step,
                // overshoot the triangle count requirement
                // since we'll simplify the final mesh anyway.
                // This leads to more uniform meshes.
                if (numBlocks > 1) trianglesPerBlock = trianglesPerBlock * 3 / 2;

                int target_count = std::min(trianglesPerBlock, static_cast<int>(Simplify::triangles[t]->size()));

                logWriter("Sampled %d faces, target is %d\n", static_cast<int>(Simplify::triangles[t]->size()), target_count);
                logWriter("Simplifying...\n");
                simplify(target_count, t);

                if (numBlocks <= 1){
                    transform(extent, t);
                    logWriter("Single block, saving to PLY\n");
                    saveFinal(OutputFile.value, t, EdgeSwapThreshold.value);
                }else{
                    logWriter("Writing to binary file...");
                    std::stringstream ss;
                    ss << OutputFile.value << "." << blockX << "-" << blockY << ".bin";
                    writeBin(ss.str(), blockSizeX + blockXPad, blockSizeY + blockYPad, t);
                    logWriter(" done!\n");
                }
            }
        }

        delete[] rasterData;
        GDALClose(dataset);

        if (numBlocks > 1){
            // Merge
            logWriter("Merge step...\n");

            Simplify::vertices[0]->clear();
            Simplify::triangles[0]->clear();

            for (int blockX = 0; blockX < subX; blockX++){
                for (int blockY = 0; blockY < subY; blockY++){
                    if (emptyTable[blockY * subX + blockX]) continue;

                    std::stringstream ss;
                    ss << OutputFile.value << "." << blockX << "-" << blockY << ".bin";
                    logWriter("Reading %s\n", ss.str().c_str());
                    readBin(ss.str(), blockX, blockY, 0);

                    if (std::remove(ss.str().c_str()) != 0){
                        logWriter("Error while deleting intermediate file: %s\n", ss.str().c_str());
                    }
                }
            }

            // Cleanup
            if (pointToVertexIdMap != nullptr){
                delete[] pointToVertexIdMap;
                pointToVertexIdMap = nullptr;
            }
            vertexToVertexMap.clear();

            logWriter("Simplifying final mesh...\n");
            int target_count = std::min(MaxVertexCount.value * 2, static_cast<int>(Simplify::triangles[0]->size()));
            simplify(target_count, 0);
            transform(extent, 0);
            saveFinal(OutputFile.value, 0, EdgeSwapThreshold.value);
        }

        Simplify::cleanup();
    }else{
        std::cerr << "Cannot open " << InputFile.value << std::endl;
    }
}
