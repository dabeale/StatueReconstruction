/* ----------------------------------------------------------------------
 * Copyright (C) 2016 Daniel Beale. All rights reserved.
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 * ---------------------------------------------------------------------- */

#include <fstream>
#include <list>
#include <vector>

#include "Image.h"
#include "Message.h"
#include "Carving.h"
#include "Matrix.h"

/**
 * @brief readNameFile
 * Read the a collection of filenames from a newline delimited file.
 * @param filename The path to the file
 * @return A list of filenames
 */
std::list< std::string > readNameFile( const std::string& filename)
{
    std::fstream fs (filename, std::fstream::in);
    std::list< std::string > outputlist;
    std::string line;
    while( !fs.eof() )
    {
        std::getline(fs, line);
        if(line.size() > 0)
            outputlist.push_back(line);
    }

    fs.close();
    return outputlist;
}

/**
 * @brief print_usage
 * Print the usage of the software
 * @param msg An interface to the output stream.
 */
void print_usage( Stream::Message& msg)
{
    msg.printerr( " Usage:" );
    msg.printerr( "    ./Carving cameras.txt images.txt pointcloud.ply outputfolder (bundlerformat=1 dopcloudoptimisation=0)" );
}

/**
 * @brief Carving
 * This executable is for computing a collection of space carving algorithms given a camera
 * collection, a point cloud and corresponding segmentations. There is documentation in the
 * Carve namespace which details how each of the functions work.
 *
 * The minimum input requirements are,
 *  1. A textfile containing the locations of the camera files ( either Bundler or Octave format)
 *     Each filename must appear on a seperate line, with no spaces.
 *  2. A textfile containing the location of the segmentations. The image files must be readable
 *     by opencv, which covers most file types, as long as they have been compiled with opencv.
 *     Files such as raw and hdr images may not work properly. There must be the same number of
 *     elements as in the camera text file, and the must correspond to one another.
 *  3. A pointcloud, which must be in the stanford ply format. This point cloud does not need any other
 *     attributes than the vertices.
 *  4. An output folder to store the output meshes. The software will output a maximum of three meshes
 *     depending on the input arguments.
 *  5. (optional) A value 0 or 1 indicating whether to use the octave or Bundler matrix formats, respectively.
 *     This value is one by default.
 *  6. (optional) A value 0 or 1 to indicate whether do do a Hoppe / Visual Hull blend, as documented in
 *     the Carve namespace.
 *
 * The software will output,
 *  - Hull.ply which is the visual hull of the object
 *  - PCloud.ply which is the original point cloud, with all of the points outside of the visual hull removed
 *    and the normals transfered from the visual hull. This allows for a Poisson reconstruction.
 *  - Mesh.ply which is the surface found using the Hoppe / Visual Hull blend.
 *
 *
 * @param argv
 * @param args
 * @return 0 on success
 */
int main(int argv, char** args)
{
    Stream::Message msg( "Visual Hull" );

    if(argv < 5 || argv > 6)
    {
        print_usage(msg);
        return 1;
    }

    std::string camtxt(args[1]);
    std::string imtxt(args[2]);
    std::string pointcloud(args[3]);
    std::string outfolder(args[4]);

    bool readBundler = true;
    bool doPcloudOpt = false;
    if(argv==6)
    {
        readBundler = atoi(args[5]) == 1;
    }
    if(argv==7)
    {
        doPcloudOpt = atoi(args[6]) == 1;
    }
    msg.print(1) << (readBundler ? "Using bundler format" : "Using octave format") << std::endl;

    auto camlist = readNameFile( camtxt );
    auto imlist = readNameFile( imtxt );

    if(camlist.size() == 0 || imlist.size() == 0)
    {
        msg.printerr( "The list contained no files" );
        return 1;
    }

    if( camlist.size() != imlist.size() )
    {
        msg.printerr( "The lists are not equal in size" );
        return 1;
    }

    std::vector< Cu::Matrix > cameras(camlist.size());
    std::vector< Image > images(camlist.size());

    auto camit = camlist.begin();
    auto imit = imlist.begin();
    uint32_t index=0;
    for( ;camit != camlist.end(); ++camit, ++imit)
    {
        Cu::Matrix P;

        if(readBundler)
        {
            P = Cu::read_matrix_bundler( *camit );
        }
        else
        {
            std::string name;
            auto Pp = Cu::read_matrix_octave( *camit, name );

            if(Pp.second.size() > 2)
            {
                msg.printerr(" P is not a matrix ");
                return 1;
            }
            P = Cu::ConvertToEigen(Pp.first.data(), Pp.second[0], Pp.second[1] );
        }

        cameras[index] = P;
        images[index].Load(*imit);

        ++index;

        msg.loadbar( static_cast<double>(index)/camlist.size() );
    }
    msg.print() << "Loading cloud" << std::endl;
    Mesh pcloud(pointcloud);

    msg.print(1) << "Removing redundant points from cloud" << std::endl;
    Carve::RemovePointsOutsideHull( pcloud, cameras, images, msg );

    msg.print() << "Computing start box" << std::endl;
    const double fractionPastPointCloud = 0.1;
    Carve::Box box = Carve::GetBoxFromCloud( cameras, images, msg, fractionPastPointCloud, pcloud.GetVertices());
    // Carve::Box box = Carve::ComputeInitialBox(cameras, images, msg);

    msg.print() << "Finding visual hull mesh" << std::endl;
    const uint32_t density = 300; // 300 is quite dense 100 is quick but does not do thin horse legs.
    Mesh mesh = Carve::FindMesh(cameras, images, box, msg, density );

    msg.print() << "Computing point cloud normals" << std::endl;
    mesh.ComputeNormalsFromPointsAndFaces();

    msg.print(1) << "Computing nearest neighbours" <<std::endl;
    NearestNeighbours nnmesh( NearestNeighbours::RandomTrees, mesh.GetVertices().data(),
                          3, mesh.GetVertices().size()/3);

    msg.print(1) << "Transferring normals" <<std::endl;
    pcloud.ApplyNormals( mesh, nnmesh );

    Mesh finalsurface;
    if(doPcloudOpt)
    {
        msg.print() << "Blending surface" << std::endl;
        const double K = 200; // Number of neighbours
        double radius = mesh.ComputeAverageEdgeLength();

        NearestNeighbours nnpcloud( NearestNeighbours::RandomTrees, pcloud.GetVertices().data(),
                              3, pcloud.GetVertices().size()/3);
        finalsurface = Carve::Blend( pcloud, &nnpcloud, mesh, &nnmesh, box, K, 50*radius, msg );
        //Mesh finalsurface = Carve::Hoppe( pcloud, &nnpcloud, box, K, 10*radius, msg );
    }

    msg.print(2) << "Write surfaces" <<std::endl;
    std::stringstream ssPcloud, ssHull, ssMesh;
    ssPcloud << outfolder << "/PCloud.ply";
    ssHull << outfolder << "/Hull.ply";

    if(doPcloudOpt)
    {
        ssMesh << outfolder << "/Mesh.ply";
        finalsurface.Write( ssMesh.str() );
    }

    mesh.Write(ssHull.str());
    pcloud.Write(ssPcloud.str());

    return 0;
}
