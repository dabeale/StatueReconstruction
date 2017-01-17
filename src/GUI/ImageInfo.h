#ifndef IMAGEINFO_H
#define IMAGEINFO_H

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

#include "PenPoint.h"
#include <QImage>
#include <QGraphicsScene>
#include <QGraphicsPixmapItem>
#include <QGraphicsBlurEffect>
#include <QPainter>

#include <vector>
#include <stdint.h>
#include <fstream>
#include <sstream>
#include <vector>

/**
 * \brief SegmentationGUI
 * A collection of classes and functions written in Qt for the graphical usr interface.
 */
namespace SegmentationGUI
{

/**
 * @brief The ImageInfo class
 * Contains all the information relevant for an image. This includes, depth map information, camera file and calibration matrices,
 * segmentations. It also stores all of the interactions that a user has made with the image in the interface, such as
 * selected foreground, background and edge points.
 */
class ImageInfo
{
public:
    /**
     * @brief ImageInfo
     * An explicit constructor for an image and its information.
     * @param fileName The path to the image.
     */
    ImageInfo( QString fileName );

    QImage m_image;     ///< A scaled version of the file on disk
    QImage m_segmentation;///< The produced segmentation
    QImage m_depthMap; ///< The depth map
    std::vector<double> m_depthMapValues; ///< Column major depth map values ( this is the size of the original file )
    std::vector<int32_t> m_indexValues; ///< The vertex indecies corresponding to the depth map
    std::vector<bool> m_pointCloudDepthVisibilities; ///< The visibilities of the point cloud in this view  ( this is the size of the original file )
    bool m_depthMapLoaded; ///< True if the depth map is loaded

    void Clear(); ///< Clear all of the point and colour buffers
    void FillColoursAroundPoint(const PenPoint& pp, std::vector< double >& pointlist , std::vector<double> &depthlist); ///< Fill in the colour buffer around the specified point

    std::vector< PenPoint > m_ForegroundPoints; ///< User sketched foreground points
    std::vector< PenPoint > m_BackgroundPoints; ///< User sketched background points
    std::vector< PenPoint > m_EdgePoints;       ///< User sketched edge points
    std::vector< double > m_ForegroundDepths; ///< User sketched foreground depths
    std::vector< double > m_BackgroundDepths; ///< User sketched background depths
    std::vector< double > m_EdgeDepths;       ///< User sketched edge depths

    std::vector< double > m_ForegroundColours; ///< A vector of user selected foreground colours (in rgb format)
    std::vector< double > m_BackgroundColours; ///< A vector of user selected background colours

    QString m_filename; ///< The name of the underlying image file

    /**
     * @brief DumpToOctave
     * Dump the user selected points to an octave readable file
     * @param namePrefix
     */
    void DumpToOctave(const std::string& namePrefix) const;

    /**
     * @brief WriteSegmentation
     * Write the segmentation to file
     * @param file The file to write to
     */
    void WriteSegmentation( std::string file );

    uint32_t getFileHeight() const;
    uint32_t getFileWidth() const;
private:
    uint32_t m_fileWidth; ///< The width of the original file on disk (the loaded image is scaled down).
    uint32_t m_fileHeight; ///< The height of the original file on disk (the loaded image is scaled down.)

    /**
     * @brief DumpList
     * Dump the entire pointlist to a file
     * @param pointList
     * @param name
     */
    void DumpList( const std::vector< PenPoint >& pointList, const std::string& name ) const;

};
}
#endif // IMAGEINFO_H
