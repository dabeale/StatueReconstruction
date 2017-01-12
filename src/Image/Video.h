#ifndef VIDEO_H
#define VIDEO_H

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

#include <vector>
#include <string>

//#define _GLIBCXX_USE_CXX11_ABI 0
#include <opencv2/opencv.hpp>
//#include <opencv2/imgproc/imgproc.hpp>
//#include <opencv2/highgui/highgui.hpp>


/**
 * @brief The Video class
 * This class uses opencv to load a video into a contiguous
 * chunk of memory. Care must be taken on the size of the input video.
 */
class Video
{
public:
    Video( const std::string filename);

    inline double* Ptr(){return m_vid.data();} ///< Return a pointer to the first point in memory
    inline std::vector<double>& data(){ return m_vid;} ///< Return a data reference

    inline uint32_t Rows(){return m_Rows;}
    inline uint32_t Cols(){return m_Cols;}
    inline uint32_t Depth(){return m_Depth;}
    inline uint32_t Frames(){return m_frames;}

private:
    std::vector<double> m_vid;
    std::vector< cv::Mat > m_cvFrames;

    uint32_t m_Rows;
    uint32_t m_Cols;
    uint32_t m_Depth;
    uint32_t m_frames;
};

#endif // VIDEO_H

