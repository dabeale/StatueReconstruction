#ifndef IMAGE_H
#define IMAGE_H

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

#include <string>
#include <vector>

#define _GLIBCXX_USE_CXX11_ABI 0
#include <opencv2/opencv.hpp>
//#include <opencv2/imgproc/imgproc.hpp>
//#include <opencv2/highgui/highgui.hpp>

// It would be possible to use CImg here instead. It is more lightweight.
// #include "CImg.h"

/**
 * @brief The Image class.
 * This class represents an image. It loads and saves using opencv, and contains and
 * opencv image and also a vector containing rgb values. The purpose is to put the data
 * in to contiguous format for certain algorithms such as fft and clustering to work properly.
 *
 * There are a few methods developed such as colour conversion routines, and resizing etc. They
 * are fairly primitive.
 *
 * The storage format is column major rgb values. Alpha is not stored.
 */
class Image
{
public:
    /**
     * @brief Construct an empty image.
     */
    Image();

    /**
     * @brief Construct an image from a filename. Opencv is used to load the image.
     * @param filename the path to the file.
     */
    Image( const std::string& filename );

    /**
     * @brief Construct an image from a vector of rgb values, between 0 and 1. This has been left open
     * to create images with a depth different than 3, although currently only 3 colours are supported,
     * with no alpha.
     * @todo Provide support for an alpha channel and also images of depth different than 3.
     * @param image The vector from which to create the image.
     * @param rows The number of rows
     * @param cols The number of columns
     * @param depth The colour depth
     */
    Image( const std::vector<double>& image, const uint32_t rows, const uint32_t cols, const uint32_t depth);

    /**
     * @brief Create an empty image with the specified number of rows and columns and depth.
     * @param rows The number of rows
     * @param cols The number of columns
     * @param depth The depth (=3)
     */
    Image( const uint32_t rows, const uint32_t cols, const uint32_t depth);
    ~Image(); ///< Destroy the Image object.

    /**
     * @brief Load an image from disk.
     * @param filename The path of the image
     */
    void Load( const std::string& filename );

    inline double * Ptr() {return m_image.data();}              ///< Get pointer to the data
    inline const double * Ptr() const {return m_image.data();}  ///< Get const pointer to the data
    inline std::vector<double>& Data() {return m_image;}        ///< Get vector of data
    inline cv::Mat& GetCVImage() {return m_imageCV;}            ///< Return the cv image that the m_image was generated from.

    inline uint32_t rows() const{return m_Rows;}                ///< Get the number of rows
    inline uint32_t cols() const {return m_Cols;}               ///< Get the number of columns
    inline uint32_t depth() const{return m_Depth;}              ///< Get the depth (should always be 3)

    Image ToGrey() const;                                       ///< Convert the image to greyscale
    Image ToLUV() const;                                        ///< Convert the image to LUV
    Image ToHSV() const;                                        ///< Convert the image to HSV

    /**
     * @brief Write the image to disk.
     * The image is written using the opencv
     * interface, and so it supports all of the same filetypes.
     * @param filename the path to the image.
     */
    void WriteToDisk( const std::string& filename);

    /**
     * @brief Resize the image
     * @param rows The new number of rows
     * @param cols The new number of columns
     */
    void Resize(uint32_t rows, uint32_t cols);

    void ScaleToHeight(uint32_t rows); ///< \brief Scale the image to a new height, while proportionally scaling the width.
    void ScaleToWidth(uint32_t cols); ///< \brief Scale the image to a new width, while proportionally scaling the height.

    /**
     * \brief GetInterpolatedPatch.
     * Get a a vector of values in a neighborhood of a point on the image. The returned patch is
     * of size 3*patchsize^2, and is stored in rgb - column major format. The inputs to this function
     * need not be integral since the colours are interpolated, to allow for optimisation in the image.
     * @param i the row number
     * @param j the column number
     * @return A vector of values
     */
    template<typename T>
    inline std::vector<T> GetInterpolatedPatch( const T i, const T j, const uint32_t patchSize) const
    {
        std::vector<T> patch( patchSize*patchSize*3);
        T r,g,b;
        for( uint32_t k=0; k<patchSize; ++k )
        for( uint32_t l=0; l<patchSize; ++l )
        {
            GetInterpolatedColour( i-(k - patchSize/2) , j-(l - patchSize/2),  r, g, b );
            patch[(k*patchSize + l)*3] = r;
            patch[(k*patchSize + l)*3 + 1] = g;
            patch[(k*patchSize + l)*3 + 2] = b;
        }
        return patch;
    }

    /**
     * \brief GetInterpolatedColour.
     * Get the interpolated colour at a particular point on the image. Bilinear interpolation is
     * used to compute the value of the image at non-integral row and column.
     * @param i the row
     * @param j the column
     * @param r a reference to the red value
     * @param g a reference to the green value
     * @param b a reference to the blue value
     */
    template<typename T>
    inline void GetInterpolatedColour( T i, T j, T& r, T& g, T& b ) const
    {
        if( i < 0 ) i=0;
        if( i > m_Rows-1 ) i = m_Rows-1;
        if( j < 0 ) j=0;
        if( j > m_Cols-1 ) j = m_Cols-1;

        // This code uses bilinear interpolation, assuming that the colours are on the corners of a unit square
        uint32_t iceil = static_cast<uint32_t>( std::ceil(i) );
        uint32_t ifloor = static_cast<uint32_t>(std::floor(i));
        uint32_t jceil = static_cast<uint32_t>( std::ceil(j) );
        uint32_t jfloor = static_cast<uint32_t>(std::floor(j));
        T ifr = i - ifloor;
        T jfr = j - jfloor;

        r = (1-jfr)*(1-ifr)*m_image[ (jfloor*m_Rows + ifloor)*m_Depth ] +
               jfr *(1-ifr)*m_image[ (jceil*m_Rows + ifloor)*m_Depth ] +
            (1-jfr)*   ifr *m_image[ (jfloor*m_Rows + iceil)*m_Depth ] +
               jfr*    ifr *m_image[ (jceil*m_Rows + iceil)*m_Depth ];
        g = (1-jfr)*(1-ifr)*m_image[ (jfloor*m_Rows + ifloor)*m_Depth + 1 ] +
               jfr *(1-ifr)*m_image[ (jceil*m_Rows + ifloor)*m_Depth + 1 ] +
            (1-jfr)*   ifr *m_image[ (jfloor*m_Rows + iceil)*m_Depth + 1] +
               jfr*    ifr *m_image[ (jceil*m_Rows + iceil)*m_Depth + 1];
        b = (1-jfr)*(1-ifr)*m_image[ (jfloor*m_Rows + ifloor)*m_Depth + 2 ] +
               jfr *(1-ifr)*m_image[ (jceil*m_Rows + ifloor)*m_Depth + 2 ] +
            (1-jfr)*   ifr *m_image[ (jfloor*m_Rows + iceil)*m_Depth + 2] +
               jfr*    ifr *m_image[ (jceil*m_Rows + iceil)*m_Depth + 2];
    }

    /**
     * @brief GetColour.
     * Get the colour at a particular row and column, with double precision, between 0 and 1.
     * @param i the row
     * @param j the column
     * @param r a reference to the red value
     * @param g a reference to the green value
     * @param b a reference to the blue value
     */
    inline void GetColour(const uint32_t i, const uint32_t j, double& r, double& g, double& b) const
    {
        r = m_image[ (j*m_Rows + i)*m_Depth ];
        g = m_image[ (j*m_Rows + i)*m_Depth + 1];
        b = m_image[ (j*m_Rows + i)*m_Depth + 2];
    }

    /**
     * @brief GetColour.
     * Get the colour at a particular row and column, with integer precision, between 0 and 255.
     * @param i the row
     * @param j the column
     * @param r a reference to the red value
     * @param g a reference to the green value
     * @param b a reference to the blue value
     */
    inline void GetColour(const uint32_t i, const uint32_t j, uint8_t& r, uint8_t& g, uint8_t& b) const
    {
        r = static_cast<uint8_t>(255*m_image[ (j*m_Rows + i)*m_Depth ]);
        g = static_cast<uint8_t>(255*m_image[ (j*m_Rows + i)*m_Depth + 1]);
        b = static_cast<uint8_t>(255*m_image[ (j*m_Rows + i)*m_Depth + 2]);
    }

    /**
     * @brief SetColour.
     * Set the colour at a particular row and column, with double precision, between 0 and 1.
     * @param i the row
     * @param j the column
     * @param r the red value
     * @param g the green value
     * @param b the blue value
     */
    inline void SetColour(const uint32_t i, const uint32_t j, const double r, const double g, const double b)
    {
        m_image[ (j*m_Rows + i)*m_Depth ] = r;
        m_image[ (j*m_Rows + i)*m_Depth + 1] = g;
        m_image[ (j*m_Rows + i)*m_Depth + 2] = b;
    }

    /**
     * @brief SetColour.
     * Set the colour at a particular row and column, with integer precision, between 0 and 1.
     * @param i the row
     * @param j the column
     * @param r the red value
     * @param g the green value
     * @param b the blue value
     */
    inline void SetColour(const uint32_t i, const uint32_t j, const uint8_t r, const uint8_t g, const uint8_t b)
    {
        m_image[ (j*m_Rows + i)*m_Depth ] = static_cast<double>(r)/255.0;
        m_image[ (j*m_Rows + i)*m_Depth + 1] = static_cast<double>(g)/255.0;
        m_image[ (j*m_Rows + i)*m_Depth + 2] = static_cast<double>(b)/255.0;
    }

    std::vector<double> GetPermutation() const; ///< \brief Convert the image from 3xMxN to MxNx3 for filtering

private:
    std::vector<double> m_image; ///< The image in contiguous space
    cv::Mat m_imageCV;           ///< The image in opencv format. This is the first thing that is loaded, the colour conversions do not apply.

    uint32_t m_Rows;             ///< The number of rows
    uint32_t m_Cols;             ///< The number of columns
    uint32_t m_Depth;            ///< The number of colours

    void FillCVImage();          ///< \brief Fill m_imageCV with m_image
    void FillArrayWithCVImage(); ///< \brief Fill m_image with m_CV image (performed on construction)
};

#endif // IMAGE_H
