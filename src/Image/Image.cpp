#include "Image.h"

Image::Image() : m_imageCV()
{

}

Image::Image(const std::string &filename) :
    m_imageCV()
{
    Load(filename);
}

Image::Image( const std::vector<double>& image, const uint32_t rows, const uint32_t cols, const uint32_t depth) :
    m_image(image), m_Rows(rows), m_Cols(cols), m_Depth(depth)
{
    FillCVImage();
}

Image::Image( const uint32_t rows, const uint32_t cols, const uint32_t depth) :
    m_image(rows*cols*depth, 0.0), m_Rows(rows), m_Cols(cols), m_Depth(depth)
{
    FillCVImage();
}

void Image::FillArrayWithCVImage()
{
    uchar r,g,b,m;
    cv::Vec3b intensity;

    m = std::numeric_limits<uchar>::max();

    for( uint32_t index=0,j=0,i; j<m_Cols; ++j  )
    for( i=0; i<m_Rows; ++i )
    {
        intensity = m_imageCV.at<cv::Vec3b>(i, j);

        b = intensity.val[0];
        g = intensity.val[1];
        r = intensity.val[2];

        m_image[index] = static_cast<double>(r) / static_cast<double>(m);
        m_image[index + 1] = static_cast<double>(g) / static_cast<double>(m);
        m_image[index + 2] = static_cast<double>(b) / static_cast<double>(m);
        for(uint32_t k=0; k<3; ++k)
        {
            if(m_image[index] > 1)m_image[index]=1;
            if(m_image[index] < 0)m_image[index]=0;
        }
        index+=3;
    }
}

void Image::Load( const std::string& filename )
{
    m_imageCV = cv::imread(filename, cv::IMREAD_COLOR);
    m_Rows = m_imageCV.rows;
    m_Cols = m_imageCV.cols;
    m_Depth = m_imageCV.channels();
    m_image = std::vector<double>(m_Rows*m_Cols*m_Depth,0.0);

    FillArrayWithCVImage();
}

void Image::FillCVImage()
{
    if(m_Depth == 3)
    {
        m_imageCV = cv::Mat(m_Rows, m_Cols, cv::DataType<cv::Vec3b>::type);
        m_imageCV.setTo(cv::Scalar::all(0.0));
        uchar m = std::numeric_limits<uchar>::max();
        for( uint32_t index=0,j=0,i; j<m_Cols; ++j  )
        for( i=0; i<m_Rows; ++i )
        {
            cv::Vec3b intensity;
            intensity.val[0] = static_cast<uchar>(m*m_image[index+2]);
            intensity.val[1] = static_cast<uchar>(m*m_image[index+1]);
            intensity.val[2] = static_cast<uchar>(m*m_image[index]);
            m_imageCV.at<cv::Vec3b>(i, j) = intensity;

            index += 3;
        }
    }
    else
    {
        m_imageCV = cv::Mat(m_Rows, m_Cols, cv::DataType<uchar>::type);
        m_imageCV.setTo(cv::Scalar::all(0.0));
        uchar m = std::numeric_limits<uchar>::max();
        for( uint32_t index=0,j=0,i; j<m_Cols; ++j  )
        for( i=0; i<m_Rows; ++i )
        {
            m_imageCV.at<uchar>(i, j) = static_cast<uchar>(m*m_image[index]);
            index ++;
        }
    }
}

void Image::WriteToDisk( const std::string& filename)
{
    FillCVImage();
    std::vector<int> compression_params;
    compression_params.push_back(CV_IMWRITE_PNG_COMPRESSION);
    compression_params.push_back(0);

    cv::imwrite(filename, m_imageCV, compression_params);
}

void Image::Resize(uint32_t rows, uint32_t cols)
{
    std::vector<double> newim(rows*cols*m_Depth);

//#pragma omp parallel for
    uint32_t index=0;
    for(uint32_t j=0; j<cols; ++j)
    {
        uint32_t tindex,i;
        for(i=0; i<rows; ++i)
        {
            tindex = m_Depth*(((j*m_Cols)/cols)*m_Rows + (i*m_Rows)/rows);

            for(uint32_t c =0; c<m_Depth; ++c)
                newim[m_Depth*index + c] = m_image[ tindex + c];

            ++index;
        }
    }
    m_image = newim;
    m_Rows = rows;
    m_Cols = cols;
}

void Image::ScaleToHeight(uint32_t rows)
{
    Resize(rows, (rows*m_Cols) / m_Rows);
}

void Image::ScaleToWidth(uint32_t cols)
{
    Resize( (cols*m_Rows) / m_Cols, cols);
}

Image Image::ToGrey() const
{
    std::vector<double> imgrey(m_Rows*m_Cols, 0.0);

#pragma omp parallel for
    for(uint32_t k=0; k<m_Rows*m_Cols; ++k)
    {
        imgrey[k] = (0.2*m_image[k*3] + 0.3*m_image[k*3 + 1] + 0.5*m_image[k*3 + 2]);
    }

    return Image(imgrey, m_Rows, m_Cols, 1);
}

Image Image::ToLUV() const
{
    std::vector<double> imluv(m_Rows*m_Cols*3, 0.0);

#pragma omp parallel for
    for(uint32_t k=0; k<m_Rows*m_Cols; ++k)
    {
        double r,g,b,x,y,z,L,u,v;
        r = m_image[k*3];
        g = m_image[k*3 + 1];
        b = m_image[k*3 + 2];

        x = 0.412453*r + 0.357580*g + 0.180423*b;
        y = 0.212671*r + 0.715160*g + 0.072169*b;
        z = 0.019334*r + 0.119193*g + 0.950227*b;


        L = ( y> 0.008856 ) ? 116* std::pow(y, 0.333333333333) : 903.3*y;
        if(x==0 && y ==0 && z==0)
        {
            u=0;
            v=0;
        }
        else
        {
            u = 4*x/(x + 15*y + 3*z);
            v = 9*y/(x + 15*y + 3*z);
        }
        u = 13*L*(u - 0.19793943);
        v = 13*L*(v - 0.46831096);

        imluv[3*k] = L;
        imluv[3*k+1] = u;
        imluv[3*k+2] = v;
    }

    return Image(imluv, m_Rows, m_Cols, 3);
}

Image Image::ToHSV() const
{
    std::vector<double> imhsv(m_Rows*m_Cols*3, 0.0);

#pragma omp parallel for
    for(uint32_t k=0; k<m_Rows*m_Cols; ++k)
    {
        double r,g,b,h,s;
        r = m_image[k*3];
        g = m_image[k*3 + 1];
        b = m_image[k*3 + 2];

        const double& v = std::max(r, std::max(g, b));
        const double& vm = std::min(r, std::min(g, b));

        if(&v == &r)
        {
            h = 60*(g-b)/(v-vm);
        }
        else if(&v == &g)
        {
            h = 120 + 60*(b-r)/(v-vm);
        }
        else
        {
            h = 240 + 60*(r-g)/(v-vm);
        }

        s = (v < 1e-13) ? 1.0 - (vm/v) : 0.0;

        imhsv[k*3] = h;
        imhsv[k*3+1] = s;
        imhsv[k*3+2] = v;
    }

    return Image(imhsv, m_Rows, m_Cols, 1);
}

std::vector<double> Image::GetPermutation() const
{
    std::vector<double> perm(m_Rows*m_Cols*3);

#pragma omp parallel for
    for(uint32_t k=0; k<m_Rows*m_Cols; ++k)
    {
        perm[k]                   = m_image[k*3];
        perm[k + m_Rows*m_Cols]   = m_image[k*3 + 1];
        perm[k + 2*m_Rows*m_Cols] = m_image[k*3 + 2];
    }
    return perm;
}

Image::~Image()
{

}
