
#include "Video.h"

Video::Video( const std::string filename) :
    m_vid(), m_cvFrames()
{
    cv::VideoCapture cap(filename); // open the default camera
    if(!cap.isOpened())  // check if we succeeded
    {
        std::cerr << "Video: Could not open video" <<std::endl;
        assert(0);
    }

    m_frames = cap.get(CV_CAP_PROP_FRAME_COUNT);
    m_Rows = cap.get(CV_CAP_PROP_FRAME_HEIGHT);
    m_Cols = cap.get(CV_CAP_PROP_FRAME_WIDTH);
    m_Depth = 3;

    const uint32_t step = m_Rows*m_Cols*m_Depth;

    m_cvFrames.resize(m_frames);
    m_vid.resize(m_frames*step);

    for(uint32_t k,n=0; n < m_frames; ++n)
    {
        cap >> m_cvFrames[n]; // get a new frame from camera
        if( m_cvFrames[n].empty())
        {
            break;
        }

        for(k=0; k<step; ++k)
        {
            m_vid[n*step + k] = static_cast<double>(m_cvFrames[n].at<uchar>(k)) / 255.0;
        }
    }
}
