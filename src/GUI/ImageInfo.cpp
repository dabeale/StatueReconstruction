#include "ImageInfo.h"

std::string GetFileFromString (const std::string& str)
{
  unsigned found = str.find_last_of("/\\");
  return str.substr(found+1);
}

void print_octave_matrix_header( std::ofstream& ostream, const std::string& name, const uint32_t rows, const uint32_t cols )
{
    ostream << " # Created by QTGui" << std::endl;
    ostream << " # name: " << GetFileFromString(name) << std::endl;
    ostream << " # type: matrix " << std::endl;
    ostream << " # rows : " << rows << std::endl;
    ostream << " # columns : " << cols << std::endl;
}

ImageInfo::ImageInfo(QString fileName ) :
        m_depthMapLoaded(false),
        m_ForegroundPoints(0),
        m_BackgroundPoints(0),
        m_EdgePoints(0),
        m_filename( fileName )
{
    QImage file(fileName);
    if(!file.isNull())
    {
        m_image = QImage( file.scaledToHeight(600));
        m_segmentation = QImage( m_image.size(), m_image.format() );
        m_depthMap = QImage( file.size(), m_image.format() );
        //m_image.setDevicePixelRatio(devicePixelRatio); This is not required when using QPainter
    }
    m_fileWidth = file.width();
    m_fileHeight = file.height();
}

void ImageInfo::Clear()
{
    m_ForegroundColours.clear();
    m_BackgroundColours.clear();
    m_ForegroundPoints.clear();
    m_BackgroundPoints.clear();
}

void ImageInfo::FillColoursAroundPoint(const PenPoint& pp , std::vector<double> &pointlist, std::vector<double> &depthlist)
{
    int32_t radius = pp.m_pen.width()/2;
    int32_t x = pp.m_point.x();
    int32_t y =pp.m_point.y();

    std::vector<double> colours;
    colours.reserve(3*4*radius*radius);
    std::vector<double> depths;
    depths.reserve(4*radius*radius);
    QPoint pt;

    for(int32_t i = x-radius; i < x+radius; ++i)
    for(int32_t j = y-radius; j < y+radius; ++j)
    if( i < m_image.width() && j < m_image.height() && i > 0 && j > 0)
    if( (i - x)*(i - x) + (j - y)*(j - y) < radius*radius )
    {
        pt.setX( i );
        pt.setY( j );
        QColor rgb(m_image.pixel( pt ));
        colours.push_back( rgb.redF() );
        colours.push_back( rgb.greenF() );
        colours.push_back( rgb.blueF() );
        if(m_depthMapLoaded)
        {
            double g = m_depthMapValues[ getFileHeight() * ( (i*getFileWidth()) / m_image.width()) + j*getFileHeight() / m_image.height() ];
            depths.push_back( g );
        }
    }

    pointlist.insert( pointlist.end(), colours.begin(), colours.end() );
    if(m_depthMapLoaded)
    {
        depthlist.insert( depthlist.end(), depths.begin(), depths.end());
    }
}

uint32_t ImageInfo::getFileHeight() const
{
    return m_fileHeight;
}

uint32_t ImageInfo::getFileWidth() const
{
    return m_fileWidth;
}

void ImageInfo::DumpToOctave(const std::string& namePrefix) const
{
    // Decompose the filename
    QStringList qsl = m_filename.split('/');
    qsl = qsl.back().split('.');

    if( m_ForegroundPoints.size() > 0 )
    {
        std::stringstream ss;
        ss << namePrefix << "Foreground" << qsl.front().toStdString();
        DumpList( m_ForegroundPoints, ss.str() );
    }

    if( m_BackgroundPoints.size() > 0 )
    {
        std::stringstream ss;
        ss << namePrefix << "Background" << qsl.front().toStdString();
        DumpList( m_BackgroundPoints, ss.str() );
    }

    if( m_EdgePoints.size() > 0 )
    {
        std::stringstream ss;
        ss << namePrefix << "Edge" << qsl.front().toStdString();
        DumpList( m_EdgePoints, ss.str() );
    }
}

inline QImage applyEffectToImage(QImage src, QGraphicsEffect *effect, int extent = 0)
{
    if(src.isNull()) return QImage();   //No need to do anything else!
    if(!effect) return src;             //No need to do anything else!
    QGraphicsScene scene;
    QGraphicsPixmapItem item;
    item.setPixmap(QPixmap::fromImage(src));
    item.setGraphicsEffect(effect);
    scene.addItem(&item);
    QImage res(src.size()+QSize(extent*2, extent*2), QImage::Format_ARGB32);
    res.fill(Qt::transparent);
    QPainter ptr(&res);
    scene.render(&ptr, QRectF(), QRectF( -extent, -extent, src.width()+extent*2, src.height()+extent*2 ) );
    return res;
}

void ImageInfo::WriteSegmentation( std::string file )
{
    QImage saveSeg = m_segmentation.scaled( m_fileWidth, m_fileHeight );
    QGraphicsBlurEffect *blur = new QGraphicsBlurEffect;
    blur->setBlurRadius(20);
    QImage result = applyEffectToImage(saveSeg, blur);
    result.save( QString(file.c_str()) );
}

void ImageInfo::DumpList( const std::vector< PenPoint >& pointList, const std::string& name ) const
{
    std::ofstream ostream;
    ostream.open(name);
    print_octave_matrix_header( ostream, name, pointList.size(), 3 );
    for( const PenPoint& pt : pointList )
    {
        ostream << (pt.m_point.x()*getFileWidth()) / m_image.width() <<
                   " " <<
                   (pt.m_point.y()*getFileHeight()) / m_image.height() <<
                   " " <<
                   pt.m_pen.widthF() << std::endl;
    }
    ostream.close();
}
