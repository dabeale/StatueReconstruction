#include "ImageViewer.hpp"

ImageViewer::ImageViewer(QString &startDirectory, QWidget *parent) :
    QLabel(parent),
    m_x(0),
    m_y(0),
    m_mousex(0),
    m_mousey(0),
    m_size(5),
    m_latch(false),
    m_paintInitialised(false),
    m_viewMode(Annotation),
    m_color(Qt::red),
    m_sketchtype(Edge),
    m_imageInfos(0, ImageInfo(".")),
    m_current(m_imageInfos.begin()),
    m_startDirectory(startDirectory),
    m_msg("ImageViewer"),
    m_PointsDrawn(false),
    m_ModelComputed(false),
    m_depthMapComputed(false),
    m_sgm(m_imageInfos, m_3dmesh, m_ProjectionMatrices),
    m_mrfparamswindow(new MRFParametersWindow(this))
{
    m_mrfparamswindow->SetAlpha( m_sgm.GetParameters().m_Alpha[0]);
    m_mrfparamswindow->SetSampleN(m_sgm.GetParameters().m_NSample);
    m_mrfparamswindow->SetPairWiseWeight(m_sgm.GetParameters().m_edgeWeight);
    m_mrfparamswindow->SetUnaryWeight(m_sgm.GetParameters().m_terminalWeight);
    connect(m_mrfparamswindow, SIGNAL(CloseButtonPushed()), this, SLOT(UpdateMRFParameters()));
}

void ImageViewer::ShowParametersWindow()
{
    m_mrfparamswindow->show();
}

void ImageViewer::UpdateMRFParameters()
{
    m_sgm.GetParameters().m_Alpha = {m_mrfparamswindow->GetAlpha(),m_mrfparamswindow->GetAlpha() };
    m_sgm.GetParameters().m_type = (m_mrfparamswindow->GetType() == MRFParametersWindow::GraphCut) ? SegmentationModelParameters::GraphCut : SegmentationModelParameters::ICM;
    m_sgm.GetParameters().m_NSample = m_mrfparamswindow->GetSampleN();
    m_sgm.GetParameters().m_edgeWeight = m_mrfparamswindow->GetPairWiseEight();
    m_sgm.GetParameters().m_terminalWeight = m_mrfparamswindow->GetUnaryWeight();
}

void ImageViewer::SetSketchType( SketchType type )
{
    m_sketchtype = type;
}

void ImageViewer::SetColor(Qt::GlobalColor color)
{
    m_color = color;
}

void ImageViewer::SetSize(uint32_t size)
{
    m_size = size;
}

void ImageViewer::CycleLeft()
{
    if(m_imageInfos.size() > 0)
    {
        CycleIteratorLeft<ImageList, ImageList::iterator>(m_current, m_imageInfos);
        RefreshImage();
    }
}

void ImageViewer::CycleRight()
{
    if(m_imageInfos.size() > 0)
    {
        CycleIteratorRight<ImageList, ImageList::iterator>(m_current, m_imageInfos);
        RefreshImage();
    }
}

void ImageViewer::SetPosition(int frame)
{
    if(m_imageInfos.size() > 0)
    {
        uint32_t position = m_current - m_imageInfos.begin();
        if(frame < position )
        {
            while( frame < position )
            {
                --m_current;
                position = m_current - m_imageInfos.begin();
                if(position == 0 || position == m_imageInfos.size()-1)
                    break;
            }
        }
        else
        {
            while( frame > position )
            {
                ++m_current;
                position = m_current - m_imageInfos.begin();
                if(position == 0 || position == m_imageInfos.size()-1)
                    break;
            }
        }
    }
    RefreshImage();
}

template<class List, class Iterator>
void ImageViewer::CycleIteratorLeft(Iterator& iterator, List &imageList)
{
    if( --iterator < imageList.begin() )
    {
        iterator = imageList.end();
        iterator--;
    }
}

template<class List, class Iterator>
void ImageViewer::CycleIteratorRight(Iterator& iterator,  List &imageList)
{
    if( ++iterator == imageList.end() )
    {
        iterator = imageList.begin();
    }
}

void ImageViewer::RepaintPointBuffer()
{
    repaint();
}

void ImageViewer::mousePressEvent(QMouseEvent *e)
{
    QLabel::mousePressEvent(e);
}


void ImageViewer::mouseReleaseEvent(QMouseEvent *e)
{
    QLabel::mouseReleaseEvent(e);
}

void ImageViewer::mouseMoveEvent(QMouseEvent *e)
{
    QLabel::mouseMoveEvent(e);

    m_mousex = e->x();
    m_mousey = e->y();
    if( e->buttons() & Qt::LeftButton &&
        (m_x - m_mousex)*(m_x - m_mousex) + (m_y - m_mousey)*(m_y - m_mousey) > (m_size*m_size) / 4)
    {
        m_x = m_mousex;
        m_y = m_mousey;
        m_latch = true;
        update();
        repaint();
    }

    /*if(m_imageInfos.size() > 0)
    {
        repaint();
    }*/
}

void ImageViewer::SetPaintInitialised(const bool initialised)
{
    m_paintInitialised = initialised;
}

void ImageViewer::DrawAnnotation( QPainter& painter )
{
    float widthRatio = (float(m_current->m_image.width()) / float(geometry().width()));
    float heightRatio = (float(m_current->m_image.height()) / float(geometry().height()));

    painter.setOpacity(1);
    if(m_paintInitialised)
    {
        for (const auto p1 : m_current->m_ForegroundPoints)
        {
            QPoint p((int)std::floor(p1.m_point.x()/widthRatio), (int)std::floor(p1.m_point.y()/heightRatio));
            painter.setPen(p1.m_pen);
            painter.drawPoint(p);
        }

        for (const auto p1 : m_current->m_BackgroundPoints)
        {
            QPoint p((int)std::floor(p1.m_point.x()/widthRatio), (int)std::floor(p1.m_point.y()/heightRatio));
            painter.setPen(p1.m_pen);
            painter.drawPoint(p);
        }

        for (const auto p1 : m_current->m_EdgePoints)
        {
            QPoint p((int)std::floor(p1.m_point.x()/widthRatio), (int)std::floor(p1.m_point.y()/heightRatio));
            painter.setPen(p1.m_pen);
            painter.drawPoint(p);
        }
    }
}

void ImageViewer::paintEvent(QPaintEvent *e)
{
    if(m_paintInitialised)
    {
        float widthRatio = (float(m_current->m_image.width()) / float(geometry().width()));
        float heightRatio = (float(m_current->m_image.height()) / float(geometry().height()));

        if(m_latch)
        {
            m_PointsDrawn = true;
            QColor qc(m_color);
            qc.setAlpha(100);

            QPen paintpen(qc);
            paintpen.setWidth(m_size);
            paintpen.setCapStyle(Qt::RoundCap);

            QPoint p1((int)std::floor(widthRatio*m_x), (int)std::floor(heightRatio*m_y));
            PenPoint pp;
            pp.m_pen = paintpen;
            pp.m_point = p1;

            switch( m_sketchtype )
            {
                case Foreground :
                    m_current->m_ForegroundPoints.push_back(pp);
                    m_current->FillColoursAroundPoint( pp, m_current->m_ForegroundColours, m_current->m_ForegroundDepths );
                break;
                case Background :
                    m_current->m_BackgroundPoints.push_back(pp);
                    m_current->FillColoursAroundPoint( pp, m_current->m_BackgroundColours, m_current->m_BackgroundDepths );
                break;
                case Edge :
                    m_current->m_EdgePoints.push_back(pp);
                break;
            }

            m_latch=false;
        }

        QPainter painter(this);

        switch(m_viewMode)
        {
            case DepthMap :
            {
                painter.drawImage(QRect(0,0,this->width(), this->height()), m_current->m_depthMap);
            } break;
            case Image :
            {
                painter.drawImage(QRect(0,0,this->width(), this->height()), m_current->m_image);
            } break;
            case Both :
            {
                painter.setOpacity(0.5);
                painter.drawImage(QRect(0,0,this->width(), this->height()), m_current->m_segmentation);
                painter.drawImage(QRect(0,0,this->width(), this->height()), m_current->m_image);
                DrawAnnotation(painter);
            } break;
            case Segmentation :
            {
                painter.setOpacity(0.5);
                painter.drawImage(QRect(0,0,this->width(), this->height()), m_current->m_segmentation);
                painter.drawImage(QRect(0,0,this->width(), this->height()), m_current->m_image);
            } break;
            case Annotation :
            {
                painter.drawImage(QRect(0,0,this->width(), this->height()), m_current->m_image);
                DrawAnnotation(painter);
            } break;
        }

        if(m_imageInfos.size() > 0)
        {
            //int dpr = devicePixelRatio();
            QPoint p((int)std::floor(m_x),
                     (int)std::floor(m_y));

            QColor qc(m_color);
            qc.setAlpha(100);

            QPen paintpen(qc);
            paintpen.setWidth(m_size);
            paintpen.setCapStyle(Qt::RoundCap);

            painter.setPen(paintpen);
            painter.drawPoint(p);
        }
    }

    QLabel::paintEvent(e);
}

void ImageViewer::SaveEdges()
{
    m_msg.print() << "Compute Edges" << std::endl;
    QString dir = QFileDialog::getExistingDirectory(this, tr("Open Directory"),
                     m_startDirectory,
                     QFileDialog::ShowDirsOnly
                     | QFileDialog::DontResolveSymlinks);
    m_startDirectory = dir;

    bool ok;
    double scale = QInputDialog::getDouble(this, tr("QInputDialog::getDouble()"),
                                        tr("Scale"), 30,
                                        -10000,
                                        10000,
                                        1, &ok);
    int perc =  0;
    if( ok && m_imageInfos.size() > 0 )
    {
        for( const ImageInfo imageinfo : m_imageInfos )
        {
            if( imageinfo.m_ForegroundPoints.size() > 0 ||
                imageinfo.m_BackgroundPoints.size() > 0 ||
                imageinfo.m_EdgePoints.size() > 0 )
            {
                // Decompose the filename
                QStringList qsl = imageinfo.m_filename.split('/');
                qsl = qsl.back().split('.');

                std::cout << (100*(perc++)) / m_imageInfos.size() << " pc complete \r";
                std::cout.flush();
            }
        }
        std::cout << "100 pc complete" << std::endl;
    }
    else
    {
        m_msg.print() << "Command exit" << std::endl;
    }
}

void ImageViewer::DumpToOctave()
{
    m_msg.print() << "Dump points to octave file " << std::endl;
    QString dir = QFileDialog::getExistingDirectory(this, tr("Open Directory"),
                     m_startDirectory,
                     QFileDialog::ShowDirsOnly
                     | QFileDialog::DontResolveSymlinks);
    m_startDirectory = dir;
    std::stringstream ss;
    ss << dir.toStdString() << "/";

    for( auto ii : m_imageInfos )
    {
        ii.DumpToOctave(ss.str());
    }
}

void  ImageViewer::SaveSegmentation()
{
    m_msg.print() << "Save segmentation " << std::endl;
    QString dir = QFileDialog::getExistingDirectory(this, tr("Open Directory"),
                     m_startDirectory,
                     QFileDialog::ShowDirsOnly
                     | QFileDialog::DontResolveSymlinks);
    m_startDirectory = dir;

    uint32_t index=0;
    for( auto ii : m_imageInfos )
    {
        std::stringstream ss;
        ss << dir.toStdString() << "/";
        ss << "Segmentation" << index++ << ".png";
        ii.WriteSegmentation( ss.str() );
    }
}

void  ImageViewer::LoadMesh()
{
    m_msg.print() << "Load mesh " << std::endl;
    QString plyfile = QFileDialog::getOpenFileName(this, tr("Open Mesh"),
                     m_startDirectory, tr("Stanford PLY (*.ply)"));

    m_3dmesh = Mesh(plyfile.toStdString());
}

void ImageViewer::LoadCameras()
{
    m_msg.print() << "Load cameras " << std::endl;
    QStringList cameras = QFileDialog::getOpenFileNames(this, tr("Open Cameras"),
                     m_startDirectory, tr("Bundler Cameras (*.txt)"));

    if(cameras.size() > m_imageInfos.size())
    {
        m_msg.printerr( "Warning:: ImageViewer::LoadCameras - There are more cameres than images");
    }
    if(cameras.size() < m_imageInfos.size())
    {
        m_msg.printerr( "Error:: ImageViewer::LoadCameras - There are fewer cameres than images");
        return;
    }

    m_ProjectionMatrices.resize(cameras.size());

    auto it = m_ProjectionMatrices.begin();
    auto itim = m_imageInfos.begin();
    for( auto& str : cameras )
    {
         /*Math::Matrix M( 3,3, {1.0, 0.0, static_cast<double>(itim->m_image.width())/static_cast<double>(itim->getFileWidth()),
                               0.0, 1.0, static_cast<double>(itim->m_image.height())/static_cast<double>(itim->getFileHeight()),
                               0.0, 0.0, 1.0});*/
        *(it++) = Utils::read_matrix_bundler( str.toStdString() );
         ++itim;
    }

}

bool ImageViewer::AddImages(const QStringList& fileNames)
{
    if( fileNames.length() > 0 )
    {
        int perc = 0;
        for ( const QString fileName : fileNames )
        {
            ImageInfo iinf(fileName);
            m_imageInfos.push_back(iinf);

            std::cout << (100*(perc++)) / fileNames.length() << " pc complete \r";
            std::cout.flush();
        }
        m_current = m_imageInfos.end();
        m_current--;
        std::cout << std::endl;
        return true;
     }
    else
    {
        return false;
    }
}

void ImageViewer::OpenFile()
{
    m_msg.print() << "Open files" << std::endl;
    QStringList fileNames = QFileDialog::getOpenFileNames(this,
        tr("Open Image"), m_startDirectory, tr("Image Files (*.png *.jpg *.bmp)"));
    fileNames.sort();

    if( AddImages(fileNames) )
    {
        m_startDirectory = QFileInfo(*fileNames.begin()).path();
        SetSketchType(Edge);
        RefreshImage();
        repaint();
        SetPaintInitialised(true);

        emit ResizeEvent();
        emit ImagesLoaded();
    }
}

void ImageViewer::GetDimensions( uint32_t& height, uint32_t& width ) const
{
    if(m_imageInfos.size() > 0)
    {
        height = m_imageInfos[0].getFileHeight();
        width = m_imageInfos[0].getFileWidth();
    }
    else
    {
        std::cerr << "ImageViewer::GetDimensions:: There are no imageinfos" << std::endl;
        height= 0;
        width = 0;
    }
}

/**
 * @brief MainWindow::RefreshImage
 * Refreshes the image - used for clearing any other drawn points
 */
void ImageViewer::RefreshImage()
{
    repaint();
}

void ImageViewer::Clear()
{
    if(m_imageInfos.size()>0)
    {
        m_current->Clear();
        RefreshImage();
    }
}

void ImageViewer::ClearAll()
{
    m_msg.print() << "Clear the whole annotation buffer" << std::endl;
    for ( auto& i : m_imageInfos )
    {
        i.Clear();
    }
    RefreshImage();
    m_PointsDrawn = false;
}

void ImageViewer::Select()
{
    SetColor( Qt::red );
    SetSketchType( Edge );
}

void ImageViewer::SelectForeground()
{
    SetColor( Qt::green );
    SetSketchType( Foreground );
}

void ImageViewer::SelectBackground()
{
    SetColor( Qt::blue );
    SetSketchType( Background );
}

uint32_t ImageViewer::GetNumberOfImages() const
{
    return m_imageInfos.size();
}

ImageViewer::ImageList& ImageViewer::GetImageInfos()
{
    return m_imageInfos;
}

void ImageViewer::InitialiseBuffer()
{
    SetSketchType( Edge );
    RefreshImage();
    repaint();
    SetPaintInitialised(true);
}

void ImageViewer::ComputeModel()
{
    m_msg.print() << "Computing Segmentation" <<std::endl;
    if(m_PointsDrawn)
    {
        m_sgm.CreateData( );
        if(m_depthMapComputed)
        {
            m_sgm.ComputeDepthModels();
            m_sgm.ComputePointProbabilities();
        }
        m_sgm.RunSegmentation(1000, 0.001);
        m_ModelComputed = true;
    }
    else
    {
        m_msg.printerr("Error :: ImageViewer::ComputeModel() - Points have not been drawn");
    }
    repaint();
    emit ModelComputed();
}

void ImageViewer::RunMRF()
{
    m_msg.print() << "Running MRF" <<std::endl;
    if(m_ModelComputed)
    {
        m_sgm.EvaluateMRF(m_msg);
    }
    else
    {
        m_msg.printerr("Error :: ImageViewer::RunMRF() - The model has not been computed");
    }
    repaint();
}

void ImageViewer::ComputeDepthMaps()
{
    m_msg.print() << "Computing depth maps" <<std::endl;
    if(m_3dmesh.GetNVerts() == 0)
    {
        m_msg.printerr("Error:: ImageViewer::ComputeDepthMaps() - there are no vertices in the mesh");
        return;
    }
    if(m_ProjectionMatrices.size() == 0)
    {
        m_msg.printerr("Error:: ImageViewer::ComputeDepthMaps() - there are no projection matrices loaded");
        return;
    }

    auto itp = m_ProjectionMatrices.begin();
    for( auto& image : m_imageInfos )
    {
        uint32_t M = image.getFileHeight();
        uint32_t N = image.getFileWidth();
        Buffer::ProjectionMatrix P(itp->GetArr());
        Buffer::ZBuffer buf({
                                &m_3dmesh.GetVertices()[0],
                                &m_3dmesh.GetFaces()[0],
                                m_3dmesh.GetNVerts(),
                                m_3dmesh.GetNFaces()
                            },
                            M, N, P);

        static const uint32_t blobsize = 5;
        buf.ComputePointCloudDepthBuffer( blobsize );
        buf.ComputePointCloudVisibilities( blobsize );

        const auto& db = buf.ReturnDepthBuffer();
        image.m_depthMapValues = db;
        image.m_pointCloudDepthVisibilities = buf.ReturnVertexVisibilities();
        image.m_indexValues = buf.ReturnIndexBuffer();

        //Cu::print_matrix_octave(image.m_indexValues.data(), 1920, 1080 , "IndexValues");

        double max=-std::numeric_limits<double>::max();
        double min=std::numeric_limits<double>::max();

        for(auto d : db)
        {
            if(d > max) max = d;
            if(d < min && d > 0) min = d;
        }

        for( uint32_t i=0; i<M; ++i)
        for( uint32_t j=0; j<N; ++j)
        {
            if(db[ M*j + i ] > 0)
            {
                double grayScale = (db[ M*j + i ] - min)/(max - min);
                uint32_t igray = static_cast<uint32_t>(std::floor(255*grayScale));
                QColor qgray( igray,igray,igray );
                image.m_depthMap.setPixelColor( QPoint(j,i ), qgray);
            }
            else
            {
                QColor qundef( 200, 200, 230 );
                image.m_depthMap.setPixelColor( QPoint(j,i ), qundef);
            }
        }

        image.m_depthMap = image.m_depthMap.scaled(image.m_image.size());
        image.m_depthMapLoaded = true;

        ++itp;
    }
    m_depthMapComputed = true;
}

void ImageViewer::ShowAnnotations()
{
    m_msg.print() << "Displaying annotations" << std::endl;
    m_viewMode = Annotation;
    RefreshImage();
}

void ImageViewer::ShowSegmentation()
{
    m_msg.print() << "Displaying segmentations" << std::endl;
    m_viewMode = Segmentation;
    RefreshImage();
}

void ImageViewer::ShowImage()
{
    m_msg.print() << "Displaying image" << std::endl;
    m_viewMode = Image;
    RefreshImage();
}

void ImageViewer::ShowBoth()
{
    m_msg.print() << "Displaying both" << std::endl;
    m_viewMode = Both;
    RefreshImage();
}

void ImageViewer::ShowDepthMap()
{
    m_msg.print() << "Displaying depth map" << std::endl;
    m_viewMode = DepthMap;
    RefreshImage();
}
