#include "MainWindow.h"


// Alows the QVector<QPoint> to be registered and used in a signal, when combined with qRegisterMetaType
Q_DECLARE_METATYPE ( QVector<QPoint> )

namespace SegmentationGUI
{
MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
//    ui(),
    m_mainToolBar(new QToolBar(this)),
    m_videoToolBar(new QToolBar(this)),
    m_groupBox(new QGroupBox(tr("&Selector"), this)),
    m_radioEdge(new QRadioButton(tr("&Edge"), this)),
    m_radioForeground(new QRadioButton(tr("&Foreground"), this)),
    m_radioBackground(new QRadioButton(tr("&Background"), this)),
    m_viewGroupBox(new QGroupBox(tr("&View"), this)),
    m_radioImage(new QRadioButton(tr("&Image"), this)),
    m_radioAnnotations(new QRadioButton(tr("&Annotations"), this)),
    m_radioSegmentation(new QRadioButton(tr("&Segmentation"), this)),
    m_radioBoth(new QRadioButton(tr("&Both"), this)),
    m_radioDepthMap( new QRadioButton(tr("&Depth Map"), this)),
    m_CycleLeft(new QToolButton(this)),
    m_CycleRight(new QToolButton(this)),
    m_Clear(new QToolButton(this)),
    m_ComputeModel(new QToolButton(this)),
    m_RunMRF(new QToolButton(this)),
    m_ComputeDepthMaps( new QToolButton(this)),
    m_widthSlider(new QSlider(Qt::Horizontal,this)),
    m_frameSlider(new QSlider(Qt::Horizontal, this)),
    m_frameLabel(new QLabel(tr("0"),this)),
    m_window(new QWidget(this)),
    m_fileMenu(new QMenu(tr("File"), this)),

    m_openAction(new QAction(tr("&Open"), this)),
    m_dumpToOctaveAction(new QAction(tr("&Octave Dump"), this)),
    m_saveSegmentation( new QAction(tr("&Save Segmentation"), this) ),
    m_loadMesh( new QAction(tr("&Load mesh") ,this)),
    m_loadCameras( new QAction(tr("&Load Cameras"), this)),
    m_parameterMenu(new QMenu(tr("Parameters"),this)),
    m_MRFParameters(new QAction(tr("&MRF Parameters"), this)),
    m_imageViewer(new ImageViewer(m_startDirectory, this)),
    m_showColors(true),
    m_SampleSigma(5.0)
{
    menuBar()->addMenu( m_fileMenu );
    m_fileMenu->addAction( m_openAction );
    m_fileMenu->addAction( m_dumpToOctaveAction );
    m_fileMenu->addAction( m_saveSegmentation );
    m_fileMenu->addAction( m_loadMesh );
    m_fileMenu->addAction( m_loadCameras );

    menuBar()->addMenu( m_parameterMenu );
    m_parameterMenu->addAction(m_MRFParameters);

    m_startDirectory = QDir::homePath();
    addToolBar(Qt::LeftToolBarArea, m_mainToolBar);
    addToolBar(Qt::BottomToolBarArea, m_videoToolBar);

    m_imageViewer->setMinimumSize(480,320); // Allows downscaling ( not obvious )


    QHBoxLayout* cl = new QHBoxLayout;
    cl->addWidget(m_imageViewer);

    m_imageViewer->setScaledContents(true);
    m_imageViewer->setMouseTracking(true);
    m_window->setLayout(cl);
    m_window->setMinimumSize(QSize(460, 380));
    setCentralWidget(m_window);

    //qApp->installEventFilter(this);
    m_CycleLeft->setText( "<-" );
    m_CycleRight->setText( "->" );
    m_Clear->setText("Clear");
    m_ComputeModel->setText("Compute Model");
    m_RunMRF->setText("Run MRF");
    m_ComputeDepthMaps->setText("Compute Depth Maps");

    // Connect listener loaded points to the load points slot
    qRegisterMetaType< QVector<QPoint> >("QVector<PenPoint>");

    m_widthSlider->setMinimum(5);
    m_widthSlider->setMaximum(60);
    m_widthSlider->installEventFilter(this);

    // Set up the main toolbar - for selection tools
    m_mainToolBar->addWidget(m_groupBox);
    QVBoxLayout *vbox = new QVBoxLayout;
    vbox->addWidget(m_radioEdge);
    vbox->addWidget(m_radioForeground);
    vbox->addWidget(m_radioBackground);
    m_groupBox->setLayout(vbox);
    m_radioEdge->setChecked(true);

    m_mainToolBar->addWidget(m_viewGroupBox);
    QVBoxLayout *vboxview = new QVBoxLayout;
    vboxview->addWidget( m_radioImage );
    vboxview->addWidget( m_radioAnnotations );
    vboxview->addWidget( m_radioSegmentation );
    vboxview->addWidget( m_radioBoth );
    vboxview->addWidget( m_radioDepthMap );
    m_viewGroupBox->setLayout( vboxview );
    m_radioAnnotations->setChecked(true);

    m_mainToolBar->addSeparator();
    m_mainToolBar->addWidget(m_Clear);
    m_mainToolBar->addSeparator();
    m_mainToolBar->addWidget(m_ComputeModel);
    m_mainToolBar->addWidget(m_RunMRF);
    m_mainToolBar->addWidget(m_ComputeDepthMaps);
    m_mainToolBar->addSeparator();

    m_mainToolBar->addWidget(new QLabel("Tool width:",this));
    m_mainToolBar->addWidget(m_widthSlider);

    m_videoToolBar->addWidget(m_CycleLeft);
    m_videoToolBar->addWidget(m_frameLabel);
    m_videoToolBar->addWidget(m_frameSlider);
    m_videoToolBar->addWidget(m_CycleRight);

    connect(m_radioForeground, SIGNAL(clicked()), m_imageViewer, SLOT(SelectForeground()));
    connect(m_radioBackground, SIGNAL(clicked()), m_imageViewer, SLOT(SelectBackground()));
    connect(m_radioEdge, SIGNAL(clicked()), m_imageViewer, SLOT(Select()));
    connect(m_CycleLeft, SIGNAL(clicked()), m_imageViewer, SLOT(CycleLeft()));
    connect(m_CycleLeft, SIGNAL(clicked()), this, SLOT(DecrementSlider()));
    connect(m_CycleRight, SIGNAL(clicked()), m_imageViewer, SLOT(CycleRight()));
    connect(m_CycleRight, SIGNAL(clicked()), this, SLOT(IncrementSlider()));
    connect(m_Clear, SIGNAL(clicked()), m_imageViewer, SLOT(Clear()));
    connect(m_ComputeModel, SIGNAL(clicked()), m_imageViewer, SLOT(ComputeModel()));
    connect(m_RunMRF, SIGNAL(clicked()), m_imageViewer, SLOT(RunMRF()));
    connect(m_ComputeDepthMaps, SIGNAL(clicked()), m_imageViewer, SLOT(ComputeDepthMaps()));
    connect(m_radioImage, SIGNAL(clicked()), m_imageViewer, SLOT(ShowImage()));
    connect(m_radioAnnotations, SIGNAL(clicked()), m_imageViewer, SLOT(ShowAnnotations()));
    connect(m_radioSegmentation, SIGNAL(clicked()), m_imageViewer, SLOT(ShowSegmentation()));
    connect(m_radioBoth, SIGNAL(clicked()), m_imageViewer, SLOT(ShowBoth()));
    connect(m_radioDepthMap, SIGNAL(clicked()), m_imageViewer, SLOT(ShowDepthMap()));

    connect(m_openAction, SIGNAL(triggered()), m_imageViewer, SLOT(OpenFile()));
    connect(m_dumpToOctaveAction, SIGNAL(triggered()), m_imageViewer, SLOT(DumpToOctave()));
    connect(m_saveSegmentation, SIGNAL(triggered()), m_imageViewer, SLOT(SaveSegmentation()));
    connect(m_loadMesh, SIGNAL(triggered()), m_imageViewer, SLOT(LoadMesh()));
    connect(m_loadCameras, SIGNAL(triggered()), m_imageViewer, SLOT(LoadCameras()));

    connect(m_imageViewer, SIGNAL(ImagesLoaded()), this, SLOT(ResizeWindow()));
    connect(m_imageViewer, SIGNAL(ImagesLoaded()), this, SLOT(SetSliderValues()));

    connect(m_frameSlider, SIGNAL(sliderMoved(int)), m_imageViewer, SLOT(SetPosition(int)));
    connect(m_frameSlider, SIGNAL(sliderMoved(int)), m_frameLabel, SLOT(setNum(int)));

     connect(m_MRFParameters, SIGNAL(triggered()), m_imageViewer, SLOT(ShowParametersWindow()));

     connect(m_imageViewer, SIGNAL(ModelComputed()), m_radioSegmentation, SLOT(click()));
}

MainWindow::~MainWindow()
{

}

void MainWindow::ResizeWindow()
{
    uint32_t height, width;
    m_imageViewer->GetDimensions(height, width);
    double heightd = (420.0 / static_cast<double>(width)) * height;
    m_window->setMinimumSize(QSize(420, static_cast<uint32_t>(heightd)));
}

void MainWindow::SetSliderValues()
{
    m_frameSlider->setMinimum(0);
    m_frameSlider->setMaximum(m_imageViewer->GetNumberOfImages()-1);
}

void MainWindow::IncrementSlider()
{
    m_frameSlider->setValue(m_frameSlider->value()+1);
    m_frameLabel->setText( QString("%1").arg(m_frameSlider->value()) );
}

void MainWindow::DecrementSlider()
{
    m_frameSlider->setValue(m_frameSlider->value()-1);
    m_frameLabel->setText( QString("%1").arg(m_frameSlider->value()) );
}

void MainWindow::ComboChanged( const QString& )
{
    m_imageViewer->SetPaintInitialised(m_showColors);
    m_imageViewer->RefreshImage();
}

bool MainWindow::eventFilter(QObject *obj, QEvent *event)
 {
    if( obj == this )
    {
        // nothing for now
    }
    else if(obj == m_widthSlider )
    {
        if(event->type() == QEvent::MouseButtonRelease)
        {
            m_imageViewer->SetSize(m_widthSlider->value());
            m_SampleSigma = static_cast<double>(m_widthSlider->value()) / 5;
            //m_tracks->SetDrawWidth(m_SampleSigma);
            std::cout << "Draw width set to: " << m_SampleSigma << std::endl;
        }
    }

    return QMainWindow::eventFilter(obj, event);
}
}
