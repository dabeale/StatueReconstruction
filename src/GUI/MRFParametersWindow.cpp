#include "MRFParametersWindow.h"

MRFParametersWindow::MRFParametersWindow(QWidget *parent ) :
    QMainWindow(parent),
    m_segtypeGroup(new QGroupBox(tr("&Segmentation type"), this)),
    m_radioICM(new QRadioButton(tr("&ICM"), this)),
    m_radioGraphCut(new QRadioButton(tr("&GraphCut"), this)),
    m_alphaparam(new QLineEdit(this)),
    m_sampleN(new QLineEdit(this)),
    m_pairwiseWeight(new QLineEdit(this)),
    m_unaryWeight(new QLineEdit(this)),
    m_mainwindow(new QWidget(this)),
    m_closeButton(new QToolButton(this))
{
    QVBoxLayout *segtypevbox = new QVBoxLayout;
    segtypevbox->addWidget(m_radioICM);
    segtypevbox->addWidget(m_radioGraphCut);
    m_radioGraphCut->setChecked(true);
    m_segtypeGroup->setLayout(segtypevbox);

    m_alphaparam->setValidator( new QDoubleValidator(0, 100, 2, this) );
    m_sampleN->setValidator( new QIntValidator(0, 100, this) );
    m_pairwiseWeight->setValidator( new QDoubleValidator(0, 100, 2, this) );
    m_unaryWeight->setValidator( new QDoubleValidator(0, 100, 2, this) );
    m_alphaparam->setText("100.0");
    m_sampleN->setText("5");
    m_pairwiseWeight->setText("10.0");
    m_unaryWeight->setText("0.1");

    QVBoxLayout *mainLayout = new QVBoxLayout;
    mainLayout->addWidget(m_segtypeGroup);
    mainLayout->addSpacing(10);
    mainLayout->addWidget( new QLabel(tr("Alpha:"),this));
    mainLayout->addWidget( m_alphaparam );
    mainLayout->addWidget( new QLabel(tr("N Samples"), this) );
    mainLayout->addWidget( m_sampleN );
    mainLayout->addWidget( new QLabel(tr("Pairwise weight"), this) );
    mainLayout->addWidget( m_pairwiseWeight );
    mainLayout->addWidget( new QLabel(tr("Unary weight"), this) );
    mainLayout->addWidget( m_unaryWeight );

    m_closeButton->setText(tr("Close"));
    mainLayout->addWidget(m_closeButton);


    m_mainwindow->setMinimumSize(200, 300);
    m_mainwindow->setLayout(mainLayout);
    setCentralWidget(m_mainwindow);

    connect( m_closeButton, SIGNAL(clicked()), this, SLOT(CloseWindow()) );

    hide();
}

void MRFParametersWindow::CloseWindow()
{
    hide();
    emit CloseButtonPushed();
}

MRFParametersWindow::~MRFParametersWindow()
{

}
