#ifndef MAINWINDOW_H
#define MAINWINDOW_H

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

#include <QMainWindow>
#include <QMouseEvent>
#include <QResizeEvent>
#include <QComboBox>
#include <QToolButton>
#include <QToolBar>
#include <QFileDialog>
#include <QStackedWidget>
#include <QInputDialog>
#include <QSlider>
#include <QVector>
#include <QMenu>
#include <QMenuBar>
#include <QHBoxLayout>
#include <QDir>
#include <QAction>
#include <QGroupBox>
#include <QRadioButton>

#include <iostream>
#include <iomanip>

#include "Matrix.h"

#include "ImageViewer.hpp"

/**
 * @brief The MainWindow class
 * This class provides the main graphical user interface.
 */
class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private slots:
    bool eventFilter(QObject *obj, QEvent *event); ///< Event filter signal
    void ComboChanged(const QString&); ///< Notify when the combo box has changed
    void ResizeWindow(); ///< Resize the window so that it fits the (scaled) image
    void SetSliderValues(); ///< Set the slider values correctly for the video
    void IncrementSlider();
    void DecrementSlider();

private:
    QToolBar *m_mainToolBar;
    QToolBar *m_videoToolBar;

    QGroupBox *m_groupBox;
    QRadioButton *m_radioEdge;
    QRadioButton *m_radioForeground;
    QRadioButton *m_radioBackground;

    QGroupBox *m_viewGroupBox;
    QRadioButton *m_radioImage;
    QRadioButton *m_radioAnnotations;
    QRadioButton *m_radioSegmentation;
    QRadioButton *m_radioBoth;
    QRadioButton *m_radioDepthMap;

    QToolButton *m_CycleLeft;
    QToolButton *m_CycleRight;
    QToolButton *m_Clear;
    QToolButton *m_ComputeModel;
    QToolButton *m_RunMRF;
    QToolButton *m_ComputeDepthMaps;
    QSlider *m_widthSlider;
    QSlider *m_frameSlider;
    QLabel *m_frameLabel;

    QHBoxLayout *m_cl;

    QWidget *m_window;

    QMenu *m_fileMenu;
    QAction *m_openAction;
    QAction *m_dumpToOctaveAction;
    QAction *m_saveSegmentation;
    QAction *m_loadMesh;
    QAction *m_loadCameras;

    QMenu *m_parameterMenu;
    QAction *m_MRFParameters;

    ImageViewer *m_imageViewer; ///< An extension of a QLabel which can also deal with painting points

    bool m_showColors; ///< If true show colors

    QString m_startDirectory;

    double m_SampleSigma;
};


#endif // MAINWINDOW_H
