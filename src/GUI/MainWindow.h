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

namespace SegmentationGUI
{
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
    bool eventFilter(QObject *obj, QEvent *event); ///< \brief Event filter signal.
    void ComboChanged(const QString&);             ///< \brief Notify when the combo box has changed
    void ResizeWindow();                           ///< \brief Resize the window so that it fits the (scaled) image
    void SetSliderValues();                        ///< \brief Set the slider values correctly for the video
    void IncrementSlider();                        ///< \brief Increment the slider by one unit
    void DecrementSlider();                        ///< \brief Decrement the slider by one unit

private:
    QToolBar *m_mainToolBar;                ///< The main toolbar containing all of the radio buttons, function buttons and tool width slider
    QToolBar *m_videoToolBar;               ///< The tool bar containing the frame selection slider

    QGroupBox *m_groupBox;                  ///< The user selection radio button group box
    QRadioButton *m_radioEdge;              ///< The 'Edge' selection tool radio button
    QRadioButton *m_radioForeground;        ///< The 'Foreground' selection tool radio button
    QRadioButton *m_radioBackground;        ///< The 'Background' selection radio button

    QGroupBox *m_viewGroupBox;              ///< The view mode radio button group
    QRadioButton *m_radioImage;             ///< 'Image' view mode radio
    QRadioButton *m_radioAnnotations;       ///< 'Annotations' view mode radio
    QRadioButton *m_radioSegmentation;      ///< 'Segmentation' view mode radio
    QRadioButton *m_radioBoth;              ///< 'Both' view mode radio
    QRadioButton *m_radioDepthMap;          ///< 'Depth Map' view mode radio

    QToolButton *m_CycleLeft;               ///< Tool button to cycle the image one to the left
    QToolButton *m_CycleRight;              ///< Tool button to cycle the image one to the right
    QToolButton *m_Clear;                   ///< Tool button to clear annotations from the current image
    QToolButton *m_ComputeModel;            ///< Tool button to compute the segmentation model
    QToolButton *m_RunMRF;                  ///< Tool button to compute only the MRF
    QToolButton *m_ComputeDepthMaps;        ///< Tool button to compute the depthmaps
    QSlider *m_widthSlider;                 ///< A slider identifying the width of the tool in pixels
    QSlider *m_frameSlider;                 ///< A slider identifying the current image number
    QLabel *m_frameLabel;                   ///< A labels specifying the current frame

    QWidget *m_window;                      ///< The main widget, containing everything.

    QMenu *m_fileMenu;                      ///< The file menu
    QAction *m_openAction;                  ///< Open a collection of images from disk
    QAction *m_dumpToOctaveAction;          ///< Dump all of the user annotations to octave files on disk
    QAction *m_saveSegmentation;            ///< Save the generated segmentations
    QAction *m_loadMesh;                    ///< Load a ply file from disk
    QAction *m_loadCameras;                 ///< Load cameras from disk

    QMenu *m_parameterMenu;                 ///< A menu relating to parameter input
    QAction *m_MRFParameters;               ///< Show the MRF parameters window

    ImageViewer *m_imageViewer;             ///< The image viewer

    bool m_showColors;                      ///< If true show the user annotations

    QString m_startDirectory;               ///< The start directory for the file input dialogues

    double m_SampleSigma;                   ///< The width of the tool in pixels
};
}


#endif // MAINWINDOW_H
