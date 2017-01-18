#ifndef IMAGEVIEWER_H
#define IMAGEVIEWER_H

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

#include "PenPoint.h"
#include "ImageInfo.h"
#include "SegmentationModel.h"
#include "Mesh.h"
#include "Utils.h"
#include "ZBuffer.h"
#include "MRFParametersWindow.h"

#include <QLabel>
#include <QPainter>
#include <QMouseEvent>
#include <QFileDialog>
#include <QInputDialog>
#include <iostream>
#include <cmath>
#include <sstream>

namespace SegmentationGUI
{
/**
 * @brief The ImageViewer class.
 * Extends QLabel into an object which can deal with images and a collection of points drawn on it.
 */
class ImageViewer : public QLabel
{
    Q_OBJECT

public:
    typedef std::vector< PenPoint > PointBuffer; ///< A vector of drawn PenPoints
    typedef std::vector< ImageInfo > ImageList;  ///< A vector of ImageInfos

    /**
     * @brief An explicit constructor for the ImageViewer.
     * @param startDirectory This directory indicates where the file input dialoges should start.
     * @param parent the parent QWidget
     */
    ImageViewer(QString& startDirectory, QWidget *parent = 0);

    void SetPaintInitialised(const bool initialised);    ///< \brief Allow points to be drawn on the ImageViewer.
    void SetColor(Qt::GlobalColor color);                ///< \brief Set the colour of the brush
    void SetSize(uint32_t size);                         ///< \brief Set the size of the brush
    void RepaintPointBuffer();                           ///< \brief Repaint all of the points in the point buffer
    bool AddImages(const QStringList& fileNames);        ///< \brief Add images to the ImageList
    void InitialiseBuffer();                             ///< \brief Initialise the object
    void RefreshImage();                                 ///< \brief Repaint
    uint32_t GetNumberOfImages() const;                  ///< \brief Return the number of images in the list.
    ImageList& GetImageInfos();                          ///< \brief Return a reference to the ImageList

    void GetDimensions( uint32_t& height, uint32_t& width ) const; ///< \brief Get the file width and height

public slots:
    void OpenFile();                ///< \brief Open a file from disk
    void DumpToOctave();            ///< \brief Dump the user selected points to disk (in octave format)
    void CycleLeft();               ///< \brief Select the neighbouring image to the left for viewing
    void CycleRight();              ///< \brief Select the neighbouring image to the right for viewing
    void SetPosition(int frame);    ///< \brief Set the viewed image to 'frame'
    void Clear();                   ///< \brief Clear the user annotations from the current image
    void ClearAll();                ///< \brief Clear all of the user annotations.
    void Select();                  ///< \brief Set the selection tool to the 'edge' tool
    void SelectForeground();        ///< \brief Set the selection tool to be the 'foreground' tool.
    void SelectBackground();        ///< \brief Set the selection tool to be the 'background' tool
    void ComputeModel();            ///< \brief Compute the statistical model given all of the annotations
    void RunMRF();                  ///< \brief Run only the MRF calculation, without computing the model.
    void ComputeDepthMaps();        ///< \brief Compute the depth maps (requires cameras and point cloud to be loaded)
    void ShowAnnotations();         ///< \brief Selects the 'Annotations' view mode.
    void ShowSegmentation();        ///< \brief Selects the 'Segmentation' view mode
    void ShowImage();               ///< \brief Selects the 'Image' view mode
    void ShowBoth();                ///< \brief Selects the 'Segmentations and Annotations' view mode
    void ShowDepthMap();            ///< \brief Selects the 'Depth Map' view mode
    void SaveSegmentation();        ///< \brief Save the segmentations to disk
    void LoadMesh();                ///< \brief Load the point cloud / mesh from file
    void LoadCameras();             ///< \brief Load the cameras from a camera file
    void ShowParametersWindow();    ///< \brief Show the MRF parameters window.

private slots:
    void mousePressEvent(QMouseEvent *e);
    void mouseReleaseEvent(QMouseEvent *e);
    void mouseMoveEvent(QMouseEvent *e);
    void paintEvent(QPaintEvent *e);
    void UpdateMRFParameters();     ///< \brief Update the MRF parameters from the window

signals:
    void ImagesLoaded();  ///< Emitted when the images have been loaded
    void ModelComputed(); ///< Emitted when the model has been computed

private:
    void DrawAnnotation(QPainter &painter); ///< \brief Draw the user annotations on to the current image

    Mesh m_3dmesh;                                  ///< A mesh or point cloud
    std::vector<Math::Matrix> m_ProjectionMatrices; ///< The projection matrices for each camera

    int m_x;         ///< Becomes m_mousex when the mouse is outside of a radius of m_x
    int m_y;         ///< Becomes m_mousey when the mouse is outside of a radius of m_y
    int m_mousex;    ///< Current dragged mouse x position
    int m_mousey;    ///< Current dragged mouse y position
    uint32_t m_size; ///< Size of the pointer

    bool m_latch;            ///< Only allow a blob to be drawn once within a radius of m_x, m_y
    bool m_paintInitialised; ///< If true the ImageViewer will start painting to the canvas

    /**
     * @brief An enumeration of view modes. It identifies what is
     * shown in the image viewer and also which user annotations are shown.
     */
    enum ViewMode
    {
        Image,
        Segmentation,
        Annotation,
        Both,
        DepthMap
    };

    ViewMode m_viewMode;     ///< The currently selected view mode.

    Qt::GlobalColor m_color; ///< Color of cursor

    ImageList m_imageInfos;        ///< A list of images and other relevant info
    ImageList::iterator m_current; ///< Pointer to the currently viewed image

    /**
     * @brief An enumeration of user selection tools.
     */
    enum SketchType
    {
        Foreground,
        Background,
        Edge
    };
    SketchType m_sketchtype; ///< The currently selected sketch type

    void SetSketchType( SketchType type ); ///< \brief Set the sketch type

    template<class List, class Iterator>
    void CycleIteratorLeft(Iterator& iterator, List& imageList); ///< Cycle the image left
    template<class List, class Iterator>
    void CycleIteratorRight(Iterator& iterator, List &imageList); ///< Cycle the image right

    QString& m_startDirectory; ///< The directory to initialise the input dialogues.

    Stream::Message m_msg;     ///< The message streaming interface

    bool m_PointsDrawn;        ///< True if any points have been drawn on the interface
    bool m_ModelComputed;      ///< True when the model has been computed
    bool m_depthMapComputed;   ///< True if the depth maps have been computed

    SegmentationModel m_sgm;    ///< The segmentation model

    MRFParametersWindow *m_mrfparamswindow; ///< The MRF parameters window
};
}

#endif // IMAGEVIEWER_H
