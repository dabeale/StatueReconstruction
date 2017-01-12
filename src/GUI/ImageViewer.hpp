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


/**
 * @brief The ImageViewer class
 * Extends QLabel into an object which can deal with
 * images and a collection of points drawn on it.
 */
class ImageViewer : public QLabel
{
    Q_OBJECT

public:
    typedef std::vector< PenPoint > PointBuffer; ///< A vector of drawn point
    typedef std::vector< ImageInfo > ImageList; ///< Typedef

    ImageViewer(QString& startDirectory, QWidget *parent = 0);

    void SetPaintInitialised(const bool initialised); ///< Allow points to be drawn on the ImageViewer.
    void SetColor(Qt::GlobalColor color); ///< Set the colour of the brush
    void SetSize(uint32_t size); ///< Set the size of the brush
    void RepaintPointBuffer(); ///< Repaint all of the points in the point buffer
    bool AddImages(const QStringList& fileNames);
    void InitialiseBuffer();
    void RefreshImage();
    uint32_t GetNumberOfImages() const;
    ImageList& GetImageInfos();

    void GetDimensions( uint32_t& height, uint32_t& width ) const;

public slots:
    void OpenFile();
    void DumpToOctave();
    void SaveEdges();
    void CycleLeft();
    void CycleRight();
    void SetPosition(int frame);
    void Clear();
    void ClearAll();
    void Select();
    void SelectForeground();
    void SelectBackground();
    void ComputeModel();
    void RunMRF();
    void ComputeDepthMaps();
    void ShowAnnotations();
    void ShowSegmentation();
    void ShowImage();
    void ShowBoth();
    void ShowDepthMap();
    void SaveSegmentation();
    void LoadMesh();
    void LoadCameras();
    void ShowParametersWindow();

private slots:
    void mousePressEvent(QMouseEvent *e);
    void mouseReleaseEvent(QMouseEvent *e);
    void mouseMoveEvent(QMouseEvent *e);
    void paintEvent(QPaintEvent *e);
    void UpdateMRFParameters();

signals:
    void ResizeEvent();
    void ImagesLoaded();
    void ModelComputed();

private:
    void DrawAnnotation(QPainter &painter); ///< Draw the annotations

    Mesh m_3dmesh; ///< A mesh or point cloud
    std::vector<Math::Matrix> m_ProjectionMatrices; ///< The projection matrices for each camera

    int m_x;  ///< Current dragged mouse x position
    int m_y;  ///< Current dragged mouse y position
    int m_mousex;
    int m_mousey;
    uint32_t m_size; ///< Size of the pointer

    bool m_latch; ///< Only allow a point to be drawn once
    bool m_paintInitialised; ///< If true the ImageViewer will start painting to the canvas

    enum ViewMode
    {
        Image,
        Segmentation,
        Annotation,
        Both,
        DepthMap
    };

    ViewMode m_viewMode; ///< Show segmentation

    Qt::GlobalColor m_color; ///< Color of cursor

    ImageList m_imageInfos; ///< A list of images and other relevant info
    ImageList::iterator m_current; ///< Pointer to the currently viewed image

    enum SketchType
    {
        Foreground,
        Background,
        Edge
    };
    SketchType m_sketchtype;

    void SetSketchType( SketchType type );

    template<class List, class Iterator>
    void CycleIteratorLeft(Iterator& iterator, List& imageList); ///< Cycle the image left
    template<class List, class Iterator>
    void CycleIteratorRight(Iterator& iterator, List &imageList); ///< Cycle the image right

    QString& m_startDirectory;

    Stream::Message m_msg;

    bool m_PointsDrawn;
    bool m_ModelComputed;
    bool m_depthMapComputed;

    SegmentationModel m_sgm;

    MRFParametersWindow *m_mrfparamswindow;
};

#endif // IMAGEVIEWER_H
