#ifndef PENPOINT_H
#define PENPOINT_H

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

#include <QPen>
#include <QPoint>

namespace SegmentationGUI
{
/**
 * @brief The PenPoint struct
 * Contains  details of the pen, such as colour and size,
 * also the position at which it was drawn on the qobject.
 */
struct PenPoint
{
public:
    /**
     * @brief PenPoint
     * The default constructor
     */
    PenPoint();

    /**
     * @brief PenPoint
     * Construct the penpoint with  a colour, an alpha value and its location
     * @param colour The colour of the drawn point
     * @param alpha The alpha value
     * @param point The location of the point
     */
    PenPoint( Qt::GlobalColor colour, uint32_t alpha, QPoint point );
    QPen m_pen;     ///< The pen information
    QPoint m_point; ///< The location of the drawn point
};
}
#endif // PENPOINT_H
