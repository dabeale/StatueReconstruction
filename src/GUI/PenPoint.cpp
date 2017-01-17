#include "PenPoint.h"

namespace SegmentationGUI
{
PenPoint::PenPoint() :
    m_pen(),
    m_point()
{

}

PenPoint::PenPoint(Qt::GlobalColor colour, uint32_t alpha, QPoint point ) :
    m_point(point)
{
    QColor qc(colour);
    qc.setAlpha(alpha);
    m_pen = QPen(qc);
}
}
