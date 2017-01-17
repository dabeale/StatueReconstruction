#include "marchingcubes.h"

namespace MC
{

void WritePly( std::string filename,
               const std::vector<Pointd3>& vertices,
               const std::vector<Pointi3>& faces)
{
    std::ofstream ofs(filename, std::ofstream::out);
    ofs << "ply" << std::endl;
    ofs << "format ascii 1.0" << std::endl;
    ofs << "comment Daniel Beale (2015)" << std::endl;
    ofs << "element vertex " << vertices.size() << std::endl;
    ofs << "property float x" << std::endl;
    ofs << "property float y" << std::endl;
    ofs << "property float z" << std::endl;
    ofs << "element face " << faces.size() << std::endl;
    ofs << "property list uchar int vertex_indices" << std::endl;
    ofs << "end_header" << std::endl;
    for(const auto& v : vertices)
    {
        ofs << v.p[0] << " " << v.p[1] << " " << v.p[2] << std::endl;
    }
    for(const auto& f : faces)
    {
        ofs << 3 << " " << f.p[0] << " " << f.p[1] << " " << f.p[2] << std::endl;
    }
    ofs.close();
}

Pointd3 VertexInterp(const double& isolevel, const Pointd3& p1, const Pointd3& p2, const double& valp1, const double& valp2)
{
    if (std::abs(isolevel-valp1) < 0.00001 || std::abs(valp1-valp2) < 0.00001)
        return(p1);
    if (std::abs(isolevel-valp2) < 0.00001)
        return(p2);

    double mu = (isolevel - valp1) / (valp2 - valp1);
    Pointd3 p;
    for(uint8_t i=0;i<3;i++)
    {
        p.p[i] = p1.p[i] + mu * (p2.p[i] - p1.p[i]);
    }

    return(p);
}

MarchingCubes::MarchingCubes( const double isoval , const double npoints, const double max[3], const double min[3]) :
    m_isoval(isoval), m_npoints(npoints),
    m_pointList(0), m_faceList(0), m_map(),
    m_renderToScreen(false),
    m_boolFunc([](const double*){return true;}),
    m_doboolfunc(false)
{
    // The faces are not necessary since the vertices are ordered in to triangles
    // they are included so that saving to ply files is easy.

    // Find the extremal values on the grid
    //ComputeExtremalValues();
    for(uint8_t i=0; i<3; i++)
    {
        m_max[i] = max[i];
        m_min[i] = min[i];
    }

    for(uint8_t i=0;i<3;i++)
    {
        m_step[i] = std::abs(m_max[i] - m_min[i])/static_cast<double>(m_npoints);
    }
}

MarchingCubes::~MarchingCubes()
{
}

void MarchingCubes::EnableBoolFunc( const std::function<bool(const double*)>& fun )
{
    m_boolFunc = fun;
    m_doboolfunc = true;
}

void MarchingCubes::SetRenderToScreen(bool render)
{
    m_renderToScreen = render;
}

void MarchingCubes::WriteToPly( const std::string& filename )
{
    const auto& f = GetFaces();
    const auto& v = GetPoints();
    WritePly(filename, v, f);
}

#ifdef RENDER_USING_GLUT
inline void cross(const double* a, const double* b, double cr[3] )
{
  cr[0] = a[1]*b[2] - a[2]*b[1];
  cr[1] = -(a[0]*b[2] - a[2]*b[0]);
  cr[2] = a[0]*b[1] - a[1]*b[0];
}
#endif

void MarchingCubes::PolygoniseCube(
        const GridCell &g, std::vector<Pointd3>& pointList,
        std::vector<Pointi3>& faceList, uint32_t cindex)
{
    /* Find the vertices where the surface intersects the cube */
    int map[12];
    //Create the triangles
    if(!m_doboolfunc || (m_doboolfunc && m_boolFunc(&g.p[0].p[0])))
    {
        uint32_t vertind = pointList.size();
        uint32_t vit=0;
        for(j=0; j<12; j++)
        {
            if (edgeTable[cindex] & pw[j])
            {
                Pointd3 pd3 = VertexInterp(m_isoval,g.p[wq[j]],g.p[ww[j]],g.val[wq[j]],g.val[ww[j]]);
                pointList.push_back(pd3);
                map[j]=vit;
                vit++;
            }
        }

        for (i=0;triTable[cindex][i]!=-1;i+=3)
        {
            Pointi3 pi3;
            pi3.p[0] = vertind + map[triTable[cindex][i  ]];
            pi3.p[1] = vertind + map[triTable[cindex][i+1]];
            pi3.p[2] = vertind + map[triTable[cindex][i+2]];
            faceList.push_back(pi3);
#ifdef RENDER_USING_GLUT //this could be much faster!!
            if(m_renderToScreen)
            {
                double a[3], b[3], cr[3], nm;
                int k,l;
                for( k=0; k<3; k++)
                {
                    for( l=0;l<3;l++)
                    {
                        a[l] = pointList[pi3.p[(k+1)%3]].p[l]-pointList[pi3.p[k]].p[l];
                        b[l] = pointList[pi3.p[(k+2)%3]].p[l]-pointList[pi3.p[k]].p[l];
                    }
                    cross(a,b,cr);
                    nm = std::sqrt(cr[0]*cr[0] + cr[1]*cr[1] + cr[2]*cr[2]);
                    glNormal3f(cr[0]/nm,cr[1]/nm,cr[2]/nm);

                    glVertex3f(pointList[pi3.p[k]].p[0],
                               pointList[pi3.p[k]].p[1],
                               pointList[pi3.p[k]].p[2]);
                }
            }
#endif
        }
    }
}


void MarchingCubes::PolygoniseCube(
        const GridCell &g,
        std::vector<Pointd3>::iterator& vertit,
        std::vector<Pointi3>::iterator &faceit,
        uint32_t& vertindex)
{
    /* Find the vertices where the surface intersects the cube */

    //Create the triangles
    if(!m_doboolfunc || (m_doboolfunc && m_boolFunc(&g.p[0].p[0])))
    {
        vi=0;
        for(j=0; j<12; j++)
        {
            if (edgeTable[cubeindex] & pw[j])
            {
                Pointd3 pd3 = VertexInterp(m_isoval,g.p[wq[j]],g.p[ww[j]],g.val[wq[j]],g.val[ww[j]]);
                *vertit = pd3;
                ++vertit;
                m_map[j]=vi;
                vi++;
            }
        }

        for (i=0;triTable[cubeindex][i]!=-1;i+=3) {
            faceit->p[0] = vertindex + m_map[triTable[cubeindex][i  ]];
            faceit->p[1] = vertindex + m_map[triTable[cubeindex][i+1]];
            faceit->p[2] = vertindex + m_map[triTable[cubeindex][i+2]];
#ifdef RENDER_USING_GLUT //this could be much faster!!
            if(m_renderToScreen)
            {
                double a[3], b[3], cr[3], nm;
                int k,l;
                for( k=0; k<3; k++)
                {
                    for( l=0;l<3;l++)
                    {
                        a[l] = m_pointList.at(faceit->p[(k+1)%3]).p[l]-m_pointList.at(faceit->p[k]).p[l];
                        b[l] = m_pointList.at(faceit->p[(k+2)%3]).p[l]-m_pointList.at(faceit->p[k]).p[l];
                    }
                    cross(a,b,cr);
                    nm = std::sqrt(cr[0]*cr[0] + cr[1]*cr[1] + cr[2]*cr[2]);
                    glNormal3f(cr[0]/nm,cr[1]/nm,cr[2]/nm);

                    glVertex3f(m_pointList.at(faceit->p[k]).p[0],
                               m_pointList.at(faceit->p[k]).p[1],
                               m_pointList.at(faceit->p[k]).p[2]);
                }
            }
#endif
            ++faceit;
        }
        vertindex+=vi;
    }
}

inline void CreateCube(GridCell& gc, const double point[3], const double step[3], const std::function<double(const double*)>& fun)
{
    gc.p[0].p[0] = point[0];
    gc.p[0].p[1] = point[1];
    gc.p[0].p[2] = point[2];
    gc.val[0] = fun(&(gc.p[0].p[0]));
    gc.p[1].p[0] = point[0] + step[0];
    gc.p[1].p[1] = point[1];
    gc.p[1].p[2] = point[2];
    gc.val[1] = fun(&(gc.p[1].p[0]));
    gc.p[2].p[0] = point[0] + step[0];
    gc.p[2].p[1] = point[1] + step[1];
    gc.p[2].p[2] = point[2];
    gc.val[2] = fun(&(gc.p[2].p[0]));
    gc.p[3].p[0] = point[0];
    gc.p[3].p[1] = point[1] + step[1];
    gc.p[3].p[2] = point[2];
    gc.val[3] = fun(&(gc.p[3].p[0]));
    gc.p[4].p[0] = point[0];
    gc.p[4].p[1] = point[1];
    gc.p[4].p[2] = point[2] + step[2];
    gc.val[4] = fun(&(gc.p[4].p[0]));
    gc.p[5].p[0] = point[0] + step[0];
    gc.p[5].p[1] = point[1];
    gc.p[5].p[2] = point[2] + step[2];
    gc.val[5] = fun(&(gc.p[5].p[0]));
    gc.p[6].p[0] = point[0] + step[0];
    gc.p[6].p[1] = point[1] + step[1];
    gc.p[6].p[2] = point[2] + step[2];
    gc.val[6] = fun(&(gc.p[6].p[0]));
    gc.p[7].p[0] = point[0];
    gc.p[7].p[1] = point[1] + step[1];
    gc.p[7].p[2] = point[2] + step[2];
    gc.val[7] = fun(&(gc.p[7].p[0]));
}

inline void CreateCube(GridCell& gc, const double point[3], const double step[3], const std::vector<double>& vec, const uint32_t i, const uint32_t j, const uint32_t k, const uint32_t N)
{
    gc.p[0].p[0] = point[0];
    gc.p[0].p[1] = point[1];
    gc.p[0].p[2] = point[2];
    gc.val[0] = vec[i*N*N + j*N + k];
    gc.p[1].p[0] = point[0] + step[0];
    gc.p[1].p[1] = point[1];
    gc.p[1].p[2] = point[2];
    gc.val[1] = vec[(i+1)*N*N + j*N + k];
    gc.p[2].p[0] = point[0] + step[0];
    gc.p[2].p[1] = point[1] + step[1];
    gc.p[2].p[2] = point[2];
    gc.val[2] = vec[(i+1)*N*N + (j+1)*N + k];
    gc.p[3].p[0] = point[0];
    gc.p[3].p[1] = point[1] + step[1];
    gc.p[3].p[2] = point[2];
    gc.val[3] = vec[i*N*N + (j+1)*N + k];
    gc.p[4].p[0] = point[0];
    gc.p[4].p[1] = point[1];
    gc.p[4].p[2] = point[2] + step[2];
    gc.val[4] = vec[i*N*N + j*N + k + 1];
    gc.p[5].p[0] = point[0] + step[0];
    gc.p[5].p[1] = point[1];
    gc.p[5].p[2] = point[2] + step[2];
    gc.val[5] = vec[(i+1)*N*N + j*N + k + 1];
    gc.p[6].p[0] = point[0] + step[0];
    gc.p[6].p[1] = point[1] + step[1];
    gc.p[6].p[2] = point[2] + step[2];
    gc.val[6] = vec[(i+1)*N*N + (j+1)*N + k + 1];
    gc.p[7].p[0] = point[0];
    gc.p[7].p[1] = point[1] + step[1];
    gc.p[7].p[2] = point[2] + step[2];
    gc.val[7] = vec[i*N*N + (j+1)*N + k + 1];
}

inline void CopyGridCell( const GridCell& gc, GridCell& gco)
{
    for(uint32_t ii=0;ii<8;ii++)
    {
        gco.p[ii] = gc.p[ii];
       // gco.n[ii] = gc.n[ii];
        gco.val[ii] = gc.val[ii];
    }
}

void MarchingCubes::KD( const GridCell& gc, uint32_t dim,
                        std::vector<Pointd3>::iterator& vertit,
                        std::vector<Pointi3>::iterator& faceit,
                        uint32_t& vertindex,
                        const std::function<double(const double*)>& fun)
{
    if( dim > 9 )
    {
        cubeindex = 0;
        if (gc.val[0] < m_isoval) cubeindex |= 1;
        if (gc.val[1] < m_isoval) cubeindex |= 2;
        if (gc.val[2] < m_isoval) cubeindex |= 4;
        if (gc.val[3] < m_isoval) cubeindex |= 8;
        if (gc.val[4] < m_isoval) cubeindex |= 16;
        if (gc.val[5] < m_isoval) cubeindex |= 32;
        if (gc.val[6] < m_isoval) cubeindex |= 64;
        if (gc.val[7] < m_isoval) cubeindex |= 128;

        // return if it is a null cube
        if (edgeTable[cubeindex] == 0) return;

        if( gc.p[6].p[0] - gc.p[0].p[0] < m_step[0] &&
            gc.p[6].p[1] - gc.p[0].p[1] < m_step[1] &&
            gc.p[6].p[2] - gc.p[0].p[2] < m_step[2] )
        {
            PolygoniseCube(gc, vertit, faceit, vertindex);
            return;
        }
    }

    // split the current cube accross dim
    double newe;
    GridCell gn[2];
    CopyGridCell( gc, gn[0]);
    CopyGridCell( gc, gn[1]);
    switch(dim%3)
    {
        case 0:
        { // (1,2,5,6) (0,3,4,7)
            newe = 0.5*(gc.p[0].p[0]+gc.p[6].p[0]);
            gn[0].p[1].p[0] = newe; gn[0].p[2].p[0] = newe; gn[0].p[5].p[0] = newe; gn[0].p[6].p[0] = newe;
            gn[0].val[1] = fun(gn[0].p[1].p);gn[0].val[2] = fun(gn[0].p[2].p);
            gn[0].val[5] = fun(gn[0].p[5].p);gn[0].val[6] = fun(gn[0].p[6].p);
            gn[1].p[0].p[0] = newe; gn[1].p[3].p[0] = newe; gn[1].p[4].p[0] = newe; gn[1].p[7].p[0] = newe;
            gn[1].val[0] = fun(gn[1].p[0].p);gn[1].val[3] = fun(gn[1].p[3].p);
            gn[1].val[4] = fun(gn[1].p[4].p);gn[1].val[7] = fun(gn[1].p[7].p);

        } break;
        case 1:
        { // (5,4,1,0) (7,6,3,2)
            newe = 0.5*(gc.p[0].p[1]+gc.p[6].p[1]);
            gn[0].p[5].p[1] = newe; gn[0].p[4].p[1] = newe; gn[0].p[1].p[1] = newe; gn[0].p[0].p[1] = newe;
            gn[0].val[5] = fun(gn[0].p[5].p);gn[0].val[4] = fun(gn[0].p[4].p);
            gn[0].val[1] = fun(gn[0].p[1].p);gn[0].val[0] = fun(gn[0].p[0].p);
            gn[1].p[7].p[1] = newe; gn[1].p[6].p[1] = newe; gn[1].p[3].p[1] = newe; gn[1].p[2].p[1] = newe;
            gn[1].val[7] = fun(gn[1].p[7].p);gn[1].val[6] = fun(gn[1].p[6].p);
            gn[1].val[3] = fun(gn[1].p[3].p);gn[1].val[2] = fun(gn[1].p[2].p);
        } break;
        case 2:
        { // (3,2,1,0) (7,6,5,4)
            newe = 0.5*(gc.p[0].p[2]+gc.p[6].p[2]);
            gn[0].p[3].p[2] = newe; gn[0].p[2].p[2] = newe; gn[0].p[1].p[2] = newe; gn[0].p[0].p[2] = newe;
            gn[0].val[3] = fun(gn[0].p[3].p);gn[0].val[2] = fun(gn[0].p[2].p);
            gn[0].val[1] = fun(gn[0].p[1].p);gn[0].val[0] = fun(gn[0].p[0].p);
            gn[1].p[7].p[2] = newe; gn[1].p[6].p[2] = newe; gn[1].p[5].p[2] = newe; gn[1].p[4].p[2] = newe;
            gn[1].val[7] = fun(gn[1].p[7].p);gn[1].val[6] = fun(gn[1].p[6].p);
            gn[1].val[5] = fun(gn[1].p[5].p);gn[1].val[4] = fun(gn[1].p[4].p);
        } break;

    }
    #pragma omp parallel for
    for(uint32_t kk=0;kk<2;kk++)
    {
        KD(gn[kk], (dim+1), vertit, faceit, vertindex, fun);
    }
}

void MarchingCubes::MarchKD(const std::function<double(const double*)>& fun)
{
    GridCell gc;

    uint32_t vertindex=0;
    m_pointList.resize(m_npoints*m_npoints*m_npoints*12);
    m_faceList.resize(m_npoints*m_npoints*m_npoints*10);

    std::vector<Pointd3>::iterator ipt = m_pointList.begin();
    std::vector<Pointi3>::iterator ifc = m_faceList.begin();

    double bigstep[3];
    for(uint32_t ll=0; ll<3; ll++)
    {
        bigstep[ll] = m_max[ll] - m_min[ll];
    }

    CreateCube(gc, m_min, bigstep, fun);
    KD( gc, 0, ipt, ifc, vertindex, fun);
    m_pointList.resize(ipt-m_pointList.begin());
    m_faceList.resize(ifc-m_faceList.begin());
}

void MarchingCubes::March(const std::function<double(const double*)>& fun)
{
    std::vector<Pointd3> pointList;
    std::vector<Pointi3> faceList;
    pointList.reserve(m_npoints*m_npoints*m_npoints*12);
    faceList.reserve(m_npoints*m_npoints*m_npoints*10);

    std::vector<double> values((m_npoints+1)*(m_npoints+1)*(m_npoints+1),0.0);
#pragma omp parallel for
    for(uint32_t i=0; i<=m_npoints; ++i)
    {
        double val[3];
        val[0] = static_cast<double>(i)*m_step[0] + m_min[0];
        for(uint32_t j=0; j<=m_npoints; ++j)
        for(uint32_t k=0; k<=m_npoints; ++k)
        {
            val[1] = static_cast<double>(j)*m_step[1] + m_min[1];
            val[2] = static_cast<double>(k)*m_step[2] + m_min[2];
            values[i*m_npoints*m_npoints + j*m_npoints + k] = fun(val);
        }
    }

#pragma omp parallel for shared(pointList, faceList)
    for(uint32_t i=0; i<m_npoints; ++i)
    for(uint32_t j=0; j<m_npoints; ++j)
    for(uint32_t k=0; k<m_npoints; ++k)
    {
        double val[3];
        GridCell gc;
        val[0] = static_cast<double>(i)*m_step[0] + m_min[0];
        val[1] = static_cast<double>(j)*m_step[1] + m_min[1];
        val[2] = static_cast<double>(k)*m_step[2] + m_min[2];

        CreateCube(gc, val, m_step, values, i, j, k, m_npoints);
        uint32_t cindex = 0; // cube index
        if (gc.val[0] < m_isoval) cindex |= 1;
        if (gc.val[1] < m_isoval) cindex |= 2;
        if (gc.val[2] < m_isoval) cindex |= 4;
        if (gc.val[3] < m_isoval) cindex |= 8;
        if (gc.val[4] < m_isoval) cindex |= 16;
        if (gc.val[5] < m_isoval) cindex |= 32;
        if (gc.val[6] < m_isoval) cindex |= 64;
        if (gc.val[7] < m_isoval) cindex |= 128;

        //Cube is in/out surface
        if (edgeTable[cindex] == 0) continue;

#pragma omp critical
        PolygoniseCube(gc, pointList, faceList, cindex);
    }

    m_pointList = pointList;
    m_faceList = faceList;
}

const std::vector<Pointd3>& MarchingCubes::GetPoints() const
{
    return m_pointList;
}

const std::vector<Pointi3>& MarchingCubes::GetFaces() const
{
    return m_faceList;
}
}
