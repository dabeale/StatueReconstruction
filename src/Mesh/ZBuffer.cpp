#include "zbuffer.h"

namespace Buffer
{
    ProjectionMatrix::ProjectionMatrix() :
        m_arr(3*4)
    {

    }

    ProjectionMatrix::ProjectionMatrix(double val) :
        m_arr(3*4, val)
    {

    }

    ProjectionMatrix::ProjectionMatrix(const BVector& arr ) :
        m_arr(arr)
    {
        assert(arr.size() == 3*4);
    }

    double ProjectionMatrix::el( const uint32_t i, const uint32_t j) const
    {
        return m_arr[j*3 + i];
    }

    double& ProjectionMatrix::el( const uint32_t i, const uint32_t j)
    {
        return m_arr[j*3 + i];
    }

    BVector ProjectionMatrix::ComputeLocation() const
    {
        double detR = el(0,0)*el(1,1)*el(2,2) +
                      el(0,1)*el(1,2)*el(2,0) +
                      el(0,2)*el(1,0)*el(2,1) -
                      el(0,2)*el(1,1)*el(2,0) -
                      el(0,1)*el(1,0)*el(2,2) -
                      el(0,0)*el(1,2)*el(2,1);
        if(std::abs(detR) < 1e-15)
        {
            std::cerr << "Error::ProjectionMatrix::ComputeLocation :: the matrix is singular" <<std::endl;
            assert(false);
        }
        double a11 = -(el(1,1)*el(2,2) - el(1,2)*el(2,1))/detR;
        double a12 = -(el(0,2)*el(2,1) - el(0,1)*el(2,2))/detR;
        double a13 = -(el(0,1)*el(1,2) - el(0,2)*el(1,1))/detR;
        double a21 = -(el(1,2)*el(2,0) - el(1,0)*el(2,2))/detR;
        double a22 = -(el(0,0)*el(2,2) - el(0,2)*el(2,1))/detR;
        double a23 = -(el(0,2)*el(1,0) - el(0,0)*el(1,2))/detR;
        double a31 = -(el(1,0)*el(2,1) - el(1,1)*el(2,0))/detR;
        double a32 = -(el(0,1)*el(2,0) - el(0,0)*el(2,1))/detR;
        double a33 = -(el(0,0)*el(1,1) - el(0,1)*el(1,0))/detR;

        return  { a11*el(0,3) + a12*el(1,3) + a13*el(2,3),
                  a21*el(0,3) + a22*el(1,3) + a23*el(2,3),
                  a31*el(0,3) + a32*el(1,3) + a33*el(2,3)};
    }

    BVector ProjectionMatrix::operator *(const BVector& X)
    {
        BVector x(3);
        x[0] = m_arr[0]*X[0] +
               m_arr[3]*X[1] +
               m_arr[6]*X[2] +
               m_arr[9]*X[3];
        x[1] = m_arr[1]*X[0] +
               m_arr[4]*X[1] +
               m_arr[7]*X[2] +
               m_arr[10]*X[3];
        x[2]=  m_arr[2]*X[0] +
               m_arr[5]*X[1] +
               m_arr[8]*X[2] +
               m_arr[11]*X[3];
        x[0] /= x[2];
        x[1] /= x[2];
        x[2] = 1;
        return x;
    }

    BVector ProjectionMatrix::operator *(const BVector& X) const
    {
        BVector x(3);
        x[0] = m_arr[0]*X[0] +
               m_arr[3]*X[1] +
               m_arr[6]*X[2] +
               m_arr[9]*X[3];
        x[1] = m_arr[1]*X[0] +
               m_arr[4]*X[1] +
               m_arr[7]*X[2] +
               m_arr[10]*X[3];
        x[2]=  m_arr[2]*X[0] +
               m_arr[5]*X[1] +
               m_arr[8]*X[2] +
               m_arr[11]*X[3];
        x[0] /= x[2];
        x[1] /= x[2];
        x[2] = 1;
        return x;
    }

    ZBuffer::ZBuffer() :
        m_vertices(NULL),
        m_faces(NULL),
        m_Nv(0),
        m_Nf(0),
        m_height(0),
        m_width(0)
    {

    }

    ZBuffer::ZBuffer(const Mesh mesh, const uint32_t height, const uint32_t width , const ProjectionMatrix &P)
        : m_vertices(mesh.vertices),
          m_faces(mesh.faces),
          m_Nv(mesh.Nv),
          m_Nf(mesh.Nf),
          m_height(height),
          m_width(width),
          m_P(P),
          m_db(height*width, -1.0),
          m_ib(height*width, -1.0),
          m_vertexDepth(m_Nv, -1.0),
          m_faceDepth(m_Nf, -1.0),
          m_visibilities(m_Nv, false)
    {

    }

    ZBuffer::~ZBuffer()
    {

    }

    void SetVertices(const double* vertices, const uint32_t* faces, BVector& a, BVector& b, BVector& c, const int& i)
    {
        for(uint32_t k = 0; k < 3; k++)
        {
            a[k] = vertices[3*faces[3*i]+k];
            b[k] = vertices[3*faces[3*i+1]+k];
            c[k] = vertices[3*faces[3*i+2]+k];
        }
    }

    // This could be optimised using element operations
    // there are a large number of copies here
    bool InsideTriangle( const BVector& a, const BVector&b, const BVector& c, const BVector& x )
    {
        double ta=0, tb=0, tc=0;
        for(uint32_t k=0; k<a.size(); ++k)
        {
            ta += (a[k] - x[k])*(b[k] - x[k]);
            tb += (b[k] - x[k])*(c[k] - x[k]);
            tc += (c[k] - x[k])*(a[k] - x[k]);
        }

        return std::abs( std::acos(ta) + std::acos(tb) + std::acos(tc) - 2*M_PI ) < 1e-1;
    }

    void ZBuffer::ComputeVisibilities()
    {
        BVector C = m_P.ComputeLocation();

//#pragma omp parallel for
        for(uint32_t i=0; i<m_Nf; i++)
        {
            BVector a(4,0), b(4,0), c(4,0);
            a[3]=1; b[3]=1; c[3]=1;

            SetVertices(m_vertices,m_faces,a,b,c,i);

            BVector pa = m_P*a;
            BVector pb = m_P*b;
            BVector pc = m_P*c;

            double maxx = std::max( {pa[0], pb[0], pc[0]});
            double maxy = std::max( {pa[1], pb[1], pc[1]} );
            double minx = std::min( {pa[0], pb[0], pc[0]} );
            double miny = std::min( {pa[1], pb[1], pc[1]} );

            // The whole triangle is inside the depth image
            if(minx>0 && miny >0 && maxx<=m_width && maxy<=m_height)
            {
                //This distance calculation could be interpolated
                double distance = m_faceDepth[i];
                double deptha = GetDepthBufferValue(pa);
                double depthb = GetDepthBufferValue(pb);
                double depthc = GetDepthBufferValue(pc);


                if( std::abs(deptha - distance) < 1e-6)
                {
                    m_visibilities[m_faces[3*i]] = true;
                }
                if( std::abs(depthb - distance) < 1e-6)
                {
                    m_visibilities[m_faces[3*i + 1]] = true;
                }
                if( std::abs(depthc - distance) < 1e-6)
                {
                    m_visibilities[m_faces[3*i + 2]] = true;
                }
            }
        }
    }

    void ZBuffer::ComputePointCloudDepthBuffer( const int32_t blobsize  )
    {
         BVector C = m_P.ComputeLocation();
         for(uint32_t i=0; i<m_Nv; i++)
         {
             BVector a(4,0);
             a[0] = m_vertices[3*i];
             a[1] = m_vertices[3*i + 1];
             a[2] = m_vertices[3*i + 2];
             a[3] = 1;

             BVector pa = m_P*a;

             // The point is within the images
             if(pa[0]>=0 && pa[1] >=0 && pa[0]<m_width && pa[1]<m_height)
             {
                 //This distance calculation could be interpolated
                 //double distance = GetDistanceFromPoint(C, a);
                 double distance =m_P.at(2,0)*a[0] + m_P.at(2,1)*a[1] + m_P.at(2,2)*a[2] + m_P.at(2,3);

                 int32_t j = static_cast<int32_t>(floorf(pa[0]));
                 int32_t k = static_cast<int32_t>(floorf(pa[1]));

                 for(int ji=j-blobsize; ji <= j+blobsize; ji++)
                 for(int ki=k-blobsize; ki <= k+blobsize; ki++)
                 if( ji>=0 && ki >=0 && ji<m_width && ki<m_height )
                 if( (ji-j)*(ji-j) + (ki-k)*(ki-k) < blobsize*blobsize )
                 if( m_db[m_height*ji+ki]<0 || distance < m_db[m_height*ji+ki])
                 {
                    m_db[m_height*ji+ki] = distance;
                    m_ib[m_height*ji+ki] = i;
                 }
             }
         }
    }

     void ZBuffer::ComputePointCloudVisibilities( const int32_t blobsize )
     {
         BVector C = m_P.ComputeLocation();
         for(uint32_t i=0; i<m_Nv; i++)
         {
             BVector a(4,0);
             a[0] = m_vertices[3*i];
             a[1] = m_vertices[3*i + 1];
             a[2] = m_vertices[3*i + 2];
             a[3] = 1;

             BVector pa = m_P*a;

             // The point is within the images
             if(pa[0]>=0 && pa[1] >=0 && pa[0]<m_width && pa[1]<m_height)
             {
                 //This distance calculation could be interpolated
                 double distance = GetDistanceFromPoint(C, a);

                 int32_t j = static_cast<int32_t>(floorf(pa[0]));
                 int32_t k = static_cast<int32_t>(floorf(pa[1]));

                 for(int ji=j-blobsize; ji <= j+blobsize; ji++)
                 for(int ki=k-blobsize; ki <= k+blobsize; ki++)
                 if( ji>=0 && ki >=0 && ji<m_width && ki<m_height )
                 if( (ji-j)*(ji-j) + (ki-k)*(ki-k) < blobsize*blobsize )
                 if( std::abs(distance - m_db[m_height*ji+ki]) < 1e-16)
                 {
                    m_visibilities[i] = true;
                    goto nextvertex;
                 }
             }

             nextvertex:;
         }
     }

    void ZBuffer::ComputeDepthBuffer()
    {
        BVector C = m_P.ComputeLocation();

//#pragma omp parallel for
        for(uint32_t i=0; i<m_Nf; i++)
        {
            BVector a(4,0), b(4,0), c(4,0);
            a[3]=1; b[3]=1; c[3]=1;
            BVector center(2,0);

            SetVertices(m_vertices,m_faces,a,b,c,i);
            BVector pa = m_P*a;
            BVector pb = m_P*b;
            BVector pc = m_P*c;

            double maxx = std::max( {pa[0], pb[0], pc[0]});
            double maxy = std::max( {pa[1], pb[1], pc[1]} );
            double minx = std::min( {pa[0], pb[0], pc[0]} );
            double miny = std::min( {pa[1], pb[1], pc[1]} );

            // The whole triangle is inside the depth image
            if(minx>=0 && miny >=0 && maxx<m_width && maxy<m_height)
            {
                //This distance calculation could be interpolated
                double distance = GetDistanceFromTriangle(C, a,b,c);
                m_faceDepth[i] = distance;


                int32_t startJ = static_cast<int32_t>(floorf(minx));
                int32_t endJ = static_cast<int32_t>(ceilf(maxx));
                int32_t startK = static_cast<int32_t>(floorf(miny));
                int32_t endK = static_cast<int32_t> (ceilf(maxy));

                for(int j=startJ; j <= endJ; j++)
                for(int k=startK; k <= endK; k++)
                if( m_db[m_height*j+k]<0 || distance < m_db[m_height*j+k])
                {
                    center[0]=static_cast<double>(j);
                    center[1]=static_cast<double>(k);

                    double dotproduct[3] = {0.0, 0.0, 0.0};
                    double norms[3] = {0.0, 0.0, 0.0};
                    for(int t=0; t<2;t++)
                    {
                        dotproduct[0]+=(pa[t]-center[t])*(pb[t]-center[t]);
                        dotproduct[1]+=(pb[t]-center[t])*(pc[t]-center[t]);
                        dotproduct[2]+=(pc[t]-center[t])*(pa[t]-center[t]);

                        norms[0]+=(pa[t]-center[t])*(pa[t]-center[t]);
                        norms[1]+=(pb[t]-center[t])*(pb[t]-center[t]);
                        norms[2]+=(pc[t]-center[t])*(pc[t]-center[t]);
                    }
                    for(int t=0; t<3;t++)
                    {
                        norms[t] = std::sqrt(norms[t]);
                    }


                    if(std::abs( std::acos(dotproduct[0]/(norms[0]*norms[1]))  +
                                 std::acos(dotproduct[1]/(norms[1]*norms[2])) +
                                 std::acos(dotproduct[2]/(norms[2]*norms[0])) - 2*M_PI ) < 1e-1 ||
                            (std::abs(pa[0]-center[0]-0.5) <= 0.5 && std::abs(pa[1]-center[1]-0.5) <= 0.5) ||
                            (std::abs(pb[0]-center[0]-0.5) <= 0.5 && std::abs(pb[1]-center[1]-0.5) <= 0.5) ||
                            (std::abs(pc[0]-center[0]-0.5) <= 0.5 && std::abs(pc[1]-center[1]-0.5) <= 0.5) )
                    {
                        m_db[m_height*j+k] = distance;
                        m_ib[m_height*j+k] = i;
                    }
                }
            }
        }
    }

    const std::vector<double>& ZBuffer::ReturnVertexDepth()
    {
        BVector C = m_P.ComputeLocation();

#pragma omp parallel for
        for(uint32_t i=0;i<m_Nv;i++)
        {
            m_vertexDepth[i] = std::sqrt(
                          (m_vertices[3*i] - C[0])*(m_vertices[3*i] - C[0]) +
                          (m_vertices[3*i+1] - C[1])*(m_vertices[3*i+1] - C[1]) +
                          (m_vertices[3*i+2] - C[2])*(m_vertices[3*i+2] - C[2]));
        }
        return m_vertexDepth;
    }

    const std::vector<double>&  ZBuffer::ReturnDepthBuffer() const
    {
        return m_db;
    }

    const std::vector<bool>& ZBuffer::ReturnVertexVisibilities() const
    {
        return m_visibilities;
    }

    std::vector< std::list<uint32_t> > FindCameraVisibilities( const ZBuffer::Mesh mesh, uint32_t height, uint32_t width, const std::vector< ProjectionMatrix >& Ps )
    {
        std::vector< std::vector<bool> > Visibilities(Ps.size(), std::vector<bool>( mesh.Nv ));

        for( uint32_t k=0; k<Ps.size(); ++k )
        {
            ZBuffer zb(mesh,  height, width, Ps[k] );
            zb.ComputeDepthBuffer();
            zb.ComputeVisibilities();
            Visibilities[k] = zb.ReturnVertexVisibilities();
        }

        std::vector< std::list<uint32_t> > VisIndexes(mesh.Nv);
        for(uint32_t k=0; k<Ps.size(); ++k)
        for(uint32_t l=0; l<mesh.Nv; ++l)
        if (Visibilities[k][l])
        {
            VisIndexes[l].push_back(k);
        }

        return VisIndexes;
    }
}
