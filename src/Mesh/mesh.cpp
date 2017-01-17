
#include "mesh.h"

Mesh::Mesh( const std::string& filename ) :
    m_NVerts(0), m_NFaces(0)
{
    Load(filename);
}

inline void GeneratePoint(const double& theta, const double& phi, double& a, double& b, double& c)
{
    a = std::sin(phi)*std::cos(theta);
    b = std::sin(phi)*std::sin(theta),
    c = std::cos(phi);
}

inline void Normalise( double& an, double& bn, double& cn, const double a, const double b, const double c  )
{
    double sum = a+b+c;
    an = a / sum;
    bn = b / sum;
    cn = c / sum;
}

Mesh::Mesh() :
    m_vertices( 0 ), m_faces(0), m_normals(0), m_colours(0),
    m_NVerts(0), m_NFaces(0), m_NRing(0)
{

}

Mesh::Mesh( const Primitive primitive )
{
    switch(primitive)
    {
        case UnitSphere :
        {
            const uint32_t N = 30;

            std::vector<double> phi(N+1);
            std::vector<double> theta(N+1);
            for(uint32_t i=0;i<=N;i++)
            {
              phi[i] = ((M_PI)/ static_cast<double>(N))*static_cast<double>(i);
              theta[i] = ((2*M_PI) / static_cast<double>(N))*static_cast<double>(i);
            }

            m_NVerts = N*N*4;
            m_NFaces =  N*N*2;
            m_vertices.resize( 3*m_NVerts, 0.0 );
            m_normals.resize( 3*m_NVerts, 0.0 );
            m_faces.resize( 3*m_NFaces, 0 );
            m_NRing.resize(m_NVerts);

            auto vit = m_vertices.begin();
            auto fit = m_faces.begin();
            auto nit = m_normals.begin();
            uint32_t oi =0;

            for(uint32_t j=0; j<N; j++)
            for(uint32_t i=0; i<N; i++)
            {
                GeneratePoint(theta[i], phi[j], *(vit), *(vit+1), *(vit+2));
                Normalise(*(nit), *(nit+1), *(nit+2), *(vit), *(vit+1), *(vit+2));
                GeneratePoint(theta[i+1], phi[j], *(vit+3), *(vit+4), *(vit+5));
                Normalise(*(nit+3), *(nit+4), *(nit+5), *(vit+3), *(vit+4), *(vit+5));
                GeneratePoint(theta[i], phi[j+1], *(vit+6), *(vit+7), *(vit+8));
                Normalise(*(nit+6), *(nit+7), *(nit+8), *(vit+6), *(vit+7), *(vit+8));
                GeneratePoint(theta[i+1], phi[j+1], *(vit+9), *(vit+10), *(vit+11));
                Normalise(*(nit+9), *(nit+10), *(nit+11), *(vit+9), *(vit+10), *(vit+11));
                vit+=12;
                nit+=12;

                *(fit++) = oi;
                *(fit++) = oi+1;
                *(fit++) = oi+2;
                *(fit++) = oi+2;
                *(fit++) = oi+1;
                *(fit++) = oi+3;
                oi+=4;
            }
        } break;
    }
}

Mesh::Mesh(const std::vector< std::vector<double> >& vertices) :
    m_vertices( vertices.size()*3 ), m_faces(0), m_normals(0), m_colours(0),
    m_NVerts(vertices.size()), m_NFaces(0), m_NRing(m_NVerts)
{
    uint32_t index =0;
    for(const auto& v : vertices )
    {
        m_vertices[3*index] = v[0];
        m_vertices[3*index + 1] = v[1];
        m_vertices[3*index + 2] = v[2];
        ++index;
    }
}

Mesh::Mesh(const std::vector<PointTypes::Pointd3> &vertices, const std::vector<PointTypes::Pointi3> &faces) :
    m_vertices( vertices.size()*3 ), m_faces(faces.size()*3), m_normals(0), m_colours(0),
    m_NVerts(vertices.size()), m_NFaces(faces.size()), m_NRing(m_NVerts)
{
    uint32_t index =0;
    for(const auto& v : vertices )
    {
        m_vertices[3*index] = v.p[0];
        m_vertices[3*index + 1] = v.p[1];
        m_vertices[3*index + 2] = v.p[2];
        ++index;
    }

    index=0;
    for(const auto& f : faces)
    {
        m_faces[3*index] = f.p[0];
        m_faces[3*index + 1] = f.p[1];
        m_faces[3*index + 2] = f.p[2];
        ++index;
    }
}

Mesh::Mesh( std::vector< double > vertices, std::vector< uint32_t > faces,
      std::vector< double > normals, std::vector< double > colours ) :
    m_vertices(vertices), m_faces(faces), m_normals(normals), m_colours(colours),
    m_NVerts(m_vertices.size()), m_NFaces(m_faces.size()), m_NRing(m_NVerts)
{

}

std::vector<double> vertices;
std::vector<double> colours;
std::vector<double> normals;
std::vector<uint32_t> faces;

static int vertex_cb( p_ply_argument argument )
{
    long length, value_index;
    ply_get_argument_property(argument, NULL, &length, &value_index);
    long instance_index,idata;
    p_ply_element element;
    void *pdata;

    ply_get_argument_element(argument, &element, &instance_index);
    ply_get_argument_user_data(argument, &pdata, &idata);

    if(idata < 3)
    {
        vertices[instance_index*3 + idata] = ply_get_argument_value(argument);
    }
    else if (idata < 6)
    {
        colours[instance_index*3 + idata - 3] = ply_get_argument_value(argument) / 255.0;
    }
    else if(idata < 9)
    {
        normals[instance_index*3 + idata - 6] = ply_get_argument_value(argument);
    }
    return 1;
}

static int face_cb( p_ply_argument argument )
{
    long length, value_index, instance_index;
    p_ply_element element;
    ply_get_argument_property(argument, NULL, &length, &value_index);
    ply_get_argument_element(argument, &element, &instance_index);
    if(value_index >=0)
        faces[instance_index*3 + value_index] = ply_get_argument_value(argument);
    return 1;
}

int read_vertices(const std::string& file)
{
    long nvertices, ncolours, nfaces, nnormals;
    p_ply ply = ply_open(file.c_str(), NULL, 0, NULL);

    if (!ply_read_header(ply))
    {
        std::cerr << "PointCloudProjection: Cannot read ply header " << std::endl;
        return 1;
    }

    nvertices = ply_set_read_cb(ply, "vertex", "x", vertex_cb, NULL, 0);
    nvertices = ply_set_read_cb(ply, "vertex", "y", vertex_cb, NULL, 1);
    nvertices = ply_set_read_cb(ply, "vertex", "z", vertex_cb, NULL, 2);
    ncolours = ply_set_read_cb(ply, "vertex", "red", vertex_cb, NULL, 3);
    ncolours = ply_set_read_cb(ply, "vertex", "green", vertex_cb, NULL, 4);
    ncolours = ply_set_read_cb(ply, "vertex", "blue", vertex_cb, NULL, 5);
    nnormals = ply_set_read_cb(ply, "vertex", "nx", vertex_cb, NULL, 6);
    nnormals = ply_set_read_cb(ply, "vertex", "ny", vertex_cb, NULL, 7);
    nnormals = ply_set_read_cb(ply, "vertex", "nz", vertex_cb, NULL, 8);

    nfaces = ply_set_read_cb(ply, "face", "vertex_indices", face_cb, NULL, 0);

    vertices = std::vector<double>(nvertices*3, 0.0);
    colours = std::vector<double>(ncolours*3, 0.0);
    normals = std::vector<double>(nnormals*3, 0.0);
    faces = std::vector<uint32_t>(nfaces*3, 0.0);

    if (!ply_read(ply))
    {
        std::cerr << "Mesh::Load::read_vertices :: Unable to read points " << std::endl;
        return 1;
    }

    ply_close(ply);

    return 0;
}

void Mesh::Load( const std::string& filename )
{
    read_vertices(filename);
    m_vertices = vertices;
    m_colours = colours;
    m_normals = normals;
    m_faces = faces;
    m_NVerts = vertices.size()/3;
    m_NFaces = faces.size()/3;
    m_NRing.resize(m_NVerts);
}

void Mesh::Write( const std::string& filename )
{
    std::ofstream f(filename, std::ofstream::out);
    f << "ply\n";
    f << "format ascii 1.0\n";
    f << "comment ELib::Mesh generated\n";
    f << "element vertex " << m_NVerts <<"\n";
    f << "property float x" << "\n";
    f << "property float y" << "\n";
    f << "property float z" << "\n";

    if( m_normals.size()/3 == m_NVerts)
    {
        f << "property float nx\n";
        f << "property float ny\n";
        f << "property float nz\n";
    }

    if( m_colours.size()/3 == m_NVerts )
    {
        f << "property uchar red\n";
        f << "property uchar green\n";
        f << "property uchar blue\n";
    }

    f << "element face " << m_NFaces << "\n";
    f << "property list uchar int vertex_indices\n";

    f << "end_header\n";

    const bool donormals =  m_normals.size()/3 == m_NVerts;
    const bool docolours = m_colours.size()/3 == m_NVerts;

    for(uint32_t i=0; i<m_NVerts; ++i)
    {
        f << m_vertices[3*i] << " " << m_vertices[3*i + 1] << " " << m_vertices[3*i + 2] << " ";

        if(donormals)
        {
            f << m_normals[3*i] << " " << m_normals[3*i + 1] << " " << m_normals[3*i + 2] << " ";
        }

        if(docolours)
        {
            f << static_cast<uint32_t>(m_colours[3*i]*255) << " " <<
                 static_cast<uint32_t>(m_colours[3*i + 1]*255)<< " " <<
                 static_cast<uint32_t>(m_colours[3*i + 2]*255) << " ";
        }

        f << "\n";
    }

    if( m_NFaces > 0 )
    {
        for(uint32_t i=0; i<m_NFaces; ++i)
        {
            f << 3 << " " << m_faces[3*i] << " " << m_faces[3*i + 1] << " " << m_faces[3*i + 2] << "\n";
        }
    }
    f.close();
}

void Mesh::FillNRing()
{
    for(uint32_t i=0; i<m_NFaces; ++i)
    {
        m_NRing[m_faces[i*3]].insert(m_faces[i*3+1]);
        m_NRing[m_faces[i*3]].insert(m_faces[i*3+2]);
        m_NRing[m_faces[i*3+1]].insert(m_faces[i*3]);
        m_NRing[m_faces[i*3+1]].insert(m_faces[i*3+2]);
        m_NRing[m_faces[i*3+2]].insert(m_faces[i*3]);
        m_NRing[m_faces[i*3+2]].insert(m_faces[i*3+1]);
    }
}

std::set<uint32_t> Mesh::GetNRing( const uint32_t i)
{
    return m_NRing[i];
}

std::vector<double> Mesh::GetVertex( const uint32_t i) const
{
    std::vector<double> vertex(3);
    assert(m_vertices.begin() + i + 3 <= m_vertices.end());
    std::copy(m_vertices.begin() + i, m_vertices.begin() + i + 3, vertex.begin());
    return vertex;
}

std::vector<uint32_t> Mesh::GetFace( const uint32_t i) const
{
    std::vector<uint32_t> face(3);
    assert(m_faces.begin() + i + 3 <= m_faces.end());
    std::copy(m_faces.begin() + i, m_faces.begin() + i + 3, face.begin());
    return face;
}

std::vector<double> Mesh::GetColour( const uint32_t i) const
{
    std::vector<double> colour(3);
    assert(m_colours.begin() + i + 3 <= m_colours.end());
    std::copy(m_colours.begin() + i, m_colours.begin() + i + 3, colour.begin());
    return colour;
}

std::vector<double> Mesh::GetNormal( const uint32_t i) const
{
    std::vector<double> normal(3);
    assert(m_normals.begin() + i + 3 <= m_normals.end());
    std::copy(m_normals.begin() + i, m_normals.begin() + i + 3, normal.begin());
    return normal;
}

void Mesh::ApplyNormals( const Mesh& mesh, const NearestNeighbours& nn)
{
    std::vector<KD::KeyPair> kp;
    m_normals.resize( m_NVerts*3 );

//#pragma omp parallel for
    for(uint32_t i=0; i<m_vertices.size()/3; ++i)
    {
        kp = nn.Search( m_vertices.data()+i*3, 1 );

        m_normals[i*3] = mesh.GetNormals()[kp[0].index*3];
        m_normals[i*3+1] = mesh.GetNormals()[kp[0].index*3+1];
        m_normals[i*3+2] = mesh.GetNormals()[kp[0].index*3+2];
    }
}

void Mesh::ComputeTangents( double* n1, double* n2, const uint32_t face)
{
    double d1=0, d2=0;
    for(uint8_t i=0; i<3;++i)
    {
        n1[i] = m_vertices[m_faces[face*3]*3 + i] - m_vertices[m_faces[face*3 + 1]*3 + i];
        n2[i] = m_vertices[m_faces[face*3]*3 + i] - m_vertices[m_faces[face*3 + 2]*3 + i];
        d1 += n1[i]*n1[i];
        d2 += n2[i]*n2[i];
    }
    d1 = std::sqrt(d1);
    d2 = std::sqrt(d2);
    for(uint8_t i=0; i<3;++i)
    {
        n1[i]/=d1;
        n2[i]/=d2;
    }
}

void Mesh::TranslateVertices( const std::vector<double>& A, const std::vector<double>& b )
{
    if(A.size() != 9)
    {
        std::cerr << "Error:: Mesh::TranslateVertices incorrect size of A" << std::endl;
        assert(false);
    }
    if(b.size() != 3)
    {
        std::cerr << "Error:: Mesh::TranslateVertices incorrect size of b" << std::endl;
        assert(false);
    }

    for(uint32_t i=0;i<m_NVerts; ++i)
    {
        double x = A[0]*m_vertices[3*i] +A[1]*m_vertices[3*i+1] + A[2]*m_vertices[3*i+2] + b[0];
        double y = A[3]*m_vertices[3*i] +A[4]*m_vertices[3*i+1] + A[5]*m_vertices[3*i+2] + b[1];
        double z = A[6]*m_vertices[3*i] +A[7]*m_vertices[3*i+1] + A[8]*m_vertices[3*i+2] + b[2];
        m_vertices[3*i] = x;
        m_vertices[3*i+1] = y;
        m_vertices[3*i+2] = z;
    }
}

inline void CrossProduct(const double n1[3], const double n2[3], double* cp)
{
    cp[0] = n1[1]*n2[2] - n2[1]*n1[2];
    cp[1] = n2[0]*n1[2] - n1[0]*n2[2];
    cp[2] = n1[0]*n2[1] - n2[0]*n1[1];
}

inline double distance(const double* a, const double* b)
{
    return std::sqrt((a[0]-b[0])*(a[0]-b[0]) +
                     (a[1]-b[1])*(a[1]-b[1]) +
                     (a[2]-b[2])*(a[2]-b[2]));
}

double Mesh::ComputeAverageEdgeLength() const
{
    double edgeLength = 0.0;
    const double* dat = m_vertices.data();
    for(uint32_t i=0; i<m_NFaces; ++i)
    {
        edgeLength += distance( &dat[3*m_faces[3*i]], &dat[3*m_faces[3*i+1]] ) / 3;
        edgeLength += distance( &dat[3*m_faces[3*i+1]], &dat[3*m_faces[3*i+2]] ) / 3;
        edgeLength += distance( &dat[3*m_faces[3*i+2]], &dat[3*m_faces[3*i]] ) / 3;
    }
    edgeLength /= m_NFaces;
    return edgeLength;
}

std::vector<double> Mesh::ComputeMean() const
{
    std::vector<double> mean(3, 0.0);
    for(uint32_t i=0; i<m_NVerts; ++i)
    {
        mean[0] += m_vertices[i*3];
        mean[1] += m_vertices[i*3 + 1];
        mean[2] += m_vertices[i*3 + 2];
    }
    mean[0] /= m_NVerts;
    mean[1] /= m_NVerts;
    mean[2] /= m_NVerts;

    return mean;
}

inline double D2(const double* a, const double* b)
{
    return ((a[0]-b[0])*(a[0]-b[0]) +
             (a[1]-b[1])*(a[1]-b[1]) +
             (a[2]-b[2])*(a[2]-b[2]));
}

std::vector<double> Mesh::ComputeClosestPoint( const std::vector<double>& pt ) const
{
    double dist = std::numeric_limits<double>::max();
    std::vector<double> minpt(3,0.0);
    const double* ptp = &pt[0];

    for(uint32_t i=0; i<m_NVerts; ++i)
    {
        double d = D2( &m_vertices[3*i], ptp );
        if(d < dist)
        {
            dist = d;
            minpt[0] += m_vertices[i*3];
            minpt[1] += m_vertices[i*3 + 1];
            minpt[2] += m_vertices[i*3 + 2];
        }
    }
    return minpt;
}

double Mesh::ComputeAverageDistance( const std::vector<double>& pt ) const
{
    double dist = 0.0;
    const double* ptp = &pt[0];

    for(uint32_t i=0; i<m_NVerts; ++i)
    {
        dist += distance( &m_vertices[3*i], ptp );
    }
    dist /= m_NVerts;
    return dist;
}

void Mesh::ComputeNormalsFromPointsAndFaces()
{
    if(m_faces.size() == 0 || m_vertices.size() == 0)
    {
        std::cerr << "Error::Mesh::ComputeNormalsFromPointsAndFaces:: There are no faces of vertices" << std::endl;
        assert(0);
    }

    m_normals = std::vector<double>(3*m_NVerts, 0.0);
    std::vector<uint32_t> perVertexNs( m_NVerts, 0);

#pragma omp parallel for
    for(uint32_t i=0; i<m_NFaces; ++i)
    {
        double n1[3],n2[3],cp[3],nm; // tangents
        ComputeTangents( n1, n2, i);
        CrossProduct( n1, n2, cp);
        nm = std::sqrt(cp[0]*cp[0] + cp[1]*cp[1] + cp[2]*cp[2]);

        if(nm == 0)
        {
#pragma omp critical
            std::cerr <<"Warning::ComputeNormalsFromPointsAndFaces:: There is a zero area face. " << std::endl;
        }

        for(uint8_t j=0; j<3; ++j)
        {
            if(nm > 0) cp[j] /= nm;
            m_normals[3*m_faces[3*i] + j] += cp[j];
            m_normals[3*m_faces[3*i+1] + j] += cp[j];
            m_normals[3*m_faces[3*i+2] + j] += cp[j];

#pragma omp critical
            {
                perVertexNs[m_faces[3*i]]++;
                perVertexNs[m_faces[3*i+1]]++;
                perVertexNs[m_faces[3*i+2]]++;
            }
        }
    }

#pragma omp parallel for
    for(uint32_t i=0; i<m_NVerts; ++i)
    {
        if(perVertexNs[i] > 0)
        {
            m_normals[3*i] /= perVertexNs[i];
            m_normals[3*i+1] /= perVertexNs[i];
            m_normals[3*i+2] /= perVertexNs[i];
            double nm = std::sqrt(m_normals[3*i]*m_normals[3*i] +
                                   m_normals[3*i+1]*m_normals[3*i+1] +
                                   m_normals[3*i+2]*m_normals[3*i+2]);
            if(nm > 0)
            {
                m_normals[3*i] /= nm;
                m_normals[3*i+1] /= nm;
                m_normals[3*i+2] /= nm;
            }
            else
            {
#pragma omp critical
                std::cerr << "Warning::ComputeNormalsFromPointsAndFaces:: There is a zero length normal." << std::endl;
            }
        }
        else
        {
#pragma omp critical
            std::cerr << "Warning::ComputeNormalsFromPointsAndFaces:: There are unreferenced vertices." << std::endl;
        }
    }
}
