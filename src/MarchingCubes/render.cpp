#include "render.h"

#include <cmath>
#include <iostream>

Render::Render(const MCParams mcp)
    : m_implicitFuncs(0),
      m_mt(Add),
      m_mc(mcp.isoval, mcp.N, mcp.max, mcp.min), m_currentType(0), m_spin(false),
      m_t(0.0),
      m_pitch(0.0),
      m_yaw(0.0)
{
    m_mc.SetRenderToScreen(true);
}

Render::~Render()
{

}

float Render::GetPitch() const
{
    return m_pitch;
}

float Render::GetYaw() const
{
    return m_yaw;
}

uint32_t Render::GetNumberOfFunctions() const
{
    return m_implicitFuncs.size();
}
void Render::IncrementTime(const double dt)
{
    if(m_spin)
    {
            m_pitch += 4.0;
            //m_yaw   += 2.5;
    }
    m_t += dt;
}

void Render::IncrementMergeType()
{
    m_currentType = (m_currentType + 1) % nvals;
    SwitchFunc(static_cast<MergeType>(m_currentType));
}

void Render::SwitchSpin()
{
    m_spin = !m_spin;
}

void Render::March()
{
    m_mc.March([&](const double* v)
    {
        return m_renderFunction(v,m_t);
    });
}

void print_line(const std::string& line)
{
    std::cout << line << std::endl;
}

void Render::SwitchFunc(MergeType mt)
{
    switch(mt)
    {
    case Add:
    {
        print_line("Merge type: Add");
        print_line("  f = f1 + f2 ");
        print_line("  This blend incorrectly estimates the radius of the ");
        print_line("  blob");
        std::cout << std::endl;

        m_renderFunction = [&](const double* val, const double t)
        {
            double out=0.0;
            for(const auto& f : m_implicitFuncs)
            {
                out += f(val,t);
            }
            return out;
        };
    } break;
    case Bernhardt:
    {
        print_line("Merge type: Clean union (Bernhardt)");
        print_line("  f = f1 + f2 - norm([f1;f2])");
        print_line("  This blend fixes the radius problem in the  ");
        print_line("  simple add blend");
        std::cout << std::endl;
        m_renderFunction = [&](const double* val, const double t)
        {
            double a1 = m_implicitFuncs.begin()->operator()(val,t);
            double a2 = (++m_implicitFuncs.begin())->operator()(val,t);
            return a1+ a2 +std::sqrt(a1*a1 + a2*a2);
            /*
            double out=0.0;
            double outsq=0.0;
            for(const auto& f : m_implicitFuncs)
            {
                out += f(val,t);
                outsq += f(val,t)*f(val,t);
            }
            return out + std::sqrt(outsq);*/
        };
    } break;
    case Barthe:
    {
        print_line("Merge type: Displacement (Pasco)");
        print_line("  f = f1 + f2 - norm([f1;f2]) + a/(1 + (f1/b)^2 + (f2/c)^2)");
        print_line("  At the cost of a slightly incorrect radius the extra");
        print_line("  term smooths the blend.");
        print_line("  There is recent research in to types of possible displacement");
        print_line("  function, allowing control over gradient and blend region.");
        std::cout << std::endl;
        const double a0 = 0.2;
        const double a = 0.1;
        m_renderFunction = [&](const double* val, const double t)
        {
            double a1 = m_implicitFuncs.begin()->operator()(val,t);
            double a2 = (++m_implicitFuncs.begin())->operator()(val,t);
            return a1 + a2 + sqrt(a1*a1 +a2*a2) + a0/((a1/a)*(a1/a) + (a2/a)*(a2/a));
            /*
            double out=0.0;
            double outsq=0.0;
            double asub=1.0;

            for(const auto& f : m_implicitFuncs)
            {
                out += f(val,t);
                outsq += f(val,t)*f(val,t);
                asub += (f(val,t) / a)*(f(val,t) / a);
            }
            return out + std::sqrt(outsq) + (a0 / asub);
            */
        };
    } break;
    case Rockwood:
    {
        print_line("Merge type: Superelliptic (Rockwood)");
        print_line("  f = 1 - (1 - f1)^t/a - (1-f2)^t/b");
        print_line("  This is early research into blending. ");
        print_line("  It can cause discontinuities on the blend.");
        std::cout << std::endl;
        const double tstar=1.3;
        const double r=2;
        m_renderFunction = [&](const double* val, const double t)
        {
            double a1 = m_implicitFuncs.begin()->operator()(val,t);
            double a2 = (++m_implicitFuncs.begin())->operator()(val,t);

            return  1.0 -
                    std::pow(std::max((1-a1)/r,0.0),tstar) -
                    std::pow(std::max((1-a2)/r,0.0),tstar);
        };
    } break;
    case Functional:
    {
        print_line("Merge type: Functional (Barthe)");
        print_line("  f = min(f1,f2) - H(|f1-f2|)");
        print_line("  H is a piecewise cubic defining the shape of the ");
        print_line("  blend. It allows the user to define how the final ");
        print_line("  blend looks, but can give unusual results. ");
        std::cout << std::endl;

        m_renderFunction = [&](const double* val, const double t)
        {
            double a1 = m_implicitFuncs.begin()->operator()(val,t);
            double a2 = (++m_implicitFuncs.begin())->operator()(val,t);
            return std::min(a1,a2) + 1*std::pow(std::abs(a1-a2),3) + 0.8*std::abs(a1-a2);
        };
    } break;
    case Union:
    {
        print_line("Merge type: Union (Barthe)");
        print_line("  f = max(f1,f2)");
        print_line("  A simple approach which can cause discontinuities on");
        print_line("  some surfaces.");
        std::cout << std::endl;
        m_renderFunction = [&](const double* val, const double t)
        {
            double a1 = m_implicitFuncs.begin()->operator()(val,t);
            double a2 = (++m_implicitFuncs.begin())->operator()(val,t);
            return std::max(a1,a2);
            /*
            double max = std::numeric_limits<double>::min();

            double v;
            for(const auto& f : m_implicitFuncs )
            {
                v = f(val,t);
                max = (v > max) ? v : max;
            }
            return max;
            */
        };
    } break;
    default:
        std::cerr << "No known merge type" <<std::endl;
        break;
    }
}

void Render::AddFunction( const ImplicitFunction ifunc )
{
    m_implicitFuncs.push_back(ifunc);
    SwitchFunc(static_cast<MergeType>(m_currentType));
}

void GLR::Resize( GLsizei iWidth, GLsizei iHeight )
{
        GLfloat fAspect, fHalfWorldSize = (1.4142135623730950488016887242097/2);

        glViewport( 0, 0, iWidth, iHeight );
        glMatrixMode (GL_PROJECTION);
        glLoadIdentity ();

        if(iWidth <= iHeight)
        {
                fAspect = (GLfloat)iHeight / (GLfloat)iWidth;
                glOrtho(-fHalfWorldSize, fHalfWorldSize, -fHalfWorldSize*fAspect,
                        fHalfWorldSize*fAspect, -10*fHalfWorldSize, 10*fHalfWorldSize);
        }
        else
        {
                fAspect = (GLfloat)iWidth / (GLfloat)iHeight;
                glOrtho(-fHalfWorldSize*fAspect, fHalfWorldSize*fAspect, -fHalfWorldSize,
                        fHalfWorldSize, -10*fHalfWorldSize, 10*fHalfWorldSize);
        }

        glMatrixMode( GL_MODELVIEW );
}

void GLR::Keyboard(unsigned char cKey, int iX, int iY)
{
    (void) iX;
    (void) iY;
    switch(cKey)
    {
        case 'w' :
        {
            GetRender().SwitchSpin();
        } break;
        case 's' :
        {
            GetRender().IncrementMergeType();
        } break;
        case 27:
        {
        } break;
    }
}


void GLR::Special(int iKey, int iX, int iY)
{
    (void) iKey;
    (void) iX;
    (void) iY;
}

void GLR::Idle()
{
        glutPostRedisplay();
}

void GLR::DrawScene()
{
        glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
        glPushMatrix();

        glTranslatef(0.0, 0.0, 0.0);
        glRotatef( 0.0,    1.0, 0.0, 0.0);
        glRotatef(-GetRender().GetPitch(), 0.0, 1.0, 0.0);
        glRotatef( GetRender().GetYaw(),   0.0, 0.0, 1.0);

        //glPushAttrib(GL_LIGHTING_BIT);
                //glDisable(GL_LIGHTING);
                //glColor3f(1.0, 1.0, 1.0);
                //glutWireCube(1.0);
        //glPopAttrib();

        glPushMatrix();
        glTranslatef(0.0, 0.0, 0.0);
        glBegin(GL_TRIANGLES);

        GetRender().March();
        GetRender().IncrementTime(0.03);

        glEnd();
        glPopMatrix();
        glPopMatrix();
        glutSwapBuffers();
}

GLenum    ePolygonMode = GL_FILL;
//static const GLfloat afAmbientWhite [] = {0.25, 0.25, 0.25, 1.00};
//static const GLfloat afAmbientRed   [] = {0.25, 0.00, 0.00, 1.00};
static const GLfloat afAmbientGreen [] = {0.00, 0.25, 0.00, 1.00};
static const GLfloat afAmbientBlue  [] = {0.00, 0.00, 0.25, 1.00};
//static const GLfloat afDiffuseWhite [] = {0.75, 0.75, 0.75, 1.00};
//static const GLfloat afDiffuseRed   [] = {0.75, 0.00, 0.00, 1.00};
static const GLfloat afDiffuseGreen [] = {0.00, 0.75, 0.00, 1.00};
static const GLfloat afDiffuseBlue  [] = {0.00, 0.00, 0.75, 1.00};
//static const GLfloat afSpecularWhite[] = {1.00, 1.00, 1.00, 1.00};
//static const GLfloat afSpecularRed  [] = {1.00, 0.25, 0.25, 1.00};
//static const GLfloat afSpecularGreen[] = {0.25, 1.00, 0.25, 1.00};
//static const GLfloat afSpecularBlue [] = {0.25, 0.25, 1.00, 1.00};

void GLR::Create(const Render::MCParams mcp)
{
    m_mcp=mcp;

    GLfloat afPropertiesAmbient [] = {0.50, 0.50, 0.50, 1.00};
    GLfloat afPropertiesDiffuse [] = {0.75, 0.75, 0.75, 1.00};
    GLfloat afPropertiesSpecular[] = {1.00, 1.00, 1.00, 1.00};

    GLsizei iWidth = 640.0;
    GLsizei iHeight = 480.0;

    int argc=0;
    char **argv=NULL;

    glutInit(&argc, argv);
    glutInitWindowPosition( 0, 0);
    glutInitWindowSize(iWidth, iHeight);
    glutInitDisplayMode( GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE );
    glutCreateWindow( "Marching Cubes" );
    glutDisplayFunc( DrawScene );
    glutIdleFunc( Idle );
    glutReshapeFunc( Resize );
    glutKeyboardFunc( Keyboard );
    glutSpecialFunc( Special );

    glClearColor( 0.0, 0.0, 0.0, 1.0 );
    glClearDepth( 1.0 );

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
    glPolygonMode(GL_FRONT_AND_BACK, ePolygonMode);

    glLightfv( GL_LIGHT0, GL_AMBIENT,  afPropertiesAmbient);
    glLightfv( GL_LIGHT0, GL_DIFFUSE,  afPropertiesDiffuse);
    glLightfv( GL_LIGHT0, GL_SPECULAR, afPropertiesSpecular);
    glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, 1.0);

    glEnable( GL_LIGHT0 );

    glMaterialfv(GL_FRONT,  GL_AMBIENT,   afAmbientGreen);
    glMaterialfv(GL_FRONT,  GL_DIFFUSE,   afDiffuseGreen);
    glMaterialfv(GL_BACK, GL_AMBIENT,   afAmbientBlue);
    glMaterialfv(GL_BACK, GL_DIFFUSE,   afDiffuseBlue);
    //glMaterialfv(GL_FRONT, GL_SPECULAR,  afSpecularWhite);
    glMaterialf( GL_FRONT, GL_SHININESS, 1.0);

    Resize(iWidth, iHeight);

    std::cout << "Marching Cubes v1.0 Daniel Beale " << std::endl;
}

void GLR::Start()
{
    if( GetRender().GetNumberOfFunctions()  > 0)
    {
        glutMainLoop();
    }
    else
    {
        std::cerr << "Render: Not enough implicit functions" << std::endl;
    }
}

void GLR::AddFunction(const Render::ImplicitFunction ifunc )
{
    GetRender().AddFunction(ifunc);
}
