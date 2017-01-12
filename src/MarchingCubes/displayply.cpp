/**
* This is a c++ program to display a ply file using 
* open gl
*/

#include <vector>
#include <iostream>
#include <functional>
#include <cmath>
#include "marchingcubes.h";
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>

GLvoid vPrintHelp();

void vIdle();
void vDrawScene(); 
void vResize(GLsizei, GLsizei);
void vKeyboard(unsigned char cKey, int iX, int iY);
void vSpecial(int iKey, int iX, int iY);
void March();

enum Type
{
  Gaussian,
  Merge,
  Blend,
  Poly
};
void switchFunc(Type f);
static uint32_t type=0;

static const GLfloat afAmbientWhite [] = {0.25, 0.25, 0.25, 1.00}; 
static const GLfloat afAmbientRed   [] = {0.25, 0.00, 0.00, 1.00}; 
static const GLfloat afAmbientGreen [] = {0.00, 0.25, 0.00, 1.00}; 
static const GLfloat afAmbientBlue  [] = {0.00, 0.00, 0.25, 1.00}; 
static const GLfloat afDiffuseWhite [] = {0.75, 0.75, 0.75, 1.00}; 
static const GLfloat afDiffuseRed   [] = {0.75, 0.00, 0.00, 1.00}; 
static const GLfloat afDiffuseGreen [] = {0.00, 0.75, 0.00, 1.00}; 
static const GLfloat afDiffuseBlue  [] = {0.00, 0.00, 0.75, 1.00}; 
static const GLfloat afSpecularWhite[] = {1.00, 1.00, 1.00, 1.00}; 
static const GLfloat afSpecularRed  [] = {1.00, 0.25, 0.25, 1.00}; 
static const GLfloat afSpecularGreen[] = {0.25, 1.00, 0.25, 1.00}; 
static const GLfloat afSpecularBlue [] = {0.25, 0.25, 1.00, 1.00}; 

static double xpos = 0.0;
static double xpos2 = 0.0;
static double t=0.0;
static bool bSpin;

GLenum    ePolygonMode = GL_FILL;

const double radius = 0.3;
static double max[3] = {0.85,0.35,0.35};
static double min[3] = {-0.85,-0.35,-0.35};
static std::vector<MC::Pointd3> vertices;
static std::vector<MC::Pointi3> faces;
static std::function<double(const double*)> func;

static MC::MarchingCubes mc(0.0001, 50, max, min);

int main(int argc, char **argv) 
{ 
        mc.SetRenderToScreen( true );
        GLfloat afPropertiesAmbient [] = {0.50, 0.50, 0.50, 1.00}; 
        GLfloat afPropertiesDiffuse [] = {0.75, 0.75, 0.75, 1.00}; 
        GLfloat afPropertiesSpecular[] = {1.00, 1.00, 1.00, 1.00}; 

        GLsizei iWidth = 640.0; 
        GLsizei iHeight = 480.0; 

        glutInit(&argc, argv);
        glutInitWindowPosition( 0, 0);
        glutInitWindowSize(iWidth, iHeight);
        glutInitDisplayMode( GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE );
        glutCreateWindow( "Marching Cubes" );
        glutDisplayFunc( vDrawScene );
        glutIdleFunc( vIdle );
        glutReshapeFunc( vResize );
        glutKeyboardFunc( vKeyboard );
        glutSpecialFunc( vSpecial );

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

        vResize(iWidth, iHeight); 

        vPrintHelp();
        switchFunc((Type)type);
        glutMainLoop(); 
        return 0;
}

GLvoid vPrintHelp()
{
        std::cout << "Marching Cubes" << std::endl;

}


void vResize( GLsizei iWidth, GLsizei iHeight ) 
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

void vKeyboard(unsigned char cKey, int iX, int iY)
{
        switch(cKey)
        {
                case 'w' :
                {
                        bSpin=!bSpin;
                } break;
                case 's' :
                {
                  type=(type+1)%3;
                  switchFunc((Type)type);
                }
       }
}


void vSpecial(int iKey, int iX, int iY)
{
        
}

void vIdle()
{
        glutPostRedisplay();
}

void vDrawScene() 
{ 
        static GLfloat fPitch = 0.0;
        static GLfloat fYaw   = 0.0;
        static GLfloat fTime = 0.0;

        glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT ); 

        glPushMatrix();

        //vSetTime(fTime);
        if(bSpin)
        {
                fPitch += 4.0;
                //fYaw   += 2.5;
        }

        glTranslatef(0.0, 0.0, 0.0);  
        glRotatef( 0.0, 1.0, 0.0, 0.0);
        glRotatef(    -fPitch, 0.0, 1.0, 0.0);
        glRotatef(    fYaw, 0.0, 0.0, 1.0);

        //glPushAttrib(GL_LIGHTING_BIT);
                //glDisable(GL_LIGHTING);
                //glColor3f(1.0, 1.0, 1.0);
                //glutWireCube(1.0); 
        //glPopAttrib(); 


        glPushMatrix(); 
        glTranslatef(0.0, 0.0, 0.0);
        glBegin(GL_TRIANGLES);
        March();
        glEnd();
        glPopMatrix(); 


        glPopMatrix(); 

        glutSwapBuffers(); 
}

void switchFunc(Type f)
{
  switch(f)
  {
    case Gaussian:
      std::cout << "Merge using simple addition." << std::endl;
      func = [=](const double* val)
      {
          return std::exp(-((val[0]-xpos)*(val[0]-xpos) + val[1]*val[1] + val[2]*val[2])/(radius*radius)) + 
                 std::exp(-((val[0]-xpos2)*(val[0]-xpos2) + val[1]*val[1] + val[2]*val[2])/(radius*radius)) - 
                 2*std::exp(-1);
      };
      break;
    case Merge:
      std::cout << "f1 +f2 + norm([f1;f2]) - Bernhardt" << std::endl;
      func = [=](const double* val)
      {
          double f1 = std::exp(-((val[0]-xpos)*(val[0]-xpos) + val[1]*val[1] + val[2]*val[2])/(radius*radius)) - std::exp(-1);
          double f2 = std::exp(-((val[0]-xpos2)*(val[0]-xpos2) + val[1]*val[1] + val[2]*val[2])/(radius*radius)) - std::exp(-1);
          return f1 + f2 + std::sqrt(f1*f1 + f2*f2);
      };
      break;
    case Blend:
      std::cout << "f1 +f2 + norm([f1;f2]) + G(f1,f2) - Barthe" << std::endl;
      const double a0 = 0.2;
      const double a1 = 0.1;
      const double a2 = 0.1;
      func = [=](const double* val)
      {
          double f1 = std::exp(-((val[0]-xpos)*(val[0]-xpos) + val[1]*val[1] + val[2]*val[2])/(radius*radius)) - std::exp(-1);
          double f2 = std::exp(-((val[0]-xpos2)*(val[0]-xpos2) + val[1]*val[1] + val[2]*val[2])/(radius*radius)) - std::exp(-1);
          return f1 + f2 + std::sqrt(f1*f1 + f2*f2) + a0/(1+(f1/a1)*(f1/a1)+(f2/a2)*(f2/a2));
      };
      break;
  }   
}

inline void cross(const double* a, const double* b, double cr[3] )
{
  cr[0] = a[1]*b[2] - a[2]*b[1];
  cr[1] = -(a[0]*b[2] - a[2]*b[0]);
  cr[2] = a[0]*b[1] - a[1]*b[0];
}

void March()
{
    double cr[3], a[3], b[3];

    mc.March(func);
    xpos = 0.5*std::sin(t);
    xpos2 = -0.5*std::sin(t);
    t+=0.01;
#ifndef RENDER_USING_GLUT
    vertices = mc.GetPoints();
    faces = mc.GetFaces();
    uint32_t k,l;
    double nm = 0.0;
    for( const auto& f : faces )
    {
      for( k=0; k<3; k++)
      { 
        for( l=0;l<3;l++)
        {
          a[l] = vertices.at(f.p[(k+1)%3]).p[l]-vertices.at(f.p[k]).p[l];
          b[l] = vertices.at(f.p[(k+2)%3]).p[l]-vertices.at(f.p[k]).p[l];
        }
        cross(a,b,cr);
        nm = std::sqrt(cr[0]*cr[0] + cr[1]*cr[1] + cr[2]*cr[2]);
        glNormal3f(cr[0]/nm,cr[1]/nm,cr[2]/nm);
        
        glVertex3f(vertices.at(f.p[k]).p[0], 
                   vertices.at(f.p[k]).p[1], 
                   vertices.at(f.p[k]).p[2]);
      }
      
      
    }
#endif
}
