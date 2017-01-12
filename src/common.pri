
QMAKE_MACOSX_DEPLOYMENT_TARGET = 10.14
QMAKE_CXX = clang-omp++
QMAKE_CC = clang-omp
QMAKE_CXXFLAGS_RELEASE -= -O
QMAKE_CXXFLAGS_RELEASE -= -O1
QMAKE_CXXFLAGS_RELEASE -= -O2
QMAKE_CXXFLAGS_RELEASE -= -O3
QMAKE_CXXFLAGS_RELEASE *= -Ofast

QMAKE_CXXFLAGS += -std=c++11
#QMAKE_CXXFLAGS += -fopenmp
#LIBS += -fopenmp

DESTDIR     = $$OUT_DIR/../bin
BUILD_DIR   = $$OUT_DIR/../build/$$MODULE_NAME

TARGET      = $$MODULENAME

CONFIG   += c++11
CONFIG += std

INCLUDEPATH +=/usr/local/include \
              /usr/local/include/libiomp/ \
              $$DIR/../Aux/ \
              $$DIR/../GUI/ \
              $$DIR/../Statistics/ \
              $$DIR/../Image/ \
              $$DIR/../Matrix/ \
              $$DIR/../Mesh/ \
              $$DIR/../MaxFlow/ \
              $$DIR/../rply/ \
              $$DIR/../MarchingCubes/ \
              $$DIR/../Carving/ \
              $$DIR/../NearestNeighbours \
              $$DIR/../KDTree \
              $$DIR/../Eigen
              
DEPENDPATH += /usr/local/include \
              /usr/local/include/libiomp/ \
              $$DIR/../Aux/ \
              $$DIR/../GUI/ \
              $$DIR/../Statistics/ \
              $$DIR/../Matrix/ \
              $$DIR/../Image/ \
              $$DIR/../Mesh/ \
              $$DIR/../MaxFlow/ \
              $$DIR/../rply/ \
              $$DIR/../MarchingCubes/ \
              $$DIR/../Carving/ \
              $$DIR/../NearestNeighbours \
              $$DIR/../KDTree \
              $$DIR/../Eigen

# Using this option with openmp fails
#DEFINES += RENDER_USING_GLUT

CONFIG(debug, debug|release): DEFINES += DEBUG
CONFIG(release, debug|release): DEFINES += RELEASE

