LIBS += -L$$OUT_DIR/../bin
LIBS += -L/usr/local/lib

LIBS += -lopencv_highgui -lopencv_core -lopencv_imgproc -lNearestNeighbours -lKDTree -lMesh -lrply
LIBS += -lGUI -lMatrix -lStatistics -lAux -lImage -lMaxFlow -lMarchingCubes -lCarving
LIBS += -framework OpenGL -framework GLUT
