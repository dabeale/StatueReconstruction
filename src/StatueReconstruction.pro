TEMPLATE = subdirs

SUBDIRS += \
    GUI \
    Statistics \
    Matrix \
    Aux \
    Image \
    Mesh \
    rply \
    MaxFlow \
    MarchingCubes \
    Carving \
    KDTree \
    NearestNeighbours \
    Executables

Executables.depends = GUI Statistics Matrix Aux Image Mesh KDTree NearestNeighbours Carving MarchingCubes MaxFlow rply
