#-------------------------------------------------
#
# Project created by QtCreator 2015-04-16T12:28:44
#
#-------------------------------------------------

QT       += core
greaterThan(QT_MAJOR_VERSION, 4): QT += widgets
MODULENAME = MarchingCubes
OUT_DIR=$$OUT_PWD
DIR=$$PWD

include(../common.pri)

#TEMPLATE = app
TEMPLATE = lib
CONFIG += staticlib

SOURCES += marchingcubes.cpp \
    mlists.cpp \
    render.cpp

HEADERS += \ 
    marchingcubes.h \
    render.h

