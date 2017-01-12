#-------------------------------------------------
#
# Project created by QtCreator 2015-04-16T12:28:44
#
#-------------------------------------------------

QT       += core gui
greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

MODULENAME = GUI
OUT_DIR=$$OUT_PWD
DIR=$$PWD

include(../common.pri)

#TEMPLATE = app
TEMPLATE = lib
CONFIG += staticlib

SOURCES += MainWindow.cpp \
  ImageViewer.cpp \
  ImageInfo.cpp \
  PenPoint.cpp \
    SegmentationModel.cpp \
    MRFParametersWindow.cpp

HEADERS += MainWindow.h \
  ImageViewer.hpp \
  ImageInfo.h \
  PenPoint.h \
    SegmentationModel.h \
    MRFParametersWindow.h
