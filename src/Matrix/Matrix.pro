#-------------------------------------------------
#
# Project created by QtCreator 2014-07-03T20:16:55
#
#-------------------------------------------------

QT       += core
QT       -= gui
MODULENAME = Matrix
OUT_DIR=$$OUT_PWD
DIR=$$PWD

include(../common.pri)

#TEMPLATE = app
TEMPLATE = lib
CONFIG += staticlib

SOURCES += Matrix.cpp \
    Utils.cpp

HEADERS += Matrix.h \
    Utils.h

