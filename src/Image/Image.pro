#-------------------------------------------------
#
# Project created by QtCreator 2015-04-16T12:28:44
#
#-------------------------------------------------

QT       += core
MODULENAME = Image
OUT_DIR=$$OUT_PWD
DIR=$$PWD

include (../common.pri)

#TEMPLATE = app
TEMPLATE = lib
CONFIG += staticlib

SOURCES += \
    Image.cpp \
    libppm.cpp \
    Video.cpp \
    Colour.cpp

HEADERS += \ 
    Image.h \
    libppm.h \
    Video.h \
    Colour.h

LIBS +=
