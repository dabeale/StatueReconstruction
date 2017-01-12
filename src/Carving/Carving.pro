#-------------------------------------------------
#
# Project created by QtCreator 2014-04-29T11:07:15
#
#-------------------------------------------------

QT       -= gui
MODULENAME = Carving
OUT_DIR=$$OUT_PWD
DIR=$$PWD

include(../common.pri)

TEMPLATE = lib
CONFIG += staticlib

SOURCES += Carving.cpp

HEADERS += Carving.h

CONFIG(debug, debug|release): DEFINES += DEBUG
CONFIG(release, debug|release): DEFINES += RELEASE
