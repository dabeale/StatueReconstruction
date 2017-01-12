#-------------------------------------------------
#
# Project created by QtCreator 2014-07-03T20:16:55
#
#-------------------------------------------------

QT       += core
QT       -= gui
MODULENAME = Aux
OUT_DIR=$$OUT_PWD
DIR=$$PWD

include(../common.pri)

#TEMPLATE = app
TEMPLATE = lib
CONFIG += staticlib

SOURCES += Message.cpp

HEADERS += Message.h \
           StringFunctions.h

