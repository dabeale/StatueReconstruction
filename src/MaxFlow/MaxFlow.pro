#-------------------------------------------------
#
# Project created by QtCreator 2015-04-16T12:28:44
#
#-------------------------------------------------

QT       += core
greaterThan(QT_MAJOR_VERSION, 4): QT += widgets
MODULENAME = MaxFlow
OUT_DIR=$$OUT_PWD
DIR=$$PWD

include(../common.pri)

#TEMPLATE = app
TEMPLATE = lib
CONFIG += staticlib

SOURCES += \ 
    graph.cpp \
    maxflow.cpp \
    evaluate.cpp

HEADERS += \ 
    block.h \
    graph.h \
    maxflow.h
