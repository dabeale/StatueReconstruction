#-------------------------------------------------
#
# Project created by QtCreator 2015-04-16T12:28:44
#
#-------------------------------------------------

QT       += core gui
greaterThan(QT_MAJOR_VERSION, 4): QT += widgets
MODULENAME = KDTree
OUT_DIR=$$OUT_PWD
DIR=$$PWD

include (../common.pri)

#TEMPLATE = app
TEMPLATE = lib
CONFIG += staticlib

SOURCES += node.cpp \
    tree.cpp

HEADERS += \
    node.h \
    tree.h
