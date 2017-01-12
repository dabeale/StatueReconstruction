#-------------------------------------------------
#
# Project created by QtCreator 2015-04-16T12:28:44
#
#-------------------------------------------------

QT       += core
greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

MODULENAME = Statistics
OUT_DIR=$$OUT_PWD
DIR=$$PWD

include(../common.pri)

CONFIG -= app_bundle
#TEMPLATE = app
TEMPLATE = lib
CONFIG += staticlib

SOURCES += GaussianMixture.cpp \
    Sample.cpp \
    InfiniteGMM.cpp \
    Distributions.cpp

HEADERS += GaussianMixture.h \
    Sample.h \
    InfiniteGMM.h \
    Distributions.h

