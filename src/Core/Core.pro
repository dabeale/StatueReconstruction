#-------------------------------------------------
#
# Project created by QtCreator 2014-04-29T11:07:15
#
#-------------------------------------------------

QT       += core
MODULENAME = Core
OUT_DIR=$$OUT_PWD
DIR=$$PWD

# Include shared options
include (../common.pri)

TEMPLATE = lib
#CONFIG += dynamiclib
CONFIG += staticlib

SOURCES += \
    Sample.cpp \
    Matrix.cpp \
    Message.cpp \
    Distributions.cpp

HEADERS += \
    Matrix.h \
    Sample.h \
    Message.h \
    Distributions.h


