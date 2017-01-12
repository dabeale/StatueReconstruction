#-------------------------------------------------
#
# Project created by QtCreator 2015-04-16T12:28:44
#
#-------------------------------------------------

QT       += core
greaterThan(QT_MAJOR_VERSION, 4): QT += widgets
MODULENAME = rply
OUT_DIR=$$OUT_PWD
DIR=$$PWD

include(../common.pri);

#TEMPLATE = app
TEMPLATE = lib
CONFIG += staticlib

SOURCES += \  
    rply.c \
    etc/convert.c \
    etc/dump.c \
    etc/sconvert.c
HEADERS += \  
    rply.h \
    rplyfile.h

#QMAKE_CXXFLAGS+= -Wgnu-array-member-paren-init

