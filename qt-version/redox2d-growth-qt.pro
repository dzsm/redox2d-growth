#-------------------------------------------------
#
# Project created by QtCreator 2015-03-11T22:23:19
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport

TARGET = growthX
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp \
    drawingarea.cpp \
    qcustomplot.cpp

HEADERS  += mainwindow.h \
    drawingarea.h \
    simulation.h \
    growth_base.h \
    qcustomplot.h

FORMS    += mainwindow.ui

CONFIG += c++11
QMAKE_CXXFLAGS += -O3


INCLUDEPATH += /usr/include/eigen3
