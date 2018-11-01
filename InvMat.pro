#-------------------------------------------------
#
# Project created by QtCreator 2018-09-04T14:48:19
#
#-------------------------------------------------

QT       += core gui xml opengl

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = InvMat
TEMPLATE = app

# The following define makes your compiler emit warnings if you use
# any feature of Qt which has been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

INCLUDEPATH += ../../graphics/InvMat/include \
../../graphics/InvMat/include/eigen3 \
../../graphics/InvMat/include/alglib

LIBPATH += ../../graphics/InvMat/libs

LIBS += -larpack_SUN4 -lQGLViewer-qt5
LIBS += -lblas -llapack -lGLEW -lGL -lGLU -lgfortran
# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

FORMS += \
    ../../graphics/InvMat/src/mainwindow.ui

HEADERS += \
    ../../graphics/InvMat/src/abqparser.h \
    ../../graphics/InvMat/src/fitmat.h \
    ../../graphics/InvMat/src/tetmesh.h \
    ../../graphics/InvMat/src/mainwindow.h \
    ../../graphics/InvMat/src/settings.h \
    ../../graphics/InvMat/src/viewer.h \
    ../../graphics/InvMat/include/alglib/alglibinternal.h \
    ../../graphics/InvMat/include/alglib/alglibmisc.h \
    ../../graphics/InvMat/include/alglib/ap.h \
    ../../graphics/InvMat/include/alglib/dataanalysis.h \
    ../../graphics/InvMat/include/alglib/diffequations.h \
    ../../graphics/InvMat/include/alglib/fasttransforms.h \
    ../../graphics/InvMat/include/alglib/integration.h \
    ../../graphics/InvMat/include/alglib/interpolation.h \
    ../../graphics/InvMat/include/alglib/linalg.h \
    ../../graphics/InvMat/include/alglib/optimization.h \
    ../../graphics/InvMat/include/alglib/solvers.h \
    ../../graphics/InvMat/include/alglib/specialfunctions.h \
    ../../graphics/InvMat/include/alglib/statistics.h \
    ../../graphics/InvMat/include/alglib/stdafx.h \
    ../../graphics/InvMat/src/shader.h

SOURCES += \
    ../../graphics/InvMat/src/abqparser.cpp \
    ../../graphics/InvMat/src/fitmat.cpp \
    ../../graphics/InvMat/src/tetmesh.cpp \
    ../../graphics/InvMat/src/main.cpp \
    ../../graphics/InvMat/src/mainwindow.cpp \
    ../../graphics/InvMat/src/viewer.cpp \
    ../../graphics/InvMat/include/alglib/alglibinternal.cpp \
    ../../graphics/InvMat/include/alglib/alglibmisc.cpp \
    ../../graphics/InvMat/include/alglib/ap.cpp \
    ../../graphics/InvMat/include/alglib/dataanalysis.cpp \
    ../../graphics/InvMat/include/alglib/diffequations.cpp \
    ../../graphics/InvMat/include/alglib/fasttransforms.cpp \
    ../../graphics/InvMat/include/alglib/integration.cpp \
    ../../graphics/InvMat/include/alglib/interpolation.cpp \
    ../../graphics/InvMat/include/alglib/linalg.cpp \
    ../../graphics/InvMat/include/alglib/optimization.cpp \
    ../../graphics/InvMat/include/alglib/solvers.cpp \
    ../../graphics/InvMat/include/alglib/specialfunctions.cpp \
    ../../graphics/InvMat/include/alglib/statistics.cpp \
    ../../graphics/InvMat/src/shader.cpp


