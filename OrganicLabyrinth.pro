#-------------------------------------------------
#
# Project created by QtCreator 2015-06-05T17:04:24
#
#-------------------------------------------------

QT       += core gui opengl svg

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = OrganicLabyrinth
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp \
    GLContainer.cpp \
    GLWidget.cpp

HEADERS  += mainwindow.h \
    GLContainer.h \
    GLWidget.h

FORMS    += mainwindow.ui
