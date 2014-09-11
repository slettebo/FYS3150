TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

LIBS += -llapack -lblas -larmadillo

SOURCES += main.cpp

