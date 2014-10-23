TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

LIBS += -llapack -lblas -larmadillo -lunittest++


SOURCES += main.cpp lib.cpp
HEADERS += lib.h

