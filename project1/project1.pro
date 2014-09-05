TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp

LIBS += -llapack -lblas -larmadillo

OTHER_FILES += \
    exec_cmd.py \
    plot.py
