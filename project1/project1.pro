TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp

LIBS += -llapack -lblas -larmadillo

OTHER_FILES += \
    plot.py \
    rapport.tex \
    plot_numerical_analytical.py \
    plot_rel_error.py
