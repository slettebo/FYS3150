TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

LIBS += -llapack -lblas -larmadillo -lunittest++

SOURCES += main.cpp \
    lib.cpp \
    System.cpp \
    electron.cpp \
    metropolis.cpp \
    TrialWavefunction.cpp \
    NonInteractingWavefunction.cpp \
    hamiltonian.cpp \
    harmonicoscillatorwithoutcoloumb.cpp \
    interactingwavefunction.cpp \
    HarmonicOscillatorWithColoumb.cpp

HEADERS += \
    lib.h \
    System.h \
    electron.h \
    metropolis.h \
    TrialWavefunction.h \
    NonInteractingWavefunction.h \
    hamiltonian.h \
    harmonicoscillatorwithoutcoloumb.h \
    interactingwavefunction.h \
    HarmonicOscillatorWithColoumb.h

