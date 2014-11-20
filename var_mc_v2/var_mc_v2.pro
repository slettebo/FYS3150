TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

LIBS += -llapack -lblas -larmadillo -lunittest++ -fopenmp

QMAKE_LFLAGS += -fopenmp
QMAKE_CXXFLAGS += -fopenmp


SOURCES += main.cpp \
    lib.cpp \
    System.cpp \
    Wavefunctions/TrialWavefunction.cpp \
    Wavefunctions/NonInteractingWavefunction.cpp \
    Wavefunctions/InteractingWavefunction.cpp \
    Hamiltonians/HarmonicOscillatorWithColoumb.cpp \
    Hamiltonians/Hamiltonian.cpp \
    Hamiltonians/HarmonicOscillatorWithoutCoulomb.cpp


HEADERS += \
    lib.h \
    System.h \
    Wavefunctions/TrialWavefunction.h \
    Wavefunctions/NonInteractingWavefunction.h \
    Wavefunctions/InteractingWavefunction.h \
    Hamiltonians/Hamiltonian.h \
    Hamiltonians/HarmonicOscillatorWithCoulomb.h \
    Hamiltonians/HarmonicOscillatorWithoutCoulomb.h


