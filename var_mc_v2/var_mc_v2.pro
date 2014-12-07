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
    Hamiltonians/HarmonicOscillatorWithColoumb.cpp \
    Hamiltonians/Hamiltonian.cpp \
    Hamiltonians/HarmonicOscillatorWithoutCoulomb.cpp \
    Random/random.cpp \
    Wavefunctions/TwoElectronsNonInteractingWavefunction.cpp \
    Wavefunctions/TwoElectronsInteractingWavefunction.cpp \
    Wavefunctions/SixElectronsNonInteractingWavefunction.cpp \
    Wavefunctions/SixElectronsInteractingWavefunction.cpp


HEADERS += \
    lib.h \
    System.h \
    Wavefunctions/TrialWavefunction.h \
    Hamiltonians/Hamiltonian.h \
    Hamiltonians/HarmonicOscillatorWithCoulomb.h \
    Hamiltonians/HarmonicOscillatorWithoutCoulomb.h \
    Random/random.h \
    TwoElectronsNonInteractingWavefunction.h \
    Wavefunctions/TwoElectronsNonInteractingWavefunction.h \
    Wavefunctions/TwoElectronsInteractingWavefunction.h \
    Wavefunctions/SixElectronsNonInteractingWavefunction.h \
    Wavefunctions/SixElectronsInteractingWavefunction.h


