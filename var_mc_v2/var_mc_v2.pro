TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    lib.cpp \
    System.cpp \
    electron.cpp \
    metropolis.cpp \
    TrialWavefunction.cpp \
    NonInteractingWavefunction.cpp

HEADERS += \
    lib.h \
    System.h \
    electron.h \
    metropolis.h \
    TrialWavefunction.h \
    NonInteractingWavefunction.h
