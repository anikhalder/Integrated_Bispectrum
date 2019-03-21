TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

LIBS += -L/home/anik/anaconda3/pkgs/gsl/lib -lgsl -lgslcblas -lm

INCLUDEPATH += $$_PRO_FILE_PWD_
INCLUDEPATH += /home/anik/anaconda3/pkgs/gsl/include

SOURCES += \
        main.cpp
