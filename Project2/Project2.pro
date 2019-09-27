TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt
LIBS += -larmadillo -lblas -llapack

SOURCES += \
        jacobi.cpp \
        main.cpp \
        test2.cpp

HEADERS += \
    jacobi.h \
    test2.h
