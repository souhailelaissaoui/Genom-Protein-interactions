QT += widgets

TARGET = genom
TEMPLATE = app

HEADERS += \
    Matrix.hpp \
    procedures.hpp \
    user_interaction.hpp \
    utility.hpp \
    window.hpp \
    matrixFromSequence.hpp \
    sequenceFromMatrix.hpp \
    resultsWindow.hpp \
    genomic_coordinates.hpp \
    askBaseProb.hpp \
    logo.hpp

SOURCES += \
    mainQt.cpp \
    Matrix.cpp \
    procedures.cpp \
    user_interactionQt.cpp \
    utility.cpp \
    window.cpp \
    matrixFromSequence.cpp \
    sequenceFromMatrix.cpp \
    resultsWindow.cpp \
    genomic_coordinates.cpp \
    askBaseProb.cpp

FORMS += \
    window.ui \
    matrixFromSequence.ui \
    sequenceFromMatrix.ui \
    resultsWindow.ui \
    askBaseProb.ui


INCLUDEPATH += ../logo

LIBS += -L$PATH_TO_CIMG_LIB -lpng -lpthread
