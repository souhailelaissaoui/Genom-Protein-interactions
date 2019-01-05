#include <QApplication>
#include <stdio.h>
#include <string>
#include <window.hpp>


#include "Matrix.hpp"
#include "procedures.hpp"
#include "user_interaction.hpp"
#include "utility.hpp"



int main(int argc, char *argv[]) {
    QApplication app(argc, argv);

    Window window;
    window.show();
    return app.exec();
}
