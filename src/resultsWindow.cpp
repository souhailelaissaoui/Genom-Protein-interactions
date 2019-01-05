#include "resultsWindow.hpp"
#include "ui_resultsWindow.h"

#include <QTextStream>
#include <iostream>

#include "window.hpp"

resultsWindow::resultsWindow(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::resultsWindow)
{
    ui->setupUi(this);

    //The output file is saved one directory before, thus the addtiion of ../
    QString filename = QString::fromStdString(Ask_Outputfile_Name()+".csv");
    QFile file(filename);

    file.open(QFile::ReadOnly | QFile::Text);
    QTextStream ReadFile(&file);
    ui->textBrowser->setText(ReadFile.readAll());
}

resultsWindow::~resultsWindow()
{
    delete ui;
}

void resultsWindow::on_buttonLeave_clicked(){
    this->close();
}
