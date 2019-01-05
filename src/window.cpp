#include "window.hpp"
#include "../build/ui_window.h"


#include "logo.hpp"

std::string matrixFilePath;
std::string fastaFilePath;
std::string output;

Window::Window(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::Window)
{
    ui->setupUi(this);
    matrixFilePath="";
    fastaFilePath="";
    output="";
}

Window::~Window() {
    delete ui;
}

void Window::on_buttonLeave_clicked() {
    qApp->quit();
}

void Window::on_buttonFasta_clicked() {
    fastaLocation = QFileDialog::getOpenFileName(this, "Open a Sequence File", fastaLocation, "Fasta files only! ( *.fasta *.fas *.ffn *.faa *.frn )");

    ui->showFastaLocation->setText(fastaLocation);
    ui->showFastaLocation->adjustSize();
}

void Window::on_buttonMatrix_clicked() {
    matrixLocation = QFileDialog::getOpenFileName(this, "Open a Matrix File", matrixLocation, "(*.mat)");
    ui->showMatrixLocation->setText(matrixLocation);
    ui->showMatrixLocation->adjustSize();
}

void Window::on_buttonSequenceFromMatrix_clicked() {
    if(fastaLocation.isEmpty()) QMessageBox::critical(this, "Error: Files not chosen", "You haven't chosen your sequence file.");
    else if(matrixLocation.isEmpty()) QMessageBox::critical(this, "Error: Files not chosen", "You haven't chosen your matrix file.");
    else if(ui->editFileName->text().isEmpty()) QMessageBox::critical(this, "Error: Files not chosen", "You havent chosen your output name.");
    else {
        matrixFilePath = matrixLocation.toStdString();
        fastaFilePath = fastaLocation.toStdString();
        output = ui->editFileName->text().toStdString();

        sequenceFromMatrix mac;
        mac.show();
        mac.exec();

        askBaseProb base;
        base.show();
        base.exec();

        enzyme_on_sequence();

        resultsWindow results;
        results.show();
        results.exec();
    }
}

void Window::on_buttonMatrixFromSequence_clicked() {
    if(fastaLocation.isEmpty()) {
        QMessageBox::critical(this, "Error: File not chosen", "You haven't chosen your Sequence file.");
    }
    else if(ui->editFileName->text().isEmpty()) {
        QMessageBox::critical(this, "Error: Files not chosen", "You haven't chosen your output name.");
    }
    else {
    fastaFilePath = fastaLocation.toStdString();
    output = ui->editFileName->text().toStdString();

    matrixFromSequence glass;
    glass.show();
    glass.exec();

    enzyme_from_sequences();
    }
}

void Window::on_buttonLogo_clicked(){
  if(matrixLocation.isEmpty()) {
      QMessageBox::critical(this, "Error: File not chosen", "You haven't chosen your Matrix file.");
  }
  else{
    matrixFilePath = matrixLocation.toStdString();
    QMessageBox::information(this, "Generation of the image", "Your logo will be saved in genom-2 with the name yourlogo.png");
    logo();
  }
}

void Window::on_buttonCorrelate_clicked(){

}

std::string Window::getFastaLocation(){
    return fastaFilePath;
}

std::string Window::getMatrixLocation(){
    return matrixFilePath;
}

std::string Window::getOutputName(){
    return "../" + output;
}
