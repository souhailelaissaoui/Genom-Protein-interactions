#include "askBaseProb.hpp"
#include "ui_askBaseProb.h"

unsigned int matrixChoiceS;
double probA;
double probC;
double probG;
double probT;
double total;
unsigned int baseChoice;

askBaseProb::askBaseProb(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::askBaseProb)
{
    ui->setupUi(this);
    matrixChoiceS =1;
    probA=0;
    probC=0;
    probG=0;
    probT=0;
    baseChoice=0;
    if (sequenceFromMatrix::getBool()==false) ui->chooseMatrix->setEnabled(false);
}

askBaseProb::~askBaseProb()
{
    delete ui;
}


void askBaseProb::on_chooseBase_currentIndexChanged(int){
    if(ui->chooseBase->currentIndex()==3){
        ui->spinA->setEnabled(true);
        ui->spinC->setEnabled(true);
        ui->spinG->setEnabled(true);
        ui->spinT->setEnabled(true);
        ui->showTotal->setEnabled(true);
        ui->label_3->setEnabled(true);
    }
    else{
        ui->spinA->setEnabled(false);
        ui->spinC->setEnabled(false);
        ui->spinG->setEnabled(false);
        ui->spinT->setEnabled(false);
        ui->showTotal->setEnabled(false);
        ui->label_3->setEnabled(false);
    }
    baseChoice = ui->chooseBase->currentIndex();
}

unsigned int askBaseProb::getBaseChoiceFinal(){
    if (baseChoice == 0) return 1;
    else if (baseChoice == 1) return 0;
    else if(baseChoice ==2) return 3;
    else if(baseChoice==3) return 2;
    else return 1;
}

double askBaseProb::getProbA() {
    return probA;
}
double askBaseProb::getProbC(){
    return probC;
}

double askBaseProb::getProbG(){
    return probG;
}

double askBaseProb::getProbT(){
    return probT;
}

void askBaseProb::on_spinA_valueChanged(double){
    probA = ui->spinA->value();
    setTotal();
}

void askBaseProb::on_spinC_valueChanged(double){
    probC = ui->spinC->value();
    setTotal();
}

void askBaseProb::on_spinG_valueChanged(double){
    probG = ui->spinG->value();
    setTotal();
}

void askBaseProb::on_spinT_valueChanged(double){
    probT = ui->spinT->value();
    setTotal();
}

void askBaseProb::on_chooseMatrix_currentIndexChanged(int){
    matrixChoiceS = ui->chooseMatrix->currentIndex()+1;
}

unsigned int askBaseProb::getMatrixChoice(){
    return matrixChoiceS;
}

void askBaseProb::setTotal(){
    total = probA + probC + probG + probT;
    ui->showTotal->setValue(total);
}

bool askBaseProb::checkCustomTotal(){
    if (total == 1) return true;
    else return false;
}

void askBaseProb::on_buttonSave_clicked(){
    if (ui->chooseBase->currentIndex()!=3){
        this->close();
    }
    else if (ui->chooseBase->currentIndex()==3 and checkCustomTotal()){
    this->close();
    }
    //QT TO DO EM IMPLEMENTATION
    else if (ui->chooseBase->currentIndex()==2) QMessageBox::information(this, "No implementation", "Please choose another base, the EM algorithm has not been implemented yet");
    else{
        QMessageBox::information(this, "Invalid Base", "Your bases must add up to 1.00 .");
    }
}
