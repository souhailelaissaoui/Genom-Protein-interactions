#include "matrixFromSequence.hpp"
#include "../build/ui_matrixFromSequence.h"

#include <QFileDialog>
#include <iostream>
#include <QMessageBox>


bool checkBox;
bool useEM2;
unsigned int matrixChoice;
unsigned int position;
unsigned int length;
unsigned int sequenceOrigin;
unsigned int sequenceSource;
matrixFromSequence::matrixFromSequence(QWidget *parent) :
    QDialog(parent), ui(new Ui::matrixFromSequence)
{
    ui->setupUi(this);
    checkBox=false;
    useEM2=false;
    matrixChoice = 1;
    position=1;
    length=0;
    sequenceOrigin=1;
    sequenceSource=1;
}

matrixFromSequence::~matrixFromSequence() {
    delete ui;
}

void matrixFromSequence::on_checkBox_stateChanged() {
    if(ui->checkBox->isChecked()){
        ui->intLength->setEnabled(true);
        checkBox = !checkBox;
    }
    else{
        ui->intLength->setEnabled(false);
        checkBox = !checkBox;
    }
}

void matrixFromSequence::on_chooseMatrix_currentIndexChanged(int) {
    matrixChoice = ui->chooseMatrix->currentIndex()+1;
}

void matrixFromSequence::on_intPosition_valueChanged(int) {
    position = ui->intPosition->value();
}

void matrixFromSequence::on_intLength_valueChanged(int) {
    length = ui->intLength->value();
}

void matrixFromSequence::on_chooseSequenceOrigin_currentIndexChanged(int){
    sequenceOrigin = ui->chooseSequenceOrigin->currentIndex()+1;
}

void matrixFromSequence::on_chooseSourceSequence_currentIndexChanged(int){
    sequenceSource = ui->chooseSourceSequence->currentIndex()+1;
    if (sequenceSource == 2){
        ui->chooseSequenceOrigin->setEnabled(false);
    }
    else if (sequenceSource ==1){
        ui->chooseSequenceOrigin->setEnabled(true);
    }
}

unsigned int matrixFromSequence::getPosition() {
    return position;
}

unsigned int matrixFromSequence::getLength() {
    return length;
}

bool matrixFromSequence::getBool() {
    return checkBox;
}

unsigned int matrixFromSequence::getSequenceOrigin(){
    return sequenceOrigin;
}

unsigned int matrixFromSequence::getSequenceSource(){
    return sequenceSource;
}

void matrixFromSequence::on_buttonSave_clicked() {
    if( length == 0 ){
        int choice = QMessageBox::question(this, "EM Algorithm", "You haven't chosen a lenght value. Press Ok if you want to use EM Algorithm. Press Cancel if you wish to modify this.", QMessageBox::Ok | QMessageBox::Cancel);

        if (choice == QMessageBox::Ok){
            useEM2=true;
            this->close();
        }
        else if (choice == QMessageBox::Cancel){}
    }
    else if( position != 0 ) this->close();
}

bool matrixFromSequence::isEM(){
    return useEM2;
}
