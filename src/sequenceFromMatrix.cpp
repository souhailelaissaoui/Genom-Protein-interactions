#include "sequenceFromMatrix.hpp"
#include "../build/ui_sequenceFromMatrix.h"

double cutoff;
bool saveM;
bool useEM;

sequenceFromMatrix::sequenceFromMatrix(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::sequenceFromMatrix)
{
    ui->setupUi(this);
    cutoff=0;
    saveM=false;
    useEM=false;
}

sequenceFromMatrix::~sequenceFromMatrix()
{
    delete ui;
}

void sequenceFromMatrix::on_spinCutoff_valueChanged(){
    cutoff = ui->spinCutoff->value();
}

double sequenceFromMatrix::getCutoff(){
    return cutoff;
}

void sequenceFromMatrix::on_buttonSave_clicked(){
    this->close();
}

void sequenceFromMatrix::on_checkBoxSaveResults_stateChanged(){
    if(ui->checkBoxSaveResults->isChecked()){
        saveM = true;
        ui->checkBoxUseEM->setEnabled(true);
    }
    else{
        saveM = false;
        useEM = false;
        ui->checkBoxUseEM->setChecked(false);
        ui->checkBoxUseEM->setEnabled(false);
    }
}

bool sequenceFromMatrix::getBool(){
    return saveM;
}

void sequenceFromMatrix::on_checkBoxUseEM_stateChanged(){
    if(ui->checkBoxSaveResults->isChecked() and ui->checkBoxUseEM->isCheckable() and ui->checkBoxUseEM->isChecked()){
        useEM = true;
    }
    else useEM = false;
}

bool sequenceFromMatrix::isEM(){
    return useEM;
}
