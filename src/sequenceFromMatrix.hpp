#ifndef SEQUENCEFROMMATRIX_HPP
#define SEQUENCEFROMMATRIX_HPP

#include <QDialog>

#include "resultsWindow.hpp"
#include "procedures.hpp"

namespace Ui {
class sequenceFromMatrix;
}

/*!
 * @class sequenceFromMatrix
 *
 * @brief Display all widgets related to the calculation of a matrix from a
 * given seuquence file.
 */

class sequenceFromMatrix : public QDialog
{
    Q_OBJECT

public:
    explicit sequenceFromMatrix(QWidget *parent = 0);
    ~sequenceFromMatrix();

    /*! @brief Change a global parameter to an interpretable value.
     *
     *
     * @return Returns the cutoff value set by the user.
     */
    static double getCutoff();

    /*!
     * @brief Allows the interpretation of a checkbox that asks the user if he
     * wants to save a matrix from the sequence obtained
     *
     *
     * @return Whether or not a resulting matrix should be saved.
     */
    static bool getBool();

    /*!
     * @brief Transfer a global value to a usable one.
     *
     *
     * @return A bool for whether or not the EM algorithm shall be used.
     */
    static bool isEM();

public slots:
    /*!
     * @brief When this button will be clicked, the program will check if the exit
     * conditions are respected and will react accordingly.
     */
    void on_buttonSave_clicked();

    /*!
     * @brief Changes the value of a global parameter to correspond to the current
     * choice, based on a spin box.
     */
    void on_spinCutoff_valueChanged();

    /*!
     * @brief Changes a global bool value on user input.
     */
    void on_checkBoxSaveResults_stateChanged();

    /*!
     * @brief Changes a global bool value on user input.
     */
    void on_checkBoxUseEM_stateChanged();

private:
    Ui::sequenceFromMatrix *ui;
};

#endif
