#ifndef ASKBASEPROB_HPP
#define ASKBASEPROB_HPP

#include <QDialog>
#include <QMessageBox>

#include "sequenceFromMatrix.hpp"

namespace Ui {
class askBaseProb;
}

/*!
 * @class askBaseProb
 *
 * @brief Display a window that will ask the user its choice of base probability,
 * and its choice in matrix type, if the choice is open, based on pre-conditions.
 */

class askBaseProb : public QDialog

{
    Q_OBJECT

public:
  /*!
   * @brief Constructor of the window.
   *
   * @param The widget affiliated to this one.
   */
    explicit askBaseProb(QWidget *parent = 0);
    ~askBaseProb();

    /*!
     * @brief Allows the interpretation of a combo box.
     *
     *
     * @return The final choice of base as an int, allowing compatibility with
     * the terminal version.
     */
    static unsigned int getBaseChoiceFinal();

    /*!
     * @brief Returns the value entered by the user in the interface, when
     * the user chose to use a custom base probability.
     *
     * @return Returns the probability of Adenine.
     */
    static double getProbA();

    /*!
     * @brief Returns the value entered by the user in the interface, when
     * the user chose to use a custom base probability.
     *
     * @return Returns the probability of Cytosine.
     */
    static double getProbC();

    /*!
     * @brief Returns the value entered by the user in the interface, when
     * the user chose to use a custom base probability.
     *
     * @return Returns the probability of Guanine.
     */
    static double getProbG();

    /*!
     * @brief Returns the value entered by the user in the interface, when
     * the user chose to use a custom base probability.
     *
     * @return Returns the probability of Thymine.
     */
    static double getProbT();

    /*!
     * @brief Returns the user matrix choice, which depends on a combo box.
     *
     * @return Returns the probability of Adenine.
     */
    static unsigned int getMatrixChoice();

    /*!
     * @brief Checks if the sum of all 4 bases is equal to one.
     *
     * @return The boolean will be true only if the sum of the 4 bases is equal.
     */
    bool checkCustomTotal();

    /*!
     * @brief Modify the interface value of the sum of the 4 bases so that
     * it is easier for the user to see the result.
     */
    void setTotal();

public slots:
    /*!
     * @brief Changes the value of a global parameter to correspond to the current
     * index choice, based on a combo box.
     */
    void on_chooseBase_currentIndexChanged(int);

    /*!
     * @brief Changes the value of a global parameter to correspond to the current
     * choice, based on a spin box.
     */
    void on_spinA_valueChanged(double);

    /*!
     * @brief Changes the value of a global parameter to correspond to the current
     * choice, based on a spin box.
     */
    void on_spinC_valueChanged(double);

    /*!
     * @brief Changes the value of a global parameter to correspond to the current
     * choice, based on a spin box.
     */
    void on_spinG_valueChanged(double);

    /*!
     * @brief Changes the value of a global parameter to correspond to the current
     * choice, based on a spin box.
     */
    void on_spinT_valueChanged(double);

    /*!
     * @brief Changes the value of a global parameter to correspond to the current
     * choice, based on a combo box.
     */
    void on_chooseMatrix_currentIndexChanged(int);

    /*!
     * @brief When this button will be clicked, the program will check if the exit
     * conditions are respected and will react accordingly.
     */
    void on_buttonSave_clicked();

private:
    Ui::askBaseProb *ui;
};

#endif
