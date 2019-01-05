#ifndef matrixFromSequence_HPP
#define matrixFromSequence_HPP

#include <QDialog>


namespace Ui {
class matrixFromSequence;
}

/*!
 * @class matrixFromSequence
 *
 * @brief Display all widgets related to the calculation of a matrix from a
 * given seuquence file.
 */

class matrixFromSequence : public QDialog
{
    Q_OBJECT

public:
  /*!
   * @brief Constructor of the window.
   *
   * @param The widget affiliated to this one.
   */
    explicit matrixFromSequence(QWidget *parent = 0);
    ~matrixFromSequence();

    /*!
     * @brief Allows the interpretation of a spin box.
     *
     *
     * @return The position of the sequence.
     */
    static unsigned int getPosition();

    /*!
     * @brief Allows the interpretation of a spin box.
     *
     *
     * @return The length of the sequence.
     */
    static unsigned int getLength();

    /*!
     * @brief Allows the interpretation of a checkbox that asks the user if he
     * knows the length of the sequence.
     *
     *
     * @return Whether or not the EM algorithm should be used.
     */
    static bool getBool();

    /*!
     * @brief Allows the interpretation of a combo box.
     *
     *
     * @return An int, interprable algorithmicallyto represent the origin of the
     * file, whether it comes from a fasta, or whatnot.
     */
    static unsigned int getSequenceOrigin();

    /*!
     * @brief Allows the interpretation of a combo box.
     *
     *
     * @return An int, interprable algorithmicallyto represent the origin of the
     * file, whether it comes from a fasta, or whatnot.
     */
    static unsigned int getSequenceSource();

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
     * index choice, based on a combo box.
     */
    void on_chooseMatrix_currentIndexChanged(int);

    /*!
     * @brief Changes the value of a global parameter to correspond to the current
     * choice, based on a spin box.
     */
    void on_intPosition_valueChanged(int);

    /*!
     * @brief Changes the value of a global parameter to correspond to the current
     * choice, based on a spin box.
     */
    void on_intLength_valueChanged(int);

    /*!
     * @brief Changes a global bool value on user input.
     */
    void on_checkBox_stateChanged();

    /*!
     * @brief Changes the value of a global parameter to correspond to the current
     * index choice, based on a combo box.
     */
    void on_chooseSequenceOrigin_currentIndexChanged(int);

    /*!
     * @brief Changes the value of a global parameter to correspond to the current
     * index choice, based on a combo box.
     */
    void on_chooseSourceSequence_currentIndexChanged(int);


private:
    Ui::matrixFromSequence *ui;


};

#endif
