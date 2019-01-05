#ifndef WINDOW_HPP
#define WINDOW_HPP

#include <QApplication>
#include <QColorDialog>
#include <QDialog>
#include <QFileDialog>
#include <QFontDialog>
#include <QInputDialog>
#include <QMessageBox>
#include <QObject>
#include <QPushButton>
#include <QSlider>
#include <QWidget>

#include "askBaseProb.hpp"
#include "procedures.hpp"
#include "matrixFromSequence.hpp"
#include "sequenceFromMatrix.hpp"
#include "utility.hpp"


namespace Ui {
  class Window;
}

/*!
 * @class Window
 *
 * @brief Initial widget of Qt.
 */

class Window : public QDialog
{
  Q_OBJECT

//-----------------------------------------------------------------------
public:
  /*!
  * @brief Constructor
  */
  explicit Window(QWidget *parent = 0);

  /*!
  * @brief Destructor
  */
  ~Window();

  /*!
  * @brief Returns the latest Fasta location specified by the user.
  */
  static std::string getFastaLocation();

  /*!
  * @brief Returns the latest Matrix location specified by the user.
  */
  static std::string getMatrixLocation();

  /*!
  * @brief Returns the latest location and name of the file the user wants
  * the results to be saved to.
  */
  static std::string getOutputName();

//-----------------------------------------------------------------------
  public slots:
  /*!
  * @brief When this button is clicked, it opens a dialog that allows the user
  * to choose a fasta file.
  */
  void on_buttonFasta_clicked();

  /*!
  * @brief When this button is clicked, it opens a dialog that allows the
  * user to choose a matrix file.
  */
  void on_buttonMatrix_clicked();

  /*!
  * @brief When this button is clicked, the program closes.
  */
  void on_buttonLeave_clicked();

  /*!
  * @brief When this button is clicked, the program will start the
  * algorithms to obtain the binding positions of a protein from a
  * probability matrix.
  */
  void on_buttonSequenceFromMatrix_clicked();

  /*!
  * @brief When this button is clicked, the program will start the
  * algorithms to obtain a probability weigth matrix from sequences
  * obtained from a fasta file.
  */
  void on_buttonMatrixFromSequence_clicked();

  /*!
  * @brief When this button is clicked, the program will lainch the
  * algorithms to create a logo based on a probability weigth matrix.
  */
  void on_buttonLogo_clicked();

  /*!
  * @brief When this button is clicked, the program will launch the
  * algorithm to  compare the score of sequences to the score
  * in genomic coordinates.
  */
  void on_buttonCorrelate_clicked();

//-----------------------------------------------------------------------
private:
  Ui::Window *ui;

  QString fastaLocation;
  QString matrixLocation;
  std::string fastaLocationFromWindow;
  std::string matrixLocationFromWindow;
};

#endif
