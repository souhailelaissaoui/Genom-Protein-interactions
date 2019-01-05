#ifndef RESULTSWINDOW_HPP
#define RESULTSWINDOW_HPP

#include <QDialog>
#include <QTextBrowser>

namespace Ui {
class resultsWindow;
}

/*!
 * @class resultsWindow
 *
 * @brief Display the final results in an independant widget.
 */

class resultsWindow : public QDialog
{
    Q_OBJECT

public:
    explicit resultsWindow(QWidget *parent = 0);
    ~resultsWindow();

public slots:
/*!
 * @brief Close the window.
 */
    void on_buttonLeave_clicked();


private:
    Ui::resultsWindow *ui;
    std::string resultsLocation;
    QTextBrowser resultsBrowser;
};

#endif
