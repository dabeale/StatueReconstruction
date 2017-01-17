#ifndef MRFPARAMETERSWINDOW_H
#define MRFPARAMETERSWINDOW_H

/* ----------------------------------------------------------------------
 * Copyright (C) 2016 Daniel Beale. All rights reserved.
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 * ---------------------------------------------------------------------- */

#include <QMainWindow>
#include <QWidget>
#include <QGroupBox>
#include <QRadioButton>
#include <QLineEdit>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QLabel>
#include <QToolButton>
#include <QDoubleValidator>

namespace SegmentationGUI
{

/**
 * @brief The MRFParametersWindow class
 * This is a QMainWindow which allows the user to input the parameters for the
 * Markov random field algorithm.
 */
class MRFParametersWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MRFParametersWindow(QWidget *parent = 0);
    ~MRFParametersWindow();

    /**
     * @brief The SegType enum
     * An enumeration of Markov random field algorithms
     */
    enum SegType
    {
        GraphCut,
        ICM
    };

    SegType GetType() {return (m_radioICM->isChecked()) ? ICM : GraphCut;  }                    ///< Get the type of computation selected
    double GetAlpha() {return atof(m_alphaparam->text().toStdString().c_str());}                ///< Get the selected Alpha hyper parameter
    double GetSampleN() {return atoi(m_sampleN->text().toStdString().c_str());}                 ///< Get the selected Sample N hyper parameter
    double GetPairWiseWeight() { return atof(m_pairwiseWeight->text().toStdString().c_str());}  ///< Get the pairwise weights parameter
    double GetUnaryWeight() { return atof(m_unaryWeight->text().toStdString().c_str()); }       ///< Get the unart weight parameter

    void SetAlpha(double alpha){ m_alphaparam->setText(QString::number(alpha)); }                ///< Set the Alpha hyper parameter
    void SetSampleN(uint32_t SampleN){ m_sampleN->setText(QString::number(SampleN));}            ///< Set the sample N hyper parameter
    void SetPairWiseWeight( double weight ){m_pairwiseWeight->setText(QString::number(weight));} ///< Set the pairwise weight parameter
    void SetUnaryWeight(double weight){m_unaryWeight->setText(QString::number(weight));}         ///< Set the unary weight parameter

private slots :
    void CloseWindow();                 ///< Close the current window

signals :
    void CloseButtonPushed();           ///< Emitted when the close button has been pushed

private:
    QGroupBox *m_segtypeGroup;          ///< The segmentation type radio button group
    QRadioButton *m_radioICM;           ///< The 'ICM' algorithm radio button
    QRadioButton *m_radioGraphCut;      ///< The 'Graph Cut' radio button

    QLineEdit *m_alphaparam;            ///< The Alpha parameter input box
    QLineEdit *m_sampleN;               ///< The sample N parameter input box
    QLineEdit *m_pairwiseWeight;        ///< The pairwise weight input box
    QLineEdit *m_unaryWeight;           ///< The unary weight input box

    QWidget *m_mainwindow;              ///< The main widget

    QToolButton *m_closeButton;         ///< The close button
};
}

#endif
