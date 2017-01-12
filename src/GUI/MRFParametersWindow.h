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

class MRFParametersWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MRFParametersWindow(QWidget *parent = 0);
    ~MRFParametersWindow();

    enum SegType
    {
        GraphCut,
        ICM
    };

    SegType GetType() {return (m_radioICM->isChecked()) ? ICM : GraphCut;  }
    double GetAlpha() {return atof(m_alphaparam->text().toStdString().c_str());}
    double GetSampleN() {return atoi(m_sampleN->text().toStdString().c_str());}
    double GetPairWiseEight() { return atof(m_pairwiseWeight->text().toStdString().c_str());}
    double GetUnaryWeight() { return atof(m_unaryWeight->text().toStdString().c_str()); }

    void SetAlpha(double alpha){ m_alphaparam->setText(QString::number(alpha)); }
    void SetSampleN(uint32_t SampleN){ m_sampleN->setText(QString::number(SampleN));}
    void SetPairWiseWeight( double weight ){m_pairwiseWeight->setText(QString::number(weight));}
    void SetUnaryWeight(double weight){m_unaryWeight->setText(QString::number(weight));}

private slots :
    void CloseWindow();

signals :
    void CloseButtonPushed();

private:
    QGroupBox *m_segtypeGroup;
    QRadioButton *m_radioICM;
    QRadioButton *m_radioGraphCut;

    QLineEdit *m_alphaparam;
    QLineEdit *m_sampleN;
    QLineEdit *m_pairwiseWeight;
    QLineEdit *m_unaryWeight;

    QWidget *m_mainwindow;

    QToolButton *m_closeButton;
};

#endif
