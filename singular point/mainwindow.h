#ifndef MAINWINDOW_H
#define MAINWINDOW_H
#include <QMainWindow>
#include <QMessageBox>
#include <QDesktopWidget>
#include <chrono>
#include <thread>
#include <cmath>
#include <fstream>
#include "secondwindow.h"
#include "thirdwindow.h"

typedef std::vector<double> dv;

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

signals:
    void signal(QString, QString);
    void signal1();

private slots:
    void on_pushButton_clicked();
    void on_pushButton_2_clicked();

private:
    Ui::MainWindow *ui;
    secondwindow *w2;
    thirdwindow *w3;
};
#endif // MAINWINDOW_H
