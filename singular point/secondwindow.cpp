#include "secondwindow.h"
#include "ui_secondwindow.h"
#include <QPixmap>

secondwindow::secondwindow(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::secondwindow)
{
    ui->setupUi(this);
    resize(QDesktopWidget().availableGeometry(this).size());
}

secondwindow::~secondwindow()
{
    delete ui;
}

void secondwindow::on_pushButton_2_clicked()
{
    this->close();
    emit mainwindow();
}


void secondwindow::slot() {
    QPixmap pix("graph.png");
    int w = ui->label->width();
    int h = ui->label->height();
    ui->label->setPixmap(pix.scaled(w,h, Qt::KeepAspectRatio));
}

void secondwindow::on_pushButton_clicked()
{
    remove("graph.png");
    remove("data.txt");
    QApplication::quit();
}
