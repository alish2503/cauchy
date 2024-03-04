#include "thirdwindow.h"
#include "ui_thirdwindow.h"
#include <QPixmap>

thirdwindow::thirdwindow(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::thirdwindow)
{
    ui->setupUi(this);
    resize(QDesktopWidget().availableGeometry(this).size());
}

thirdwindow::~thirdwindow()
{
    delete ui;
}

void thirdwindow::on_pushButton_2_clicked() {
    this->close();
    emit mainwindow();
}

void thirdwindow::slot(QString name, QString coors) {
    ui->label_3->setText(name);
    ui->label_4->setText(coors);
    QPixmap pix("graph.png");
    int w = ui->label_5->width();
    int h = ui->label_5->height();
    ui->label_5->setPixmap(pix.scaled(w,h, Qt::KeepAspectRatio));
}

void thirdwindow::on_pushButton_clicked()
{
    remove("graph.png");
    remove("data.txt");
    QApplication::quit();
}
