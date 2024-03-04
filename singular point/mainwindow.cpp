#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QDebug>

double A, B;

double f1(double y1, double y2) {return A - y1 * y2;}
double f2(double y1, double y2) {return y1 * y2 - y2 / (B + y2);}

void runge_kutta_method(dv &, dv &, dv &);
void calc_next(double, double, double steps[], double &, double &);
void grid(dv &t, size_t n, double step) {for (size_t i = 1; i < n; i++) t[i] = step * i;}
const char* singular_point_name(double, double, double);
void write_data(const dv&, const dv&);
FILE * call_gnuplot(const char *, const char *);

MainWindow::MainWindow(QWidget *parent): QMainWindow(parent), ui(new Ui::MainWindow) {
    ui->setupUi(this);
    resize(QDesktopWidget().availableGeometry(this).size());
    w2 = new secondwindow(); w3 = new thirdwindow();
    connect(w2, &secondwindow::mainwindow, this, &MainWindow::show);
    connect(this, &MainWindow::signal1, w2, &secondwindow::slot);
    connect(w3, &thirdwindow::mainwindow, this, &MainWindow::show);
    connect(this, &MainWindow::signal, w3, &thirdwindow::slot);
}

MainWindow::~MainWindow() {delete ui;}

void MainWindow::on_pushButton_clicked() {
    size_t n = ui->lineEdit->text().toInt();
    double step = ui->lineEdit_2->text().toDouble();
    double eps = ui->lineEdit_3->text().toDouble();
    A = ui->lineEdit_4->text().toDouble();
    B = ui->lineEdit_5->text().toDouble();
    bool error = false;
    if (!ui->lineEdit->text().isEmpty() && !ui->lineEdit_2->text().isEmpty() &&
        !ui->lineEdit_3->text().isEmpty() && !ui->lineEdit_4->text().isEmpty() &&
        !ui->lineEdit_5->text().isEmpty()) {
        if (A == 1.0 || !B || eps <= 0.0 || n <= 0 || step <= 0.0) {
            QMessageBox::warning(this, " ", "Ошибка: Некорректные значения параметров");
            error = true;
        }
        double Y1 = (1.0 - A) / B, Y2 = A * B / (1.0 - A);
        double y10 = Y1 + eps, y20 = Y2 + eps;
        if (y10 < 0.0 || y20 < 0.0) {
            QMessageBox::warning(this, " ", "Ошибка: начальная точка вне 1-ой четверти");
            error = true;
        }
        if (Y1 < 0.0 || Y2 < 0.0) {
            QMessageBox::warning(this, " ", "Ошибка: стационарная точка вне 1-ой четверти");
            error = true;
        }
        if (!error) {
            dv t(n, 0.0), y1(n, y10), y2(n, y20);
            grid(t, n, step);
            runge_kutta_method(y1, y2, t);
            FILE *gp = NULL;
            if (ui->radioButton->isChecked() || ui->radioButton_2->isChecked()) {
                if (ui->radioButton->isChecked()) {
                    write_data(t, y1);
                    gp = call_gnuplot("set xlabel 't'\n", "set ylabel 'y1'\n");
                }
                else {
                    write_data(t, y2);
                    gp = call_gnuplot("set xlabel 't'\n", "set ylabel 'y2'\n");
                }
                if (gp) {
                    std::this_thread::sleep_for(std::chrono::nanoseconds(4));
                    w2->show();
                    emit signal1();
                    this->close();
                }
            }
            else if (ui->radioButton_3->isChecked()) {
                write_data(y1, y2);
                gp = call_gnuplot("set xlabel 'y1'\n", "set ylabel 'y2'\n");
                if (gp) {
                    std::this_thread::sleep_for(std::chrono::nanoseconds(4));
                    w3->show();
                    QString name = singular_point_name(1.0, -((1.0 - A) * A / B + A * B / (A - 1.0)), A * (1.0 - A));
                    QString coors = "(" + QString::number(Y1) + ", " + QString::number(Y2) + ")";
                    emit signal(name, coors);
                    this->close();
                }
            }
        }
    }
}

void runge_kutta_method(dv &y1, dv &y2, dv &t) {
    size_t n = t.size();
    double steps[3] = {t[1], t[1] / 2.0, t[1] / 6.0};
    for (size_t i = 1; i < n; i++) {
        calc_next(y1[i - 1], y2[i - 1], steps, y1[i], y2[i]);
        double d11 = y1[i - 1] / y1[i], d12 = y2[i - 1] / y2[i];
        double d21 = fabs(y1[i] - y1[i - 1]), d22 = fabs(y2[i] - y2[i - 1]);
        if (d11 > 10e2 || d12 > 10e2 || (d21 < 10e-10 && d22 < 10e-10)) {
            y1.resize(i); y2.resize(i); t.resize(i);
            break;
        }
    }
}

void calc_next(double y1, double y2, double steps[], double &y1_next, double &y2_next) {
    double k1 = f1(y1, y2), m1 = f2(y1, y2);
    double p1 = y1 + steps[1] * k1, p2 = y2 + steps[1] * m1;
    double k2 = f1(p1, p2), m2 = f2(p1, p2);
    p1 = y1 + steps[1] * k2; p2 = y2 + steps[1] * m2;
    double k3 = f1(p1, p2), m3 = f2(p1, p2);
    p1 = y1 + steps[0] * k3; p2 = y2 + steps[0] * m3;
    double k4 = f1(p1, p2), m4 = f2(p1, p2);
    y1_next = y1 + steps[2] * (k1 + 2 * (k2 + k3) + k4);
    y2_next = y2 + steps[2] * (m1 + 2 * (m2 + m3) + m4);
}

const char* singular_point_name(double a, double b, double c) {
    double d = b * b - 4.0 * a * c;
    if (d > 0.0) {
        d = sqrt(d);
        double root1 = (-b + d) / (2.0 * a), root2 = (-b - d) / (2.0 * a);
        if (root1 < 0.0 && root2 < 0.0) return "Асимптотически устойчивый узел";
        else if (root1 > 0.0 && root2 > 0.0) return "Неустойчивый узел";
        else if (root1 > 0.0 && root2 < 0.0) return "Седло";
    }
    double root = - b / (2.0 * a);
    if (root < 0.0) {
        if (!d) return "Асимптотически устойчивый вырожденный узел";
        return "Асимптотически устойчивый фокус";
    }
    else {
        if (!d) return "Неустойчивый вырожденный узел";
        return "Неустойчивый фокус";
    }
}

void write_data(const dv& v1, const dv& v2) {
    size_t n = v1.size();
    std::ofstream ofs("data.txt");
    for (size_t i = 0; i < n; i++) ofs << v1[i] << ' ' << v2[i] << std::endl;
    ofs.close();
}

FILE* call_gnuplot(const char *str1, const char *str2) {
    const char *Usage = NULL;
    FILE* gp = popen("gnuplot","w");
    fprintf(gp, str1, Usage);
    fprintf(gp, str2, Usage);
    fprintf(gp, "set terminal png size 1200,430\n");
    fprintf(gp, "set output 'graph.png'\n");
    fprintf(gp, "plot 'data.txt' notitle\n");
    pclose(gp);
    return gp;
}


void MainWindow::on_pushButton_2_clicked()
{
    remove("graph.png");
    remove("data.txt");
    QApplication::quit();
}
