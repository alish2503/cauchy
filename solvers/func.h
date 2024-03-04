#include <math.h>

double u_a1(double x, double y) {return sin(M_PI * x) * cos(M_PI * y);}
double fun1(double x, double y) {return 2 * M_PI * M_PI * sin(M_PI * x) * cos(M_PI * y);}
double u_a2(double x, double y) {return exp(2 * x) * sin(2 * y);}
double fun2(double x, double y) {return 0;}
double u_a3(double x, double y) {return y * pow(x, 3) + x * y * y;}
double fun3(double x, double y) {return -2 * x * (3 * y + 1);}
double u_a4(double x, double y) {return x * x + y * y;}
double fun4(double x, double y) {return -4;}
double north(double x) {return u_a1(x, 1.0);}
double south(double x) {return u_a1(x, 0.0);}
double east(double y) {return u_a1(1.0, y);}
double west(double y) {return u_a1(0.0, y);}
