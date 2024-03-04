double f1(double z) {
    return z;
}
double f2(double y, double z) {
    return -0.2*z-10*std::sin(y);
}
double f(double y, double t) {
    return t * y;
}
double df(double t) {
    return t;
}
double df1dy() {
    return 0.0;
}
double df1dz(){
    return 1.0;
}
double df2dy(double y) {
    return -10*std::cos(y);
}
double df2dz(double z) {
    return -0.2;
}