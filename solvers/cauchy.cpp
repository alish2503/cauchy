#include <iostream>
#include <vector>
#include <cmath>
#include "cauchy.h"

typedef std::vector<double> dv;
typedef std::vector<dv> dm;

int n;
double T = 0.3;
double tau = 1.e-6;
double eps = 1.e-8; 

void explEuler(dv&);
void impl(dv&, int);
void explMidpoint(dv&);
void Newton(dm&, dv&, bool);
void RadauIA(dv&);
void RadauIAfrozen(dv&);
double norm2(const dv&, const dv&);
void F(dv&, const dv&, double);
void dF(dv&, const dv&, double);
void EulerSLAE(dm&, const dv&, const dv&, double);
void MidPointSLAE(dm&, const dv&, const dv&, double);
void TrapezoidalSLAE(dm&, const dv&, const dv&, double);
void Adams_SLAE(dm&, const dv&, const dm&, const dv&, double);
void Gear_SLAE(dm&, const dv&, const dm&, double);
void RadauIA_Jacobimatrix(dm&, const dv&, const dm&, const dv&, double);
void RadauIA_F(dm&, const dv&, const dm&, const dv&, double);
void Gauss(dm&);
void LU(dm&);
void explAdams(dv&);
void implAdams(dv&);
void Gear(dv&);
void solveLU(dm&);

int main(void) {
    dv y = {0.3, 0.0};
    n = y.size();
    explAdams(y);  
}

void explEuler(dv &y) {
    dv f(n);
    double t = tau;
    while(t < T) {
        F(f, y, t-tau);
        for (int i = 0; i < n; i++) y[i] = y[i] + tau * f[i];
        t += tau;
    }
}
void impl(dv &y, int param) {
    double t = tau;
    while(t < T) {
        dv y_next = y;
        do {
            dv y_prev = y_next;
            dm J;
            if (!param) EulerSLAE(J, y_prev, y, t);
            else if(param==1) MidPointSLAE(J, y_prev, y, t);
            else if(param==2) TrapezoidalSLAE(J, y_prev, y, t);
            Newton(J, y_next, false);
            if (norm2(y_prev, y_next) <= eps) break;
        } while(true);
        y = y_next;
        std::cout << y[0] << ' ' << y[1] << std::endl;
        t += tau;
    }
}
void explAdams(dv &y) {
    double t = tau;
    dm f(4, dv(n));
    for (int i = 0; i < 3; i++) {
        F(f[i], y, t-tau);
        for (int j = 0; j < n; j++) y[j] += tau * f[i][j];
        t += tau;
    }
    F(f[3], y, t-tau);
    while (true) {     
        for (int j = 0; j < n; j++) y[j] += tau/24*(55*f[3][j]-59*f[2][j]+37*f[1][j]-9*f[0][j]);
        std::cout << y[0] << ' ' << y[1] << std::endl;
        t += tau;
        if (t >= T) break;
        for (int k = 0; k < 3 ; k++) f[k] = f[k+1];
        F(f[3], y, t - tau);
    }
}
void implAdams(dv &y) {
    double t = tau;
    dm f(4, dv(n));
    for (int i = 0; i < 2; i++) {
        F(f[i], y, t-tau);
        for (int j = 0; j < n; j++) y[j] += tau * f[i][j];
        t += tau;
    }
    F(f[2], y, t-tau);
    while (true) { 
        dv y_next = y;  
        do {
            dv y_prev = y_next;
            F(f[3], y_prev, t);
            dm J;    
            Adams_SLAE(J, y_prev, f, y, t);
            Newton(J, y_next, false);
            if (norm2(y_prev, y_next) <= eps) break;
        } while(true);  
        y = y_next;
        std::cout << y[0] << ' ' << y[1] << std::endl;
        t += tau;
        if (t >= T) break;
        for (int k = 0; k < 3 ; k++) f[k] = f[k+1];
        F(f[3], y, t);
    }
}
void Adams_SLAE(dm &J, const dv &y_prev, const dm &f, const dv& y, double t) {
    dv df(n*n);
    dF(df, y_prev, t);
    J = {
        {1.0 - 9.0/24*tau*df[0], -9.0/24*tau*df[1], tau/24*(9*f[3][0]+19*f[2][0]-5*f[1][0]+f[0][0]) - y_prev[0] + y[0]},
        {-9.0/24*tau*df[2], 1.0 - 9.0/24*tau*df[3], tau/24*(9*f[3][1]+19*f[2][1]-5*f[1][1]+f[0][1]) - y_prev[1] + y[1]}
    };
}
void Gear(dv &y) {
    double t = tau;
    dv f(n);
    dm Y(3, y);
    for (int i = 1; i < 3; i++) {
        F(f, Y[i - 1], t-tau);
        for (int j = 0; j < n; j++) Y[i][j] += tau * f[j];
        t += tau;
    }
    while (true) {   
        dv y_next = Y[2];
        do {  
            dv y_prev = y_next;
            dm J;    
            Gear_SLAE(J, y_prev, Y, t);
            Newton(J, y_next, false);
            if (norm2(y_prev, y_next) <= eps) break;
        } while(true);  
        y = y_next;
        std::cout << y[0] << ' ' << y[1] << std::endl;
        t += tau;
        if (t >= T) break;
        for (int k = 0; k < 2; k++) Y[k] = Y[k+1];
        Y[2] = y;
    }
}
void Gear_SLAE(dm &J, const dv &y_prev, const dm &Y, double t) {
    dv df(n*n), f(n);
    F(f, y_prev, t), dF(df, y_prev, t);
    J = {
        {1.0 - 6.0/11*tau*df[0], -6.0/11*tau*df[1], (18*Y[2][0]-9*Y[1][0]+2*Y[0][0]+6*tau*f[0])/11 - y_prev[0]},
        {-6.0/11*tau*df[2], 1.0 - 6.0/11*tau*df[3], (18*Y[2][1]-9*Y[1][1]+2*Y[0][1]+6*tau*f[1])/11 - y_prev[1]}
    };
}
void RadauIA(dv &y) {
    double t = tau;
    dv k_next = {0.0, 0.0, 0.0, 0.0};
    while (t < T) {
        do {
            dv k_prev = k_next;
            dm J;
            dm args = {
                {y[0] + 0.25*tau*(k_next[0] - k_next[1]), y[1] + 0.25*tau*(k_next[2] - k_next[3])},
                {y[0] + 0.25*tau*(k_next[0] + 5.0/3*k_next[1]), y[1] + 0.25*tau*(k_next[2] + 5.0/3*k_next[3])} 
            };
            RadauIA_Jacobimatrix(J, k_next, args, y, t);
            RadauIA_F(J, k_next, args, y, t);
            Newton(J, k_next, false);
            if (norm2(k_prev, k_next) <= eps) break;
        } while(true);
        y[0] += 0.25*tau*(k_next[0] + 3*k_next[1]);
        y[1] += 0.25*tau*(k_next[2] + 3*k_next[3]);
        t += tau; 
    }           
    std::cout << y[0] << "   " << y[1] << std::endl;  
}
void RadauIAfrozen(dv &y) {
    double t = tau;
    dv k_next = {0.0, 0.0, 0.0, 0.0};
    while (t < T) {
        dm J;
        dm args = {
            {y[0] + 0.25*tau*(k_next[0] - k_next[1]), y[1] + 0.25*tau*(k_next[2] - k_next[3])},
            {y[0] + 0.25*tau*(k_next[0] + 5.0/3*k_next[1]), y[1] + 0.25*tau*(k_next[2] + 5.0/3*k_next[3])} 
        };
        RadauIA_Jacobimatrix(J, k_next, args, y, t - tau);
        LU(J);
        do {
            dv k_prev = k_next;
            RadauIA_F(J, k_next, args, y, t);
            Newton(J, k_next, true);
            if (norm2(k_prev, k_next) <= eps) break;
            args = {
                {y[0] + 0.25*tau*(k_next[0] - k_next[1]), y[1] + 0.25*tau*(k_next[2] - k_next[3])},
                {y[0] + 0.25*tau*(k_next[0] + 5.0/3*k_next[1]), y[1] + 0.25*tau*(k_next[2] + 5.0/3*k_next[3])} 
            };
        } while(true);
        y[0] += 0.25*tau*(k_next[0] + 3*k_next[1]);
        y[1] += 0.25*tau*(k_next[2] + 3*k_next[3]);
        std::cout << y[0] << "   " << y[1] << std::endl;
        t += tau; 
    }           
      
}
void explMidpoint(dv &y) {
    dv f(n);
    double t = tau;
    while(t < T) {
        F(f, y, t-tau);
        y[0] += tau*f1(y[1] + 0.5*tau*f[1]);
        y[1] += tau*f2(y[0] + 0.5*tau*f[0], y[1] + 0.5*tau*f[1]);
        std::cout << y[0] << ' ' << y[1] << std::endl;
        t += tau;
    }
}
void Newton(dm &J, dv &next, bool frozen) {
    int n = J.size();
    if (frozen) solveLU(J);
    else Gauss(J);
    for (int i = 0; i < n; i++) next[i] += J[i].back();
}
void EulerSLAE(dm &J, const dv &y_next, const dv& y, double t) {
    dv f(n), df(n*n);
    F(f, y_next, t); dF(df, y_next, t);
    J = {
        {1.0 - tau*df[0], -tau*df[1], tau*f[0] - y_next[0] + y[0]},
        {-tau*df[2], 1.0 - tau*df[3], tau*f[1] - y_next[1] + y[1]}
    };
}
void MidPointSLAE(dm &J, const dv &y_next, const dv& y, double t) {
    dv f(n), df(n*n);
    dv yMid = {0.5*(y_next[0]+y[0]), 0.5*(y_next[1]+y[1])};
    F(f, yMid, t); dF(df, yMid, t);
    J = {
        {1.0 - 0.5*tau*df[0], -0.5*tau*df[1], tau*f[0] - y_next[0] + y[0]},
        {-0.5*tau*df[2], 1.0 - 0.5*tau*df[3], tau*f[1] - y_next[1] + y[1]}
    };
}
void TrapezoidalSLAE(dm &J, const dv &y_next, const dv& y, double t) {
    dv f(n), df(n*n);
    F(f, y_next, t); dF(df, y_next, t);
    dv sum = f;
    F(f, y, t); 
    sum[0] += f[0]; sum[1] += f[1];
    J = {
        {1.0 - 0.5*tau*df[0], -0.5*tau*df[1], 0.5*tau * sum[0] - y_next[0] + y[0]},
        {-0.5*tau*df[2], 1.0 - 0.5*tau*df[3], 0.5*tau * sum[1] - y_next[1] + y[1]}
    };
}
void RadauIA_Jacobimatrix(dm &J, const dv &k_next, const dm &args, const dv& y, double t) {
    dv df1(n*n), df2(n*n);
    dF(df1, args[0], t);
    dF(df2, args[1], t + 2.0/3*tau);
    J = {
        {1.0 - 0.25*tau*df1[0], 0.25*tau*df1[0], -0.25*tau*df1[1], 0.25*tau*df1[1], 0.0},
        {-0.25*tau*df2[0], 1.0 - 5.0/12*tau*df2[0], -0.25*tau*df2[1], -5.0/12*tau*df2[1], 0.0},
        {-0.25*tau*df1[2], 0.25*tau*df1[2], 1.0 - 0.25*tau*df1[3], 0.25*tau*df1[3], 0.0},
        {-0.25*tau*df2[2], -5.0/12*tau*df2[2], -0.25*tau*df2[3], 1.0 - 5.0/12*tau*df2[3], 0.0}
    };
}
void RadauIA_F(dm &J, const dv &k_next, const dm &args, const dv &y, double t) {
    dv f1(n), f2(n);
    F(f1, args[0], t); 
    F(f2, args[1], t + 2.0/3*tau); 
    J[0].back() = f1[0] - k_next[0];
    J[1].back() = f2[0] - k_next[1];
    J[2].back() = f1[1] - k_next[2];
    J[3].back() = f2[1] - k_next[3];
}
void Gauss(dm &A) {
    int n = A.size();  
    for (int i = 0; i < n - 1; ++i) {
        int maxRowIndex = i;
        for (int k = i + 1; k < n; ++k) {
            if (std::abs(A[k][i]) > std::abs(A[maxRowIndex][i])) {
                maxRowIndex = k;
            }
        }
        std::swap(A[i], A[maxRowIndex]);
        for (int k = i + 1; k < n; ++k) {
            double factor = A[k][i] / A[i][i];
            for (int j = i; j <= n; ++j) {
                A[k][j] -= factor * A[i][j];
            }
        }
    }              
    for (int i = n - 1; i >= 0; --i) {
        for (int j = i + 1; j < n; ++j) {
            A[i][n] -= A[i][j] * A[j][n];
        }
        A[i][n] /= A[i][i];
    }    
}
void LU(dm& A) {
    int n = A.size();
    for (int i = 0; i < n; i++) {
        for (int k = i + 1; k < n; k++) {
            if (abs(A[k][i]) > abs(A[i][i])) {
                swap(A[i], A[k]);
            }
            if (A[i][i] == 0.0) {
                std::cerr << "LU decomposition failed: Matrix is singular." << std::endl;
                return;
            }
            double factor = A[k][i] / A[i][i];
            for (int j = i; j < n; j++) {
                A[k][j] -= factor * A[i][j];
            }
            A[k][i] = factor;
        }
    }
}
void solveLU(dm& A) {
    int n = A.size();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < i; j++) {
            A[i][n] -= A[i][j] * A[j][n];
        }
    }
    for (int i = n - 1; i >= 0; i--) {
        for (int j = i + 1; j < n; j++) {
            A[i][n] -= A[i][j] * A[j][n];
        }
        A[i][n] /= A[i][i];
    }
}
void F(dv &f, const dv &y, double t) {
    f[0] = f1(y[1]);
    f[1] = f2(y[0], y[1]);
}
void dF(dv &df, const dv &y, double t) {
    df[0] = df1dy();
    df[1] = df1dz();
    df[2] = df2dy(y[0]);
    df[3] = df2dz(y[1]);
}

double norm2(const dv &y1, const dv &y2) {
    int n = y1.size();
    double sum = 0.0;
    for (int i = 0; i < n; i++) sum += (y1[i] - y2[i])*(y1[i] - y2[i]);
    return std::sqrt(sum);
}
