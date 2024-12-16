#include "boundSolver.h"

#include <cmath>
#include <fstream>

boundSolver::boundSolver(std::function<double(double)>& v, const int l,
                         const double R) {
    this->v = v;
    this->l = l;
    this->R = R;
    h = R / n;
    mat A = mat_A(n);
    mat B = mat_B(n);

    eig_gen(eval, evec, B);
    assert_real(evec);

    mat V(n, n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            V(i, j) = evec(i, j).real();
        }
    }
    mat D = V.t() * B * V;

    vec v0(n);
    for (int i = 0; i < n; i++) {
        v0(i) = 1 / D(i, i);
    }

    mat D_reversed = diagmat(v0);

    mat C = D_reversed * V.t() * A * V;

    eig_gen(eval, evec, C);
    assert_real(evec);
    
}

// f(r,l) =  l(l+1)/r^2 + v(r)
double boundSolver::F(double r) {
    return v(r) + l*(l+1) / r / r;
}

// E = h^2 /12 f(r+i) + 10/12 h^2 f(r) + h^2 /12 f(r-i)
mat boundSolver::mat_A(int n) {
    mat V(n, 3);
    for (int i = 1; i < n + 1; i++) {
        V(i - 1, 0) = -1 + h*h/12 * F(i*h);
        V(i - 1, 1) = 2 + h*h*10/12 * F(i*h);
        V(i - 1, 2) = -1 + h*h/12 * F(i*h);
    }
    ivec temp = {-1, 0, +1};
    mat A = diags(V, temp, n, n);
    return A;
}

// 3diag matr with diags 1 10 1
mat boundSolver::mat_B(int n) {
    mat V(n, 3);
    for (int i = 0; i < n; i++) {
        V(i, 0) = h*h/12;
        V(i, 1) = h*h*10/12;
        V(i, 2) = h*h/12;
    }
    ivec temp = {-1, 0, +1};
    mat B = diags(V, temp, n, n);
    return B;
}

void boundSolver::assert_real(const cx_mat& mat) {
    for (int i = 0; i < mat.n_rows; i++) {
        for (int j = 0; j < mat.n_cols; j++) {
            if (mat(i, j).imag() != 0) {
                cout << i << ' ' << j << ' ' << mat(i, j).imag() << endl;
                throw 1;
            }
        }
    }
}

void boundSolver::printE() const {
    for (auto elem : eval) {
        if (elem.real() < 0) {
            cout << elem.real() << endl;
        }
    }
}


void boundSolver::writeWF(const string ener) const {
    ofstream myfile;
    myfile.open(ener);
    evec.print(myfile);
    myfile.close();
}

boundSolver::~boundSolver() {}