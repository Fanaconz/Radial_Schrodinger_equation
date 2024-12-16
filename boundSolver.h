#ifndef BOUNDSOLVER_H
#define BOUNDSOLVER_H

#include <iostream>
#include <functional>
#include <armadillo>

using namespace std;
using namespace arma;

class boundSolver {
public:
    boundSolver( \
            std::function< double(double) >& v, const int l, \
            const double R);
    // l - orbital moment
    // R - the point of setting the asymptotic boundary condition
    // v - potential

    double F(double);
    mat mat_A(int);
    mat mat_B(int);
    void assert_real(const cx_mat&);

    void printE() const;

    // displays the energies of the associated states
    void writeWF(const string ener) const;
    virtual ~boundSolver();

private:
    int n = 1000;
    std::function < double(double) > v;
    int l;
    double R;
    double h;
    cx_vec eval;
    cx_mat evec;
};

#endif /* BOUNDSOLVER_H */
