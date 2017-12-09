#define HG_LEVERRIER_CPP
#ifdef HG_LEVERRIER_CPP

#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>

const double eps = 1e-6;

double trace(const std::vector<std::vector<double> > &a);

std::vector<std::vector<double> > operator*(double p, std::vector<std::vector<double> > &a);

std::vector<std::vector<double> > operator+(std::vector<std::vector<double> > &a, std::vector<std::vector<double> > &b);

std::vector<std::vector<double> > operator-(std::vector<std::vector<double> > &a, std::vector<std::vector<double> > &b);

std::vector<std::vector<double> > operator*(std::vector<std::vector<double> > &a, std::vector<std::vector<double> > &b);

double f(double x, const std::vector<double> &p);

double chordMeth(double a, double b, const std::vector<double> &p);

void countPolynomCoefs(std::vector<std::vector<double> > &A, std::vector<double> &pcoef);

#endif // HG_LEVERRIER_CPP
