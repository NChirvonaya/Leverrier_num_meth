#include "Leverrier.h"

double trace(const std::vector<std::vector<double> > &a) 
{
  double res(0.0);
  for (size_t i(0); i < a.size(); ++i) {
    res += a[i][i];
  }
  return res;
}

std::vector<std::vector<double> > operator*(double p, std::vector<std::vector<double> > &a)
{
  std::vector<std::vector<double> >res(a);
  for (size_t i(0); i < res.size(); ++i) {
    for (size_t j(0); j < res[i].size(); ++j) {
      res[i][j] *= p;
    }
  }
  return res;
}

std::vector<std::vector<double> > operator+(std::vector<std::vector<double> > &a, std::vector<std::vector<double> > &b)
{
  try {
    throw a.size() != b.size() || a[0].size() != b[0].size();
  }
  catch (bool f) {
    if (f) {
      std::cout << "Wrong matrixes size!" << std::endl;
    }
  }
  std::vector<std::vector<double> > res(a);
  for (size_t i(0); i < res.size(); ++i) {
    for (size_t j(0); j < res[i].size(); ++j) {
      res[i][j] += b[i][j];
    }
  }
  return res;
}

std::vector<std::vector<double> > operator-(std::vector<std::vector<double> > &a, std::vector<std::vector<double> > &b)
{
  return a + (-1.0)*b;
}

std::vector<std::vector<double> > operator*(std::vector<std::vector<double> > &a, std::vector<std::vector<double> > &b)
{
  try {
    throw a[0].size() != b.size();
  }
  catch (bool f) {
    if (f)
      std::cout << "Wrong matrixes size!" << std::endl;
  }

  std::vector<std::vector<double> > res(a.size(), std::vector<double>(b[0].size()));

  for (size_t i(0); i < res.size(); ++i) {
    for (size_t j(0); j < res[i].size(); ++j) {
      double el(0.0);
      for (size_t k(0); k < b.size(); ++k) {
        el += a[i][k] * b[k][j];
      }
      res[i][j] = el;
    }
  }
  return res;
}

double f(double x, const std::vector<double> &p)
{
  double res(0.0);
  double d(1.0);
  for (size_t i(0); i < p.size(); ++i) {
    res += p[i] * d;
    d *= x;
  }

  return res;
}

double chordMeth(double a, double b, const std::vector<double> &p)
{

  if (fabs(f(a, p)) < eps) {
    return a;
  }
  double x0, x;

  if (f(a, p) > f(b, p)) {
    x0 = b;
    x = x0 - (x0 - a) * f(x0, p) / (f(x0, p) - f(a, p));
    while (fabs(x - x0) > eps) {
      x0 = x;
      x = x0 - (x0 - a) * f(x0, p) / (f(x0, p) - f(a, p));
    }
  }
  else {
    x0 = a;
    x = x0 - (b - x0) * f(x0, p) / (f(b, p) - f(x0, p));
    while (fabs(x - x0) > eps) {
      x0 = x;
      x = x0 - (b - x0) * f(x0, p) / (f(b, p) - f(x0, p));
    }
  }

  return x;
}

void countPolynomCoefs(std::vector<std::vector<double>>& A, std::vector<double>& pcoef)
{
  int n(A.size());
  std::vector<std::vector<double> > I(n, std::vector<double>(n, 0.0));
  for (size_t i(0); i < I.size(); ++i) {
    I[i][i] = 1.0;
  }
  std::vector<double> p;
  std::vector <std::vector<double> > Bi(A);
  p.push_back(trace(Bi));
  int i(1);
  while (i < n) {
    i++;
    Bi = A * (Bi - p.back() * I);
    p.push_back(trace(Bi) / (1.0 * i));
  }

  for (int i(p.size() - 1); i >= 0; i--) {
    pcoef.push_back(-p[i]);
  }
  pcoef.push_back(1);
}
