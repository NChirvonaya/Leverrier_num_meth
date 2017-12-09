#define _CRT_SECURE_NO_WARNINGS

#include "Leverrier.h"

using namespace std;

int main()
{
  freopen("input.txt", "r", stdin);

  int n;
  cin >> n;
  vector<vector<double> > A(n, vector<double>(n, 0.0));
  vector<double> sum(n, 0.0);

  for (int i(0); i < n; ++i) {
    for (int j(0); j < n; ++j) {
      cin >> A[i][j];
      sum[i] += abs(A[i][j]);
    }
  }

  vector<double> pcoef;
  countPolynomCoefs(A, pcoef);
  
  double c(0.0);
  double rad(0.0);
  double l(1000.0), r(-1000.0);
  vector<pair<double, double> > br;
  for (size_t i(0); i < A.size(); ++i) {
    rad = -abs(A[i][i]) + sum[i];
    c = A[i][i];
    l = min(l, c - rad);
    r = max(r, c + rad);
  }

  vector<double> eigenvalues;
  double dx(1e-3);
  int m(1.0 * (r - l) / dx);
  double prev(0.0), cur(0.0);
  for (int i(0); i <= m; ++i) {
    double t(l + i * dx);
    cur = f(t, pcoef);
    if (fabs(cur) < eps) {
      br.push_back({ l + (i - 1) * dx, l + i * dx + eps });
      prev *= (-1);
    }
    else {
      if (prev * cur < 0) {
        br.push_back({ l + (i - 1) * dx, l + i * dx });
      }
      prev = cur;
    }
  }

  for (size_t i(0); i < br.size(); ++i) {
    eigenvalues.push_back(chordMeth(br[i].first, br[i].second, pcoef));
  }

  for (size_t i(0); i < eigenvalues.size(); ++i) {
    cout << eigenvalues[i] << endl;
  }

  return 0;
}