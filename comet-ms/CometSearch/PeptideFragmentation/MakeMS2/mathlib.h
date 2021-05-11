#ifndef MATHLIB_H
#define MATHLIB_H

#include <vector>

double nCr(int n, int r);
double binomialP(int n, int k, float p);



std::vector<float> matrixMultiply(std::vector<float> x, std::vector<std::vector<float> > A);

void vectorMultiply(std::vector<float>& x, const std::vector<float>& y);
void vectorMultiply(std::vector<float>& x, float y);
void vectorDivide(std::vector<float>& x, float y);
void vectorAdd(std::vector<float>& x, const std::vector<float>& y);

#endif // MATHLIB_H
