
#include <vector>
#include <math.h>
#include <stdexcept>

#include "mathlib.h"

using namespace std;

double nCr(int n, int r) {
    // needs to be type safe
    double retVal = 1;
    for (int i = 1; i <= n; ++i) {
        retVal *= i;
        if (i <= r) {
            retVal /= i;
        }
        if (i <= (n - r)) {
            retVal /= i;
        }
    }

    return retVal;
}

double binomialP(int n, int k, float p) {
    return (nCr(n, k) * pow(p, k) * pow((1.0 - p), (n - k)));
}

vector<float> matrixMultiply(vector<float> x, vector<vector<float> > A) {
    if (x.size() != A.size()) {
        throw std::runtime_error("invalid dimensions for matrix multiply");
    }
    int length = A.size();
    int columns = A.at(0).size();

    vector<float> output(columns);
    for (auto i(0); i < length; ++i) {
        int currentColumns = A.at(i).size();
        for (auto j(0); j < columns && j < currentColumns; ++j) {
            output[j] += A.at(i).at(j) * x.at(i);
        }
    }
    return output;
}

void vectorMultiply(vector<float>& x, const vector<float>& y) {
    int length = x.size();
    for (auto i(0); i < length; ++i) {
        x[i] *= y.at(i);
    }
}

void vectorMultiply(vector<float>& x, float y) {
    int length = x.size();
    for (auto i(0); i < length; ++i) {
        x[i] *= y;
    }
}

void vectorDivide(vector<float>& x, float y) {
    int length = x.size();
    for (auto i(0); i < length; ++i) {
        x[i] /= y;
    }
}

void vectorAdd(vector<float>& x, const vector<float>& y) {
    int length = x.size();
    for (auto i(0); i < length; ++i) {
        x[i] += y.at(i);
    }
}

