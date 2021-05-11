
#include "fft.h"

#include <iostream>
#include <eigen3/unsupported/Eigen/FFT>

using namespace std;

const double PI = 3.141592653589793238460;

vector<complex<float> > fft(vector<float>& real) {
    Eigen::FFT<float> eigenFFT;
    vector<complex<float> > result(real.size());
    eigenFFT.fwd(result, real);
    return result;
}

void ifft(vector<complex<float> >& in, vector<float>& out) {
    Eigen::FFT<float> eigenFFT;
    eigenFFT.inv(out, in);
}
