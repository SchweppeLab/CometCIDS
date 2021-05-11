#ifndef FFT_H
#define FFT_H

#include <complex>
#include <valarray>
#include <vector>

typedef std::complex<float> Complex;
typedef std::valarray<Complex> CArray;

std::vector<std::complex<float> > fft(std::vector<float> &r);

void ifft(std::vector<std::complex<float> > &in, std::vector<float>& out);

#endif // FFT_H
