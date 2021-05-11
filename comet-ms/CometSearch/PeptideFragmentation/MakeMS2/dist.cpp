
#include "dist.h"

#include "fft.h"
#include "mathlib.h"

#ifdef __APPLE__
#include "eigen3/Eigen/Dense"
#else
#include <Eigen/Dense>
#endif

#include <iostream>
#include <map>
#include <random>
#include <vector>

using namespace std;

map<string, int> n_swaps = {
    {"A", 4},
    {"R", 3},
    {"N", 2},
    {"D", 2},
    {"C", 2},
    {"Cam", 2},
    {"Q", 4},
    {"E", 4},
    {"G", 2},
    {"H", 3},
    {"I", 1},
    {"L", 1},
    {"X", 1},
    {"K", 1},
    {"M", 1},
    {"F", 0},
    {"P", 3},
    {"S", 3},
    {"T", 0},
    {"W", 0},
    {"Y", 0},
    {"V", 1},
    {"O", 0},
    {"U", 0},
    {"*", 0}
};

vector<string> elements = {
    "C",
    "H",
    "N",
    "O",
    "E",
    "S"
};

// Composition using (C, H, N, O, E, S)
map<string, vector<int>> cforms {
   {"A",   { 3,  7,  1, 2, 0, 0} },
   {"R",   { 6,  14, 4, 2, 0, 0} },
   {"N",   { 4,  8,  2, 3, 0, 0} },
   {"D",   { 4,  7,  1, 4, 0, 0} },
   {"C",   { 5,  10, 2, 3, 0, 1} }, // All C treated as Cam
//   {"C",   { 3,  7,  1, 2, 0, 1} },
//   {"Cam", { 5,  10, 2, 3, 0, 1} },
   {"Q",   { 5,  10, 2, 3, 0, 0} },
   {"E",   { 5,  9,  1, 4, 0, 0} },
   {"G",   { 2,  5,  1, 2, 0, 0} },
   {"H",   { 6,  9,  3, 2, 0, 0} },
   {"I",   { 6,  13, 1, 2, 0, 0} },
   {"L",   { 6,  13, 1, 2, 0, 0} },
   {"X",   { 6,  13, 1, 2, 0, 0} },
   {"K",   { 6,  14, 2, 2, 0, 0} },
   {"M",   { 5,  11, 1, 2, 0, 1} },
   {"F",   { 9,  11, 1, 2, 0, 0} },
   {"P",   { 5,  9,  1, 2, 0, 0} },
   {"S",   { 3,  7,  1, 3, 0, 0} },
   {"T",   { 4,  9,  1, 3, 0, 0} },
   {"W",   { 11, 13, 1, 2, 0, 0} },
   {"Y",   { 9,  11, 1, 2, 0, 0} },
   {"V",   { 5,  11, 1, 2, 0, 0} },
   {"O",   { 12, 21, 3, 3, 0, 0} },
   {"U",   { 3,  7,  1, 2, 1, 0} },
   {"*",   { 0,  0,  0, 1, 0, 0} }
};

// isotope abundances of (C, H, N, O, E, S, D)
vector<vector<float> > elementDist = {
    {0.989, 0.011},
    {0.99985, .00015},
    {.9963, .0037},
    {.99762, .00038, .0020},
    {.009, 0, .09, .076, .235, 0, .496, 0, .094},
    {.9502, .0075, .0421, .0002},
    {0, 1, 0, 0},
};

// https://stackoverflow.com/questions/24518989/how-to-perform-1-dimensional-valid-convolution
// Equivalent to convolve(p1, rev(p2), type = "open") in R
template<typename T>
std::vector<T>
conv(std::vector<T> const &f, std::vector<T> const &g) {
  int const nf = f.size();
  int const ng = g.size();
  int const n  = nf + ng - 1;

  //Avoids out-of-range error, but may be hiding a bigger problem
  if (n <= 0) {
      return std::vector<T>(0);
  }

  std::vector<T> out(n, T());
  for(auto i(0); i < n; ++i) {
    int const jmn = (i >= ng - 1)? i - (ng - 1) : 0;
    int const jmx = (i <  nf - 1)? i            : nf - 1;
    for(auto j(jmn); j <= jmx; ++j) {
      out[i] += (f[j] * g[i - j]);
    }
  }
  return out;
}

void removeZeroTail(vector<float>& v) {
    if (!v.size()) {
        return;
    }
    int lastValue = 0;
    for (int i = v.size() - 1; i >= 0; --i) {
        lastValue = i;
        if (v.at(i) > 2e-08) {
            break;
        }
    }
    v.resize(lastValue + 1);
}

/**
 * Calculates the isotope distribution given a vector of elements
 * e.g. ["H", "H", "O"]
 *
 * Useful for calculating the isotope dist. for a single AA or molecule.
 *
 * distAA is a bit of a misnomer since it takes elements as input.
 *
 * @param x
 */
vector<float> distAA(const vector<string>& x) {
    vector<float> out = {1};
    for (string s : x) {
        if (s == "C") {
            out = conv(out, elementDist[0]);
        }
        else if(s == "H") {
            out = conv(out, elementDist[1]);
        }
        else if(s == "O") {
            out = conv(out, elementDist[3]);
        }
        else if(s == "N") {
            out = conv(out, elementDist[2]);
        }
        else if(s == "E") {
            out = conv(out, elementDist[4]);
        }
        else if(s == "S") {
            out = conv(out, elementDist[5]);
        }
        else if(s == "D") {
            out = conv(out, elementDist[6]);
        }
        else {
            throw new std::runtime_error("Unsupported element");
        }
    }
    return out;
}

/**
 * Takes in a composition from above, e.g. { 3,  7,  1, 2, 0, 0}
 * and returns a vector of elements, e.g. ["C", "C", "C", "H", ...]
 * @param composition
 * @return
 */
vector<string> expandElements(const vector<int>& composition) {
    vector<string> output;
    for (auto i(0u); i < composition.size(); ++i) {
        int count = composition.at(i);
        while (count > 0) {
            output.push_back(elements.at(i));
            --count;
        }
    }
    return output;
}

vector<string> swapD(vector<string> in, int n) {
    vector<string> output;
    for (auto s : in) {
        if (s == "H" && n > 0) {
            output.push_back("D");
            --n;
        }
        else {
            output.push_back(s);
        }
    }
    return output;
}

void zapSmall(vector<float>& v) {
    for (auto& x : v) {
        if (x < 0.00000001) {
            x = 0;
        }
    }
}

void removeDist(vector<float>& x, vector<float> y) {
    int xSize = x.size();

    y.resize(xSize, 0);

    auto resultX = fft(x);
    auto resultY = fft(y);

    for (auto i(0); i < xSize; ++i) {
        resultX[i] /= resultY.at(i);
    }
    ifft(resultX, x);
    zapSmall(x);
}

/**
 * Returns the isotopic distribution of an amino acid
 * using the natural abundance isotope distributions.
 */
vector<float> naturalDist(string AA) {
    auto elements = expandElements(cforms[AA]);
    auto v = distAA(elements);
    zapSmall(v);
    removeZeroTail(v);
    return v;
}

/**
 * Inputs: An amino acid and a probability of uptake.
 * Returns the isotopic distribution
 */
vector<float> deuterateDist(string AA, float p) {
    int nLabel = n_swaps[AA];
    auto elements = expandElements(cforms[AA]);

    // generate distributions for 0 .. nLabel D uptake
    vector<vector<float> > distributions;
    for (auto i(0); i <= nLabel; ++i) {
        auto v = distAA(swapD(elements, i));
        zapSmall(v);
        removeZeroTail(v);
        distributions.push_back(v);
    }

    vector<float> binomials;
    for (auto i(0); i <= nLabel; ++i) {
        binomials.push_back(binomialP(nLabel, i, p));
    }

    return matrixMultiply(binomials, distributions);
}

void deconvolute(vector<float>& sumVec, const vector<float>& subVec) {
    Eigen::MatrixXf subMat(sumVec.size(), sumVec.size());
    subMat.setZero(sumVec.size(), sumVec.size());
    for(unsigned int i(0); i < subMat.rows(); ++i) {
        unsigned int j(0);
        for(; j < subVec.size(); ++j) {
            if (j + i >= subMat.rows()) {
                break;
            }
            subMat(j + i, i) = subVec.at(j);
        }
    }

    // Map sumVec to v.  output goes to sum vec
    Eigen::Map<Eigen::VectorXf> v(sumVec.data(), sumVec.size());
    v = subMat.inverse() * v;
}

void removeDistDensity(map<string, vector<float> >& dist, vector<float> remove) {
    for (auto& item : dist) {
        removeDist(item.second, remove);
    }
}

#define FRAGMENT_PEPTIDE 1
#define FRAGMENT_B 2
#define FRAGMENT_Y 3
vector<float> peptideDist(string peptide, const map<string, vector<float> >& AAdensity, int fragmentType) {
    vector<float> peptideDist = {1};
    for (auto c : peptide) {
        string aa(1, c);
        peptideDist = conv(peptideDist, AAdensity.at(aa));
    }

    auto length = peptide.size();

    // composition is counts of C H N O E S
    vector<int> composition(6, 0);
    if (fragmentType == FRAGMENT_PEPTIDE) {
        composition[1] = 2 * (length - 1);
        composition[3] = length - 1;
    }
    else if(fragmentType == FRAGMENT_B) {
        composition[1] = (2 * (length - 1)) + 1;
        composition[3] = length;
    }
    else if(fragmentType == FRAGMENT_Y) {
        composition[1] = (2 * (length - 1)) + 1;
        composition[3] = length - 1;
    }

    auto waterLoss = distAA(expandElements(composition));

    removeDist(peptideDist, waterLoss);

    return peptideDist;
}

ionDistOutput ionDist(
        string peptide,
        const map<string, vector<float> >& AAdensity,
        const vector<vector<float> >& bDeficit,
        const vector<vector<float> >& yDeficit,
        int minSize)
{
    ionDistOutput output;
    vector<float> bIonResult = {1};
    for (auto c : peptide) {
        string aa(1, c);
        bIonResult = conv(bIonResult, AAdensity.at(aa));
        zapSmall(bIonResult);
        removeZeroTail(bIonResult);
        output.bIons.push_back(bIonResult);
    }

    for (int i = 0; i < output.bIons.size(); ++i) {
        output.bIons[i] = conv(output.bIons[i], bDeficit.at(i));
        if (output.bIons[i].size() <= minSize) {
            output.bIons[i].resize(minSize + 1);
        }
    }

    vector<float> yIonResult = {1};
    for (auto it = peptide.rbegin(); it != peptide.rend(); ++it) {
        string aa(1, *it);
        yIonResult = conv(yIonResult, AAdensity.at(aa));
        zapSmall(yIonResult);
        removeZeroTail(yIonResult);
        output.yIons.push_back(yIonResult);
    }

    for (int i = output.yIons.size() - 1; i >= 0; --i) {
        output.yIons[i] = conv(output.yIons[i], yDeficit.at(i));
        if (output.yIons[i].size() <= minSize) {
            output.yIons[i].resize(minSize + 1);
        }
    }

    output.peptide = conv(output.bIons.at(0),
                          output.yIons.at(output.yIons.size() - 2));
    return output;
}

map<string, vector<float> > AANaturalIsotopicDists() {
    map<string, vector<float> > out;
    for (auto& item : cforms) {
        auto dist = naturalDist(item.first);
        zapSmall(dist);
        removeZeroTail(dist);
        out[item.first] = dist;
    }
    return out;
}

map<string, vector<float> > AADeuteratedIsotopicDists(float p) {
    map<string, vector<float> > out;
    for (auto& item : cforms) {
        auto dist = deuterateDist(item.first, p);
        zapSmall(dist);
        removeZeroTail(dist);
        out[item.first] = dist;
    }
    return out;
}

vector<float> readDist(string in) {
    std::stringstream ss(in);
    float x;
    vector<float> output;
    while (ss >> x)
    {
        output.push_back(x);
        if (ss.peek() == ',') {
            ss.ignore();
        }
    }

    return output;
}
