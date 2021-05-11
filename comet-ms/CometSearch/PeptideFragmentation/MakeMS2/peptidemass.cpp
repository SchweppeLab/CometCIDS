
#include "peptidemass.h"

#include <map>
#include <vector>

using namespace std;

map<char, double> aaMasses = {
    {'A', 71.03711},
    {'R', 156.10111},
    {'N', 114.04293},
    {'D', 115.02694},
    {'C', 160.030654}, // 57.021464 original: 103.00919,  assume all C are modified Cam
    {'E', 129.04259},
    {'Q', 128.05858},
    {'G', 57.02146},
    {'H', 137.05891},
    {'I', 113.08406},
    {'L', 113.08406},
    {'K', 128.09496},
    {'M', 131.04049},
    {'F', 147.06841},
    {'P', 97.05276},
    {'S', 87.03203},
    {'T', 101.04768},
    {'W', 186.07931},
    {'Y', 163.06333},
    {'V', 99.06841},
    {'U', 150.95363},
};

float masstomz(float mass, int charge) {
    return (mass + (charge-1) * 1.00727646688)/charge;
}

//TODO: Phil's new way - problems? issues?
//float masstomz(float mass, int charge) {
//    return ((mass + charge * 1.00782503224) - charge*0.00054856587)/charge;
//}

vector<double> generateRev(const string peptide) {
    vector<double> output;
    double mass = 1.00733 + 18.01078;
    int peptideLength = peptide.size();
    for (int i = peptideLength - 1; i >= 0; --i) {
        mass += aaMasses[peptide.at(i)];
        output.push_back(mass);
    }
    return output;
}

vector<double> generateFwd(const string peptide) {
    vector<double> output;
    double mass = 1.00733;
    for (int i = 0; i < peptide.size(); ++i) {
        mass += aaMasses[peptide.at(i)];
        output.push_back(mass);
    }
    return output;
}

void generateFragments(const string peptide, int fragmentType, vector<double>& output) {
    auto length = peptide.length();
    output.reserve(output.size() + length);
    vector<double> masses(length);
    for (auto i(0); i < length; ++i) {
        masses[i] = aaMasses.at(peptide.at(i));
    }

    if (fragmentType == MASSES_B) {
        for (auto end(0); end < length - 1; ++end) {
            double mass = 1.00733;
            for (auto i(0); i <= end; ++i) {
                mass += masses[i];
            }
            output.push_back(mass);
        }
    }
    else if(fragmentType == MASSES_Y) {
        for (int start(length - 1); start >= 0; --start) {
            double mass = 1.00733 + 18.01078;
            for (auto i(start); i < length; ++i) {
                mass += masses[i];
            }
            output.push_back(mass);
        }
    }
}
