
#include "SpectrumGenerator.h"

#include "mathlib.h"
#include "peptidemass.h"

#include <algorithm>
#include <iostream>

using namespace std;

void SpectrumGenerator::initialize(FragmentModelOptions::Model model) {
    options.isotopes.initialize();
    options.type = model;
    fragmenter = FragmentModelFactory::get(options);
}

void SpectrumGenerator::initialize() {
    initialize(FragmentModelOptions::ISOTOPE);
}

static inline void ltrim(string &s) {
    s.erase(s.begin(), find_if(s.begin(), s.end(), [](int ch) {
        return !isspace(ch);
    }));
}

static inline void rtrim(string &s) {
    s.erase(find_if(s.rbegin(), s.rend(), [](int ch) {
        return !isspace(ch);
    }).base(), s.end());
}

static inline void trim(string &s) {
    ltrim(s);
    rtrim(s);
}

bool readBoolValue(const string& value) {
    if (!value.size()) {
        return false;
    }
    return (bool) stoi(value);
}

void SpectrumGenerator::setOption(string key, string value) {
    trim(value);
    if (key == "relative_minimum_intensity" && value.size()) {
        options.intensityCutoff = stof(value);
    }
    else if(key == "obs_charge_offset" && value.size()) {
        options.obsChargeOffset = stoi(value);
    }
    else if(key == "normalization") {
        options.normalize = readBoolValue(value);
    }
    else if(key == "model_sqrt") {
        options.doSqrt = readBoolValue(value);
    }
    else if(key == "model_windowed_scaling") {
        options.doWindowedScaling = readBoolValue(value);
    }
    else if (key == "isotope_sqrt") {
        options.isUseSqrtNormalization = readBoolValue(value);
    }
    else if (key == "use_mz_range") {
        options.useMzRange = readBoolValue(value);
    }
    else {
        options.isotopes.setOption(key, value);
    }
}

void SpectrumGenerator::printOptions() {
    cout << "=======================================================" << endl;
    options.isotopes.printOptions();
    cout << endl;
    cout << "relative_minimum_intensity = " << options.intensityCutoff << endl;
    cout << "obs_charge_offset = " << options.obsChargeOffset << endl;
    cout << "normalization = " << options.normalize << endl;
    cout << "model_sqrt = " << options.doSqrt << endl;
    cout << "model_windowed_scaling = " << options.doWindowedScaling << endl;
    cout << "isotope_sqrt = " << options.isUseSqrtNormalization << endl;
    cout << "use_mz_range = " << options.useMzRange << endl;
    cout << "=======================================================" << endl;
}

vector< vector<FragmentIon> > SpectrumGenerator::getChargeFragments(
        string peptide,
        FragmentModelData inputData) {
    return fragmenter->run(peptide, inputData);
}

