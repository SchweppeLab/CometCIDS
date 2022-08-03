
#include "FragmentCalculator.h"

// for Comet Data Internal
#include <cstring>
#include <string>
#include "Common.h"
using namespace std;

#include "CometDataInternal.h"

#include <boost/algorithm/string.hpp>

#include <set>

FragmentCalculator::FragmentCalculator() {

}

FragmentCalculator& FragmentCalculator::instance() {
    static FragmentCalculator calculator;
    return calculator;
}

void FragmentCalculator::setOption(string key, string value) {
    boost::trim_left(value);
    boost::trim_right(value);
    if (key == "model") {
        FragmentModelOptions::Model model;
        if (value == "charge_isotope") {
            model = FragmentModelOptions::CHARGE_ISOTOPE;
        } else if (value == "charge") {
            model = FragmentModelOptions::CHARGE;
        } else if (value == "isotope") {
            model = FragmentModelOptions::ISOTOPE;
        } else if (value == "mobile_proton") {
            model = FragmentModelOptions::MOBILE_PROTON;
        } else if (value == "mobile_proton_isotope") {
            model = FragmentModelOptions::MOBILE_PROTON_ISOTOPE;
        } else if (value == "mobile_proton2") {
            model = FragmentModelOptions::MOBILE_PROTON2;
        } else if (value == "mobile_proton_isotope2") {
            model = FragmentModelOptions::MOBILE_PROTON_ISOTOPE2;
        } else if (value == "mobile_proton3") {
            model = FragmentModelOptions::MOBILE_PROTON3;
        } else if (value == "mobile_proton_isotope3") {
            model = FragmentModelOptions::MOBILE_PROTON_ISOTOPE3;
        } else if (value == "hyper_frag") {
            model = FragmentModelOptions::HYPER_FRAG;
        } else if (value == "") {
            model = FragmentModelOptions::PLAIN;
        } else {
            throw "Invalid fragment model type";
        }
        generator.initialize(model);
    }
    else {
        generator.setOption(key, value);
    }
}

std::vector<string> FragmentCalculator::validOptionChars() {
    return {
        "A",
        "R",
        "N",
        "D",
        "C",
        // "Cam",
        "Q",
        "E",
        "G",
        "H",
        "I",
        "L",
        "X",
        "K",
        "M",
        "F",
        "P",
        "S",
        "T",
        "W",
        "Y",
        "V",
        "O",
        "U",
        "*"
    };
}

void FragmentCalculator::printOptions() {
    generator.printOptions();
}

void FragmentCalculator::generateFragments(
        char* peptide,
        int peptideLength,
        FragmentModelData inputData,
        vector<vector<FragmentIon> >& output)
{
    if(peptideLength > 99) {
        return;
    }
    string p(peptide, peptideLength);
    output = generator.getChargeFragments(p, inputData);
}

double FragmentCalculator::getFragmentIonMass(int iWhichIonSeries,
                        int i,
                        int ctCharge,
                        double *_pdAAforward,
                        double *_pdAAreverse) {
    double dFragmentIonMass = 0.0;

    switch (iWhichIonSeries)
    {
    case ION_SERIES_B:
        dFragmentIonMass = _pdAAforward[i];
        break;
    case ION_SERIES_Y:
        dFragmentIonMass = _pdAAreverse[i];
        break;
    case ION_SERIES_A:
        dFragmentIonMass = _pdAAforward[i] - g_staticParams.massUtility.dCO;
        break;
    case ION_SERIES_C:
        dFragmentIonMass = _pdAAforward[i] + g_staticParams.massUtility.dNH3;
        break;
    case ION_SERIES_Z:
        dFragmentIonMass = _pdAAreverse[i] - g_staticParams.massUtility.dNH2;
        break;
    case ION_SERIES_X:
        dFragmentIonMass = _pdAAreverse[i] + g_staticParams.massUtility.dCOminusH2;
        break;
    }

    return masstomz(dFragmentIonMass, ctCharge);
}

double FragmentCalculator::masstomz(double mass, int charge) {
    return (mass + (charge-1)*PROTON_MASS)/charge;
}

double FragmentCalculator::mztomass(double mz, int charge) {
    return (mz * charge) - ((charge - 1) * PROTON_MASS);
}
