
#include "MobileProtonModel.h"
#include "PeptideFragmentation/MakeMS2/peptidemass.h"
#include "PeptideFragmentation/MakeMS2/SpectrumGenerator.h"

#include <iostream>

using namespace std;

MobileProtonModel::PeptideInfo MobileProtonModel::readPeptideInfo(string peptide, int charge) {
    PeptideInfo output;

    output.countR = 0;
    output.countK = 0;
    output.countH = 0;
    for (const auto &c : peptide) {
        if (c == 'R') {
            ++output.countR;
        }
        else if (c == 'K') {
            ++output.countK;
        }
        else if (c == 'H') {
            ++output.countH;
        }
    }

    if (output.countR + output.countK + output.countH + 1 < charge) {
        return output;
    }

    if (charge < output.countR) {
        output.mobileAA = 'R';
        output.sequesteredAA = "R";
    }
    else if(charge <= output.countR + output.countK) {
        output.mobileAA = 'K';
        output.sequesteredAA = "R";
    }
    else {
        output.mobileAA = 'H';
        output.sequesteredAA = "RK";
    }

    output.countMobile = 0;
    for (const auto &c : peptide) {
        if (c == output.mobileAA) {
            ++output.countMobile;
        }
    }
    return output;
}

MobileProtonModel::FragmentCharges MobileProtonModel::getFragmentCharges(string peptide, PeptideInfo peptideInfo) {

    auto peptideLength = peptide.size();

    FragmentCharges output;
    output.b.resize(peptideLength);
    output.y.resize(peptideLength);

    // Figure out charge ranges for the b-ions.
    bool firstTarget = true;
    int numMobile = 0;
    Range range;

    // Handle charge at n-term
    if(!(peptideInfo.mobileAA == 'H' && peptideInfo.countMobile)) {
        // n-term may or may not hold a proton if mobile AA == H
        // otherwise, this will count as another target.
        firstTarget = false;
    }
    range.max = 1;

    // Skip peptideLength - 1 since its the full-length peptide.
    for (auto i = 0u; i < peptideLength - 1; ++i) {
        if (peptide[i] == peptideInfo.mobileAA) {
            // Mobile proton from this AA may or may not add to the charge.

            // if mobileAA == H, then were already
            // counting one from the n-term
            if (peptideInfo.mobileAA == 'H' || numMobile > 0) {
                ++range.min;
            }
            if (numMobile < peptideInfo.countMobile - 1) {
                ++range.max;
            }
            ++numMobile;
        }
        else if(peptideInfo.sequesteredAA.find(peptide[i]) != string::npos) {
            // Proton is sequestered at this position,
            // so it will always add to the charge.
            ++range.min;
            ++range.max;
        }
        else if(firstTarget) {
            // A single charge could be anywhere on these AAs,
            // so increase the range by 1 as soon as we start seeing them.
            ++range.max;
            firstTarget = false;
        }

        output.b[i] = range;
    }

    // y-ions
    firstTarget = true;
    numMobile = 0;
    range.min = 0;
    range.max = 0;
    // Skip zero since it represents the full length peptide.
    for (auto i = peptideLength - 1; i > 0; --i) {
        if (peptide[i] == peptideInfo.mobileAA) {
            // Mobile proton from this AA may or may not add to the charge.
            if (numMobile > 0) {
                ++range.min;
            }
            // Also considering n-term for 'H'
            if (peptideInfo.mobileAA == 'H' || numMobile < peptideInfo.countMobile - 1) {
                ++range.max;
            }
            ++numMobile;
        }
        else if(peptideInfo.sequesteredAA.find(peptide[i]) != string::npos) {
            // Proton is sequestered at this position,
            // so it will always add to the charge.
            ++range.min;
            ++range.max;
        }
        else if(firstTarget) {
            // A single charge could be anywhere on these AAs,
            // so increase the range by 1 as soon as we start seeing them.
            ++range.max;
            firstTarget = false;
        }

        output.y[i] = range;
    }

    return output;
}

vector< vector<FragmentIon> > MobileProtonModel::run(string peptide, const FragmentModelData inputData) {
    auto peptideLength = peptide.size();
    vector< vector<FragmentIon> > output(inputData.maxCharge + 1);

    auto peptideInfo = readPeptideInfo(peptide, inputData.obsCharge);
    if (peptideInfo.countH + peptideInfo.countK + peptideInfo.countR + 1 < inputData.obsCharge) {
        return output;
    }

    auto charges = getFragmentCharges(peptide, peptideInfo);

    // B Ions.
    // Use length - 1 to avoid full-length peptide.
    for (auto i = 0u; i < peptideLength - 1; ++i) {
        const auto &chargeRange = charges.b.at(i);
        for (auto k = max(chargeRange.min, 1); k <= chargeRange.max && k <= inputData.maxCharge; ++k) {
            double mz = masstomz(inputData.aaFwdMasses[i], k);
            if (!options.useMzRange || (mz > inputData.minMz && mz < inputData.maxMz)) {
                FragmentIon b = { mz, 1 };
                output[k].push_back(b);
            }
        }
    }

    // Y Ions.
    for (int i = peptideLength - 1; i > 0; --i) {
        const auto &chargeRange = charges.y.at(i);
        for (auto k = max(chargeRange.min, 1); k <= chargeRange.max && k <= inputData.maxCharge; ++k) {
            double mz = masstomz(inputData.aaRevMasses[peptideLength - (i + 1)], k);
            if (!options.useMzRange || (mz > inputData.minMz && mz < inputData.maxMz)) {
                FragmentIon b = { mz, 1 };
                output[k].push_back(b);
            }
        }
    }

    if(options.normalize) {
        this->normalize(output);
    }
    return output;
}
