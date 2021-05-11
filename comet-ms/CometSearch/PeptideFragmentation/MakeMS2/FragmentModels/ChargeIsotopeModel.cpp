#include "ChargeIsotopeModel.h"
#include "PeptideFragmentation/MakeMS2/dist.h"
#include "PeptideFragmentation/MakeMS2/hyper.h"
#include "PeptideFragmentation/MakeMS2/peptidemass.h"
#include "PeptideFragmentation/MakeMS2/SpectrumGenerator.h"

using namespace std;

vector< vector<FragmentIon> >
ChargeIsotopeModel::run(string peptide, const FragmentModelData inputData) {
    auto peptideLength = peptide.size();
    auto dists = options.isotopes.makeFragmentDist(peptide, inputData.nHeavy, options.isUseSqrtNormalization, options.isDebug);

    vector<int> bCharges(peptideLength);
    int charges = 1;
    for(int i = 0; i < peptideLength; ++i) {
        const char c = peptide.at(i);
        if (c == 'H' || c == 'K' || c == 'R') {
            charges += 1;
        }
        bCharges[i] = charges;
    }

    vector<int> yCharges(peptideLength);
    charges = 0;
    for(int i = peptideLength - 1; i >= 0; --i) {
        const char c = peptide.at(i);
        if (c == 'H' || c == 'K' || c == 'R') {
            charges += 1;
        }
        yCharges[i] = charges;
    }

    vector< vector<FragmentIon> > output(inputData.maxCharge + 1);
    int appliedObsCharge = inputData.obsCharge + options.obsChargeOffset;

    // Use length - 1 to avoid full-length peptide.
    for (int i = 0; i < peptideLength - 1; ++i) {
        for (int j = 0; j <= inputData.nHeavy; ++j) {
             for (int k = 1; k <= inputData.maxCharge && k <= appliedObsCharge; ++k) {
                 // B ions.
                 auto fraction = dhyper(k, bCharges[i], yCharges[i + 1], appliedObsCharge);
                 if (fraction > options.intensityCutoff) {
                     double mz = masstomz(inputData.aaFwdMasses[i] + (j * DEUTERIUM_MASS_DIFF), k);
                     if (!options.useMzRange || (mz > inputData.minMz && mz < inputData.maxMz)) {
                         FragmentIon b = { mz, dists.bIons.at(j).at(i) };
                         output[k].push_back(b);
                     }
                 }
             }
        }
    }

    for (int i = 0; i < peptideLength - 1; ++i) {
        for (int j = 0; j <= inputData.nHeavy; ++j) {
             for (int k = 1; k <= inputData.maxCharge && k <= appliedObsCharge; ++k) {
                 // Y ions.
                 auto fraction = dhyper(k, yCharges[peptideLength - (1 + i)], bCharges[peptideLength - (2 + i)], appliedObsCharge);
                 if (fraction > options.intensityCutoff) {
                     double mz = masstomz(inputData.aaRevMasses[i] + (j * DEUTERIUM_MASS_DIFF), k);
                     if (!options.useMzRange || (mz > inputData.minMz && mz < inputData.maxMz)) {
                         FragmentIon y = { mz, dists.yIons.at(j).at(i) };
                         output[k].push_back(y);
                     }
                 }
             }
        }
    }

    if(options.normalize) {
        this->normalize(output);
    }
    return output;
}

