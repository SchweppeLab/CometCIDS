
#include "IsotopeModel.h"
#include "PeptideFragmentation/MakeMS2/dist.h"
#include "PeptideFragmentation/MakeMS2/peptidemass.h"
#include "PeptideFragmentation/MakeMS2/SpectrumGenerator.h"

using namespace std;

vector< vector<FragmentIon> > IsotopeModel::run(string peptide, const FragmentModelData inputData) {

    auto dists = options.isotopes.makeFragmentDist(peptide, inputData.nHeavy, options.isUseSqrtNormalization, false);

    auto peptideLength = peptide.size();
    vector< vector<FragmentIon> > output(inputData.maxCharge + 1);

    // Use length - 1 to avoid full-length peptide.
    for (int i = 0; i < peptideLength - 1; ++i) {
        for (int j = 0; j <= inputData.nHeavy; ++j) {
             for (int k = 1; k <= max(min(inputData.maxCharge, inputData.obsCharge-1),1); ++k) {
                 // B ions.
                 double mz = masstomz(inputData.aaFwdMasses[i] + (j * options.isotopeMassDiff), k);
                 if (!options.useMzRange || (mz > inputData.minMz && mz < inputData.maxMz)) {
                     FragmentIon b = { mz, dists.bIons.at(j).at(i) };
                     output[k].push_back(b);
                 }
             }
        }
    }
    for (int i = 0; i < peptideLength - 1; ++i) {
        for (int j = 0; j <= inputData.nHeavy; ++j) {
             for (int k = 1; k <= max(min(inputData.maxCharge, inputData.obsCharge-1),1); ++k) {
                 // Y ions.
                 double mz = masstomz(inputData.aaRevMasses[i] + (j * options.isotopeMassDiff), k);
                 if (!options.useMzRange || (mz > inputData.minMz && mz < inputData.maxMz)) {
                     FragmentIon y = { mz, dists.yIons.at(j).at(i) };
                     output[k].push_back(y);
                 }
             }
        }
    }

    if(options.normalize) {
        this->normalize(output);
    }
    return output;
}

