#include "HyperFragModel.h"
#include "PeptideFragmentation/MakeMS2/peptidemass.h"

using namespace std;

vector< vector<FragmentIon> > HyperFragModel::run(string peptide, const FragmentModelData inputData) {

    auto peptideLength = peptide.size();
    vector< vector<FragmentIon> > output(inputData.maxCharge + 1);

    HyperFragLookupTable& lookupTable = HyperFragLookupTable::instance();

    vector<vector<double>> intensityDists = lookupTable.getIntensities(inputData.nHeavy, peptideLength);

    // Use length - 1 to avoid full-length peptide.
    for (unsigned int i = 0; i < peptideLength - 1; ++i) {

        vector<double> intensities = intensityDists.at(i);

        for (int j = 0; j <= inputData.nHeavy; ++j) {
             for (int k = 1; k <= max(min(inputData.maxCharge, inputData.obsCharge-1),1); ++k) {
                 // B ions.
                 double mz = masstomz(inputData.aaFwdMasses[i] + (j * options.isotopeMassDiff), k);
                 if (!options.useMzRange || (mz > inputData.minMz && mz < inputData.maxMz)) {
                     FragmentIon b = { mz, intensities.at(j) };
                     output[k].push_back(b);
                 }
             }
        }
    }

    for (unsigned int i = 0; i < peptideLength - 1; ++i) {

        vector<double> intensities = intensityDists.at(i);

        for (int j = 0; j <= inputData.nHeavy; ++j) {
             for (int k = 1; k <= max(min(inputData.maxCharge, inputData.obsCharge-1),1); ++k) {
                 // Y ions.
                 double mz = masstomz(inputData.aaRevMasses[i] + (j * options.isotopeMassDiff), k);
                 if (!options.useMzRange || (mz > inputData.minMz && mz < inputData.maxMz)) {
                     FragmentIon y = { mz, intensities.at(j) };
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
