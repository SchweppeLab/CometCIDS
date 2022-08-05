#include "HyperFragModel.h"
#include "PeptideFragmentation/MakeMS2/peptidemass.h"

#include "../../../../../extern/Rmath/dhyper.c"
#include "../../../../../extern/Rmath/dbinom.c"
#include "../../../../../extern/Rmath/stirlerr.c"
#include "../../../../../extern/Rmath/bd0.c"

using namespace std;

int R_finite(double x){
    return isfinite(x);
}

vector<double> HyperFragModel::hyperFrag(int ion, int peptideLength, int numDeuteria){

    vector<double> intensities(static_cast<unsigned long>(numDeuteria+1));

    for (unsigned int i = 0; i < intensities.size(); i++) {
        auto intensity = dhyper(i, numDeuteria, (peptideLength-numDeuteria), ion, false);
        if (isfinite(intensity)) {
            intensities[i] = intensity;
        }
    }

    return intensities;
}

vector< vector<FragmentIon> > HyperFragModel::run(string peptide, const FragmentModelData inputData) {

    auto peptideLength = peptide.size();
    vector< vector<FragmentIon> > output(inputData.maxCharge + 1);

    // Use length - 1 to avoid full-length peptide.
    for (int i = 0; i < peptideLength - 1; ++i) {

        vector<double> intensities = hyperFrag(i, static_cast<int>(peptideLength), inputData.nHeavy);

        for (int j = 0; j <= inputData.nHeavy; ++j) {
             for (int k = 1; k <= max(min(inputData.maxCharge, inputData.obsCharge-1),1); ++k) {
                 // B ions.
                 double mz = masstomz(inputData.aaFwdMasses[i] + (j * DEUTERIUM_MASS_DIFF), k);
                 if (!options.useMzRange || (mz > inputData.minMz && mz < inputData.maxMz)) {
                     FragmentIon b = { mz, intensities.at(j) };
                     output[k].push_back(b);
                 }
             }
        }
    }

    for (int i = 0; i < peptideLength - 1; ++i) {

        vector<double> intensities = hyperFrag(i, static_cast<int>(peptideLength), inputData.nHeavy);

        for (int j = 0; j <= inputData.nHeavy; ++j) {
             for (int k = 1; k <= max(min(inputData.maxCharge, inputData.obsCharge-1),1); ++k) {
                 // Y ions.
                 double mz = masstomz(inputData.aaRevMasses[i] + (j * DEUTERIUM_MASS_DIFF), k);
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
