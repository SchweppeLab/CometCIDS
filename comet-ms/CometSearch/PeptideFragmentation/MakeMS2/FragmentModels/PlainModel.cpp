
#include "PlainModel.h"
#include "PeptideFragmentation/MakeMS2/peptidemass.h"

using namespace std;

vector< vector<FragmentIon> > PlainModel::run(string peptide, const FragmentModelData inputData) {
    auto peptideLength = peptide.size();
    vector< vector<FragmentIon> > output(inputData.maxCharge + 1);

    // Use length - 1 to avoid full-length peptide.
    for (int i = 0; i < peptideLength - 1; ++i) {
        for (int k = 1; k <= max(min(inputData.maxCharge, inputData.obsCharge-1),1); ++k) {
            // B ions.
            double mz = masstomz(inputData.aaFwdMasses[i], k);
            if (!options.useMzRange || (mz > inputData.minMz && mz < inputData.maxMz)) {
                FragmentIon b = { mz, 1 };
                output[k].push_back(b);
            }
        }
    }
    for (int i = 0; i < peptideLength - 1; ++i) {
        for (int k = 1; k <= max(min(inputData.maxCharge, inputData.obsCharge-1),1); ++k) {
            // Y ions.
            double mz = masstomz(inputData.aaRevMasses[i], k);
            if (!options.useMzRange || (mz > inputData.minMz && mz < inputData.maxMz)) {
                FragmentIon y = { mz, 1 };
                output[k].push_back(y);
            }
        }
    }
    return output;
}

