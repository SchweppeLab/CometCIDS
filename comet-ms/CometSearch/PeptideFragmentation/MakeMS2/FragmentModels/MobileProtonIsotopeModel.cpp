
#include "MobileProtonIsotopeModel.h"
#include "PeptideFragmentation/MakeMS2/peptidemass.h"

using namespace std;

vector< vector<FragmentIon> > MobileProtonIsotopeModel::run(string peptide, const FragmentModelData inputData) {
    vector< vector<FragmentIon> > output(inputData.maxCharge + 1);
    auto peptideInfo = readPeptideInfo(peptide, inputData.obsCharge);
    if (peptideInfo.countH + peptideInfo.countK + peptideInfo.countR + 1 < inputData.obsCharge) {
        return output;
    }

    auto charges = getFragmentCharges(peptide, peptideInfo);

    auto peptideLength = peptide.size();
    auto dists = options.isotopes.makeFragmentDist(peptide, inputData.nHeavy, options.isUseSqrtNormalization, options.isDebug);

    // B Ions.
    // Use length - 1 to avoid full-length peptide.
    for (auto i = 0u; i < peptideLength - 1; ++i) {
        for (int j = 0; j <= inputData.nHeavy; ++j) {
            const auto &chargeRange = charges.b.at(i);
            for (auto k = max(chargeRange.min, 1); k <= chargeRange.max && k <= inputData.maxCharge; ++k) {
                double mz = masstomz(inputData.aaFwdMasses[i] + (j * DEUTERIUM_MASS_DIFF), k);
                if (!options.useMzRange || (mz > inputData.minMz && mz < inputData.maxMz)) {
                    FragmentIon b = { mz, dists.bIons.at(j).at(i) };
                    output[k].push_back(b);
                }
            }
        }
    }

    // Y Ions.
    // Use length - 1 to avoid full-length peptide.
    for (int i = peptideLength - 1; i > 0; --i) {
        for (int j = 0; j <= inputData.nHeavy; ++j) {
            const auto &chargeRange = charges.y.at(i);
            for (auto k = max(chargeRange.min, 1); k <= chargeRange.max && k <= inputData.maxCharge; ++k) {
                double mz = masstomz(inputData.aaRevMasses[peptideLength - (i + 1)] + (j * DEUTERIUM_MASS_DIFF), k);
                if (!options.useMzRange || (mz > inputData.minMz && mz < inputData.maxMz)) {
                    FragmentIon b = { mz, dists.yIons.at(j).at(peptideLength - (i + 1)) };
                    output[k].push_back(b);
                }
            }
        }
    }

    if(options.normalize) {
        this->normalize(output);
    }
    return output;
}
