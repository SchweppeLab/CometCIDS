#include "ChargeConditionalFragmentModel.h"

using namespace std;

std::vector< std::vector<FragmentIon> > ChargeConditionalFragmentModel::run(std::string peptide, const FragmentModelData inputData){
    if (inputData.obsCharge <= chargeCutoff) {
        return lowChargeModel->run(peptide, inputData);
    } else {
        return highChargeModel->run(peptide, inputData);
    }
}
