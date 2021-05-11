#ifndef CHARGECONDITIONALFRAGMENTMODEL_H
#define CHARGECONDITIONALFRAGMENTMODEL_H

#include <memory>
#include "FragmentModel.h"

class ChargeConditionalFragmentModel : public FragmentModel {
public:
       explicit ChargeConditionalFragmentModel(FragmentModelOptions& options,
                                   std::unique_ptr<FragmentModel> lowChargeModelPtr,
                                   std::unique_ptr<FragmentModel> highChargeModelPtr,
                                   int chargeCutOffVal) : FragmentModel(options) {

        lowChargeModel = std::move(lowChargeModelPtr);
        highChargeModel = std::move(highChargeModelPtr);
        chargeCutoff = chargeCutOffVal;
    }

    std::vector< std::vector<FragmentIon> > run(std::string peptide, const FragmentModelData inputData);

private:
    int chargeCutoff;
    std::unique_ptr<FragmentModel> lowChargeModel;
    std::unique_ptr<FragmentModel> highChargeModel;
};

#endif // CHARGECONDITIONALFRAGMENTMODEL_H
