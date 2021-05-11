#ifndef CHARGE_MODEL_H
#define CHARGE_MODEL_H

#include <memory>

#include "FragmentModel.h"

class ChargeModel : public FragmentModel {
public:
   ChargeModel(FragmentModelOptions& options) : FragmentModel(options) {}
   std::vector< std::vector<FragmentIon> > run(std::string peptide, const FragmentModelData inputData);
};

#endif // CHARGE_MODEL_H
