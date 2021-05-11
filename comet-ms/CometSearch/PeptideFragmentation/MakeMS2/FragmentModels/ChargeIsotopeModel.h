#ifndef CHARGE_ISOTOPE_MODEL_H
#define CHARGE_ISOTOPE_MODEL_H

#include <memory>

#include "FragmentModel.h"

class ChargeIsotopeModel : public FragmentModel {
public:
   ChargeIsotopeModel(FragmentModelOptions& options) : FragmentModel(options) {}
   std::vector< std::vector<FragmentIon> > run(std::string peptide, const FragmentModelData inputData);
};

#endif // CHARGE_ISOTOPE_MODEL_H
