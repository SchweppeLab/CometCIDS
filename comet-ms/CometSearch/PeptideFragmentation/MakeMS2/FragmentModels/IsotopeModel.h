#ifndef ISOTOPE_MODEL_H
#define ISOTOPE_MODEL_H

#include <memory>

#include "FragmentModel.h"

class IsotopeModel : public FragmentModel {
public:
   IsotopeModel(FragmentModelOptions& options) : FragmentModel(options) {}
   std::vector< std::vector<FragmentIon> > run(std::string peptide, const FragmentModelData inputData);
};

#endif // ISOTOPE_MODEL_H
