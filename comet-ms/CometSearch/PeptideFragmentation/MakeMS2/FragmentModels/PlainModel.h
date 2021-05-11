#ifndef PLAIN_MODEL_H
#define PLAIN_MODEL_H

#include <memory>

#include "FragmentModel.h"

class PlainModel : public FragmentModel {
public:
   PlainModel(FragmentModelOptions& options) : FragmentModel(options) {}
   std::vector< std::vector<FragmentIon> > run(std::string peptide, const FragmentModelData inputData);
};

#endif // PLAIN_MODEL_H
