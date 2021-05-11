#ifndef MOBILEPROTONISOTOPEMODEL_H
#define MOBILEPROTONISOTOPEMODEL_H

#include "FragmentModel.h"
#include "MobileProtonModel.h"

#include <set>
#include <string>
#include <vector>

class MobileProtonIsotopeModel : public MobileProtonModel {
public:
   MobileProtonIsotopeModel(FragmentModelOptions& options) : MobileProtonModel(options) {}

   std::vector< std::vector<FragmentIon> > run(std::string peptide, const FragmentModelData inputData);
};

#endif // MOBILEPROTONISOTOPEMODEL_H
