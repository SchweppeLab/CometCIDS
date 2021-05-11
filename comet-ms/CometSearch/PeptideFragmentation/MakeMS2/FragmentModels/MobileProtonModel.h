#ifndef MOBILEPROTONMODEL_H
#define MOBILEPROTONMODEL_H

#include "FragmentModel.h"

#include <set>
#include <string>
#include <vector>

class MobileProtonModel : public FragmentModel {
public:

   MobileProtonModel(FragmentModelOptions& options) : FragmentModel(options) {}

   std::vector< std::vector<FragmentIon> > run(std::string peptide, const FragmentModelData inputData);

protected:

   struct Range {
       int min = 0;
       int max = 0;
   };

   struct FragmentCharges{
       std::vector<Range> b;
       std::vector<Range> y;
   };

   struct PeptideInfo {
       int countR = 0;
       int countK = 0;
       int countH = 0;
       int countMobile = 0;
       char mobileAA;
       std::string sequesteredAA;
   };

   FragmentCharges getFragmentCharges(std::string peptide, PeptideInfo peptideInfo);

   PeptideInfo readPeptideInfo(std::string peptide, int charge);

};


#endif // MOBILEPROTONMODEL_H
