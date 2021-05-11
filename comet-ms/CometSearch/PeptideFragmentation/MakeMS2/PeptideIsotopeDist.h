#ifndef PEPTIDE_ISOTOPE_DIST_H
#define PEPTIDE_ISOTOPE_DIST_H

#include "dist.h"

#include <map>
#include <string>
#include <vector>

class PeptideIsotopeDist {
public:
   void initialize();

   void setOption(std::string key, std::string value);

   void printOptions();

   ionDistOutput makeFragmentDist(std::string peptide, int nHeavy,  bool isUseSqrtNormalization, bool debug);

private:
   /**
    * This becomes true when AADensity0 is modified and
    * enables calculation of the fragment ion distribution with AADensity0.
    */
   bool hasNewDist = false;

   std::map<std::string, std::vector<float> > AADist;
   std::map<std::string, std::vector<float> > AADensity0;

   std::vector<float> HDist;
   std::vector<float> HODist;
   std::vector<float> H2ODist;
   std::vector<float> H3O2Dist;

   std::vector<std::vector<float> > bDeficit;
   std::vector<std::vector<float> > yDeficit;

   std::vector< std::vector<float> > averageDist(
       const std::vector< std::vector<float> > aOld,
       const std::vector< std::vector<float> > aNew,
       const std::vector< std::vector<float> > bOld,
       const std::vector< std::vector<float> > bNew,
       float p,
       float q,
       int nHeavy);

   static std::vector<std::vector<float>> normalizeIons(std::vector<std::vector<float>> unnormalizedIons,bool debug);
   static void printIonsData(std::vector<std::vector<float>> ions);
};

#endif // PEPTIDE_ISOTOPE_DIST_H
