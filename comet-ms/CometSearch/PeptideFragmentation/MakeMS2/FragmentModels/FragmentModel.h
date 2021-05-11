#ifndef FRAGMENT_MODEL_H
#define FRAGMENT_MODEL_H

#include "FragmentIon.h"
#include "FragmentModelOptions.h"
#include "FragmentModelData.h"
#include "PeptideFragmentation/MakeMS2/PeptideIsotopeDist.h"

#include <string>
#include <vector>

class FragmentModel {
public:
   FragmentModel(FragmentModelOptions& opt) : options(opt) {}

   /**
    * @brief run
    *
    * @param peptide The peptide sequence
    * @param inputData Data from each spectrum or scoring event needed to generate fragments.
    */
   virtual std::vector< std::vector<FragmentIon> > run(std::string peptide, const FragmentModelData inputData) = 0;

   /**
    * Scale and filter fragment ions.
    *
    * @param fragments
    */
   void normalize(std::vector< std::vector<FragmentIon> > &fragments);

protected:
   FragmentModelOptions& options;
};

#endif // FRAGMENT_MODEL_H
