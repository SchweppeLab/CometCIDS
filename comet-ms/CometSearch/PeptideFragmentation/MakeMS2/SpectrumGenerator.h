#ifndef SPECTRUM_GENERATOR_H
#define SPECTRUM_GENERATOR_H

#include "dist.h"
#include "FragmentModels/FragmentIon.h"
#include "FragmentModels/FragmentModel.h"
#include "FragmentModels/FragmentModelFactory.h"
#include "PeptideIsotopeDist.h"

#include <map>
#include <memory>
#include <string>
#include <vector>

class SpectrumGenerator
{
public:
    void initialize();

    void initialize(FragmentModelOptions::Model model);

    void setOption(std::string key, std::string value);

    /**
     * Same as above, but calculates mz for different charge states
     *
     * @param peptide sequence without flanking residues
     * @param inputData Spectrum and Scoring data needed to generate fragments.
     */
    std::vector< std::vector<FragmentIon> > getChargeFragments(std::string peptide, FragmentModelData inputData);

    /**
     * Writes some inputs to stdout for debugging.
     */
    void printOptions();

private:
    FragmentModelOptions options;
    std::unique_ptr<FragmentModel> fragmenter;
};

#endif // SPECTRUM_GENERATOR_H
