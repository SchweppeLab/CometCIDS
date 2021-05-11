#ifndef FRAGMENT_CALCULATOR_H
#define FRAGMENT_CALCULATOR_H

#include <map>
#include <string>
#include <vector>

#include "MakeMS2/SpectrumGenerator.h"

class FragmentCalculator {
public:

    /**
     * Get instance of singleton.
     *
     * @return FragmentCalculator
     */
    static FragmentCalculator& instance();

    void setOption(std::string key, std::string value);

    void printOptions();

    std::vector<std::string> validOptionChars();

    void generateFragments(char* peptide, int peptideLength, FragmentModelData inputData, std::vector<std::vector<FragmentIon> >& output);
    SpectrumGenerator& getGenerator() { return generator; }
private:
    /**
     * Make the constructor private since this is a singleton.
     */
    FragmentCalculator();

    /**
     * Disable the copy constructor to make sure
     * the singleton isnt accidentally copied.
     */
    FragmentCalculator(const FragmentCalculator&) = delete;

    double getFragmentIonMass(int iWhichIonSeries,
                            int i,
                            int ctCharge,
                            double *_pdAAforward,
                            double *_pdAAreverse);
    double masstomz(double mass, int charge);
    double mztomass(double mz, int charge);

    SpectrumGenerator generator;
};

#endif // FRAGMENT_CALCULATOR_H
