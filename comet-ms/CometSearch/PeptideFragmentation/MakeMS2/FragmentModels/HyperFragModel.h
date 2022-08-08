#ifndef HYPERFRAGMODEL_H
#define HYPERFRAGMODEL_H

#include "FragmentModel.h"

class HyperFragModel : public FragmentModel {
public:
    HyperFragModel(FragmentModelOptions& options) : FragmentModel(options) {}
    std::vector< std::vector<FragmentIon> > run(std::string peptide, const FragmentModelData inputData);

    static std::vector<double> hyperFrag(int ion, int peptideLength, int numDeuteria);

};

class HyperFragLookupTable {

public:
    //     <numHeavyIsotopes, peptideLength>           //vector of ions
    static std::map<std::pair<unsigned int, unsigned int>, std::vector<std::vector<double>>> hyperFragLookupTable;

    static void init_map() {
        HyperFragLookupTable::hyperFragLookupTable = std::map<std::pair<unsigned int, unsigned int>, std::vector<std::vector<double>>>();
        HyperFragLookupTable::hyperFragLookupTable.insert(std::make_pair(std::make_pair(8, 17), std::vector<std::vector<double>>{{0.5294118,0.4705882},{0.2647059,0.5294118,0.2058824},{0.1235294,0.4235294,0.3705882,0.08235294},{0.05294118,0.2823529,0.4235294,0.2117647,0.02941176},{0.02036199,0.1628959,0.3800905,0.3257919,0.10181,0.009049774},{0.00678733,0.08144796,0.2850679,0.3800905,0.2036199,0.04072398,0.002262443},{0.00185109,0.03455368,0.1814068,0.3628137,0.3023447,0.103661,0.01295763,0.0004113534},{0.000370218,0.01184698,0.09675031,0.2902509,0.3628137,0.1935006,0.04146442,0.002961744,4.113534e-05},{4.113534e-05,0.002961744,0.04146442,0.1935006,0.3628137,0.2902509,0.09675031,0.01184698,0.000370218,0},{0,0.0004113534,0.01295763,0.103661,0.3023447,0.3628137,0.1814068,0.03455368,0.00185109,0,0},{0,0,0.002262443,0.04072398,0.2036199,0.3800905,0.2850679,0.08144796,0.00678733,0,0,0},{0,0,0,0.009049774,0.10181,0.3257919,0.3800905,0.1628959,0.02036199,0,0,0,0},{0,0,0,0,0.02941176,0.2117647,0.4235294,0.2823529,0.05294118,0,0,0,0,0},{0,0,0,0,0,0.08235294,0.3705882,0.4235294,0.1235294,0,0,0,0,0,0},{0,0,0,0,0,0,0.2058824,0.5294118,0.2647059,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0.4705882,0.5294118,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0}}));
    }
};

#endif // HYPERFRAGMODEL_H
