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

    std::vector<std::vector<double>> getIntensities(unsigned int numDeuteria, unsigned int peptideLength) {

        if (numDeuteria > peptideLength) {
            //when there are more deuteria than amino acids, assume one deuterium per AA.
            return hyperFragLookupTable.at(std::make_pair(peptideLength, peptideLength));
        }

        return hyperFragLookupTable.at(std::make_pair(numDeuteria, peptideLength));
    }

    static HyperFragLookupTable& instance();

private:

    //     <numHeavyIsotopes, peptideLength>           //vector of ions
    std::map<std::pair<unsigned int, unsigned int>, std::vector<std::vector<double>>> hyperFragLookupTable;

    void init_map() {

        hyperFragLookupTable = std::map<std::pair<unsigned int, unsigned int>, std::vector<std::vector<double>>>();




    }

    //private instantiation, since this is a singleton.
    HyperFragLookupTable(){
        init_map();
    }

};

#endif // HYPERFRAGMODEL_H
