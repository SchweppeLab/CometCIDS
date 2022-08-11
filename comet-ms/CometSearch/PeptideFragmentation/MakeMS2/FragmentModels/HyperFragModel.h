#ifndef HYPERFRAGMODEL_H
#define HYPERFRAGMODEL_H

#include "FragmentModel.h"
#include <iostream>

class HyperFragModel : public FragmentModel {
public:
    HyperFragModel(FragmentModelOptions& options) : FragmentModel(options) {}
    std::vector< std::vector<FragmentIon> > run(std::string peptide, const FragmentModelData inputData);

    static std::vector<double> hyperFrag(int ion, int peptideLength, int numDeuteria);

};

class HyperFragLookupTable {

private:

    //     <numHeavyIsotopes, peptideLength>           //vector of ions
    std::map<std::pair<unsigned int, unsigned int>, std::vector<std::vector<double>>> hyperFragLookupTable;

public:

    std::vector<std::vector<double>> getIntensities(unsigned int numDeuteria, unsigned int peptideLength) {

        if (numDeuteria > peptideLength) {
            //when there are more deuteria than amino acids, assume one deuterium per AA.
            return hyperFragLookupTable.at(std::make_pair(peptideLength, peptideLength));
        }

        std::pair<int, int> key = std::make_pair(numDeuteria, peptideLength);

        if (hyperFragLookupTable.find(key) != hyperFragLookupTable.end()) {
            return hyperFragLookupTable.at(key);
        } else {
            std::cerr << " INVALID (numDeuteria, peptideLength) key: (" << numDeuteria << ", " << peptideLength << ")" << std::endl;
            abort();
        }


    }

    static HyperFragLookupTable& instance();
    int getSize() {return hyperFragLookupTable.size();}

private:

    void init_map();
    static std::vector<std::vector<double>> decodeFragDist(std::string encodedFragDist);
    static void split(const std::string& s, std::string delim, std::vector<std::string>& v);

    //private instantiation, since this is a singleton.
    HyperFragLookupTable(){
        init_map();
    }

};

#endif // HYPERFRAGMODEL_H
