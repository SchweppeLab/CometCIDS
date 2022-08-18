#include "HyperFragModel.h"

#include <fstream>
#include <iostream>

//#include <Resources.h>

#ifndef RESOURCE_PATH
#define RESOURCE_PATH "data/"
#endif

using namespace std;

void HyperFragLookupTable::split(const string& s, string delimiter, vector<string>& v) {

    string::size_type posPrevious = 0;
    string::size_type posCurrent = s.find(delimiter, posPrevious);

    while (posCurrent != string::npos) {
        string val = s.substr(posPrevious, posCurrent-posPrevious);
        posPrevious = posCurrent + delimiter.length();

        if (!val.empty()) {
            v.push_back(val);
        }

        posCurrent = s.find(delimiter, posPrevious);
    }

    //last entry
    if (posCurrent == string::npos && posPrevious != s.length()) {
        string val = s.substr(posPrevious, s.length()-posPrevious+1);
        if (!val.empty()) {
            v.push_back(val);
        }
    }
}

void HyperFragLookupTable::init_map() {
    
    hyperFragLookupTable = std::map<std::pair<unsigned int, unsigned int>, std::vector<std::vector<double>>>();

    const string LOOKUP_TABLE_FILE = "hyperfrag_dist_lookup_table.txt";
    ifstream lookupTableFileStream(RESOURCE_PATH + LOOKUP_TABLE_FILE);

    cout <<"Started hyper frag lookup table parsing at: " << LOOKUP_TABLE_FILE << endl;

    string line;
    if (lookupTableFileStream.is_open()) {

        unsigned int counter = 0;
        unsigned int N = 100;
        while(getline(lookupTableFileStream, line)) {
            //split by \t
            vector<string> bits{};

            HyperFragLookupTable::split(line, ";", bits);

            int numHeavy = stoi(bits[0]);
            int pepLength = stoi(bits[1]);

            vector<vector<double>> fragDist = HyperFragLookupTable::decodeFragDist(bits[2]);

            hyperFragLookupTable.insert(make_pair(make_pair(numHeavy, pepLength), fragDist));

            counter++;

            if (counter % N == 0) {
                cout << "Processed " << counter << " HyperFrag Distributions." << endl;
            }
        }
    }
    lookupTableFileStream.close();
    cout << hyperFragLookupTable.size() << " HyperFrag Distributions imported into lookup table." << endl;
}

vector<vector<double>> HyperFragLookupTable::decodeFragDist(string encodedFragDist) {

    string encodedFragDistNoEnds = encodedFragDist.substr(2, encodedFragDist.size()-4); // remove {{ and }} at beginning and end

    vector<string> bits{};
    HyperFragLookupTable::split(encodedFragDistNoEnds, "},{", bits);

    vector<vector<double>> output = vector<vector<double>>(bits.size());

    for (unsigned int i = 0; i < bits.size(); i++) {
        vector<string> dist{};
        HyperFragLookupTable::split(bits[i], ",", dist);

        vector<double> distAsDouble(dist.size());
        for (unsigned int j = 0; j < dist.size(); j++) {
            distAsDouble[j] = stod(dist[j]);
        }

        output[i] = distAsDouble;
    }

    return output;
}

vector<vector<double>> HyperFragLookupTable::getIntensities(unsigned int numDeuteria, unsigned int peptideLength) {

    pair<int, int> key = make_pair(numDeuteria, peptideLength);

    if (hyperFragLookupTable.find(key) != hyperFragLookupTable.end()) {
        return hyperFragLookupTable.at(key);
    } else {
        cerr << "HyperFragLookupTable::getIntensities() INVALID (numDeuteria, peptideLength) key: (" << numDeuteria << ", " << peptideLength << ")" << endl;
        cerr << "Peptides must be between 1-63 AA, number of heavy isotopes limited to no more than 63." << endl;
        abort();
    }
}

HyperFragLookupTable& HyperFragLookupTable::instance() {
    static HyperFragLookupTable lookupTable = HyperFragLookupTable();
    return lookupTable;
}
