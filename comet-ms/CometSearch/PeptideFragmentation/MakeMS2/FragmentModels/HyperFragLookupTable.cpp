#include "HyperFragModel.h"

#include <fstream>
#include <iostream>

using namespace std;

void HyperFragLookupTable::split(const string& s, string delimiter, vector<string>& v) {

    unsigned long posPrevious = 0;
    unsigned long posCurrent = s.find(delimiter, posPrevious);

    while (posCurrent != string::npos) {

        string val = s.substr(posPrevious, posCurrent-posPrevious);
        posPrevious = posCurrent + delimiter.length();

        if (!val.empty()) {
            v.push_back(val);
        }

        posCurrent = s.find(delimiter, posPrevious);

        //last entry
        if (posCurrent == string::npos && posPrevious != s.length()) {
            string val = s.substr(posPrevious, s.length()-posPrevious+1);
            if (!val.empty()) {
                v.push_back(val);
            }
        }
    }
}

void HyperFragLookupTable::init_map() {
    
    hyperFragLookupTable = std::map<std::pair<unsigned int, unsigned int>, std::vector<std::vector<double>>>();
    
    //TODO
    string hyperFragLookupTableFile = "/Users/phillipseitzer/Projects/CometCIDS-issues/Issue7-DHyper/lookup_table.tsv";

    ifstream lookupTableFileStream(hyperFragLookupTableFile);
    string line;
    if (lookupTableFileStream.is_open()) {

        unsigned int counter = 0;
        unsigned int N = 100;
        while(getline(lookupTableFileStream, line)) {

            //split by \t
            vector<string> bits{};

            HyperFragLookupTable::split(line, "\t", bits);

            int numHeavy = stoi(bits[0]);
            int pepLength = stoi(bits[1]);

            if (numHeavy == 8 && pepLength == 17) {
                cout << "breakpoint!";
            }
            vector<vector<double>> fragDist = HyperFragLookupTable::decodeFragDist(bits[2]);

            hyperFragLookupTable.insert(make_pair(make_pair(numHeavy, pepLength), fragDist));

            counter++;

            if (counter % N == 0) {
                cout << "Processed " << counter << " HyperFrag Distributions." << endl;
            }
        }
    }
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
