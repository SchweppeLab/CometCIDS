#include "HyperFragModel.h"

#include <fstream>
#include <iostream>

#ifdef __linux__ 
#include <Resources.h>
#endif

#ifndef RESOURCE_PATH
#define RESOURCE_PATH "../../../../data"
#endif

#ifdef _WIN32
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
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
    
    const string path_sep =
#ifdef _WIN32
        "\\";
#else
        "/";
#endif

    hyperFragLookupTable = std::map<std::pair<unsigned int, unsigned int>, std::vector<std::vector<double>>>();

    const string LOOKUP_TABLE_FILE = path_sep + "hyperfrag_dist_lookup_table.txt";

    cout << "File name: \'" << LOOKUP_TABLE_FILE << "\'" << endl;

    string LOOKUP_TABLE_FULL_PATH = RESOURCE_PATH + LOOKUP_TABLE_FILE;
    ifstream lookupTableFileStream(LOOKUP_TABLE_FULL_PATH);

    if (!lookupTableFileStream.is_open()) {

        cout << "Unable to find hyperfrag_dist_lookup_table.txt file: \'" <<  LOOKUP_TABLE_FULL_PATH << "\'" << endl;

        char buff[FILENAME_MAX]; //create string buffer to hold path
        GetCurrentDir( buff, FILENAME_MAX );
        string current_working_dir(buff);
        LOOKUP_TABLE_FULL_PATH = current_working_dir + path_sep + "data" + LOOKUP_TABLE_FILE;

        cout << "Looking for hyperfrag_dist_lookup_table.txt file: \'" << LOOKUP_TABLE_FULL_PATH << "\'" << endl;
        lookupTableFileStream.open(LOOKUP_TABLE_FULL_PATH);

    }

    string line;
    if (lookupTableFileStream.is_open()) {

        cout << "Successfully found hyper frag dist lookup table file: \'" << LOOKUP_TABLE_FULL_PATH << "\'" << endl;

        unsigned int counter = 0;
        unsigned int N = 100;
        while(getline(lookupTableFileStream, line)) {

            //split by \t
            vector<string> bits{};

            HyperFragLookupTable::split(line, "\t", bits);

            int numHeavy = stoi(bits[0]);
            int pepLength = stoi(bits[1]);

            vector<vector<double>> fragDist = HyperFragLookupTable::decodeFragDist(bits[2]);

            hyperFragLookupTable.insert(make_pair(make_pair(numHeavy, pepLength), fragDist));

            counter++;

            if (counter % N == 0) {
                cout << "Processed " << counter << " HyperFrag Distributions." << endl;
            }
        }
    } else {
        cerr << "Unable to locate hyper frag dist lookup table file.\n"
                "Please ensure that file data/hyperfrag_dist_lookup_table.txt exists and is readable."
             << endl;
        abort();
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
