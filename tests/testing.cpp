
#include "testing.h"

#include <cmath>
#include <limits>
#include <sys/stat.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <unistd.h>



using namespace std;

void serializeDist(map<string, vector<float> > dist, string filename) {
    ofstream out;
    out.open(filename);
    for (auto aa : dist) {
        out << "dist_" << aa.first << " = ";
        bool first = true;
        for (auto x : aa.second) {
            if (!first) {
                out << ",";
            }
            first = false;
            out << x;
        }
        out << endl;
    }
}

/**
 * @brief unserializeDist
 *
 * This is suited to parse the distribution from the comet params.
 */
void unserializeDist(string filename) {
//    ifstream infile(filename);

//    while (infile)
//    {
//        string s;
//        if (!getline( infile, s )) {
//            break;
//        }

//        istringstream ss(s);
//        vector <string> record;

//        while (ss)
//        {
//            string s;
//            if (!getline( ss, s, ',')) {
//                break;
//            }
//            record.push_back( s );
//        }

//        // data.push_back( record );
//    }
}

std::vector<FragmentIon>
flattenFragments(const std::vector<std::vector<FragmentIon> > &in) {
    std::vector<FragmentIon> out;
    for (auto part : in) {
        out.insert(out.end(), part.begin(), part.end());
    }
    return out;
}

void writeFragments(const std::vector<FragmentIon> &in, string filename) {
    ofstream out;
    out.open(filename);
    out << serializeFragments(in);
    out.close();
}

string serializeFragments(const vector<FragmentIon> &in) {
    stringstream ss;
    for (auto fragment : in) {
        ss << setprecision(8) << fragment.mz << "\t" << fragment.intensity << "\n";
    }
    return ss.str();
}

inline bool file_exists (const std::string& name) {
    struct stat buffer;
    return (stat (name.c_str(), &buffer) == 0);
}

vector<FragmentIon> unserializeFragments(string filename) {
    vector<FragmentIon> result;
    ifstream in(filename);
    if (in.is_open()) {
        FragmentIon fragment;
        while(in >> fragment.mz >>fragment.intensity) {
            result.push_back(fragment);
        }
    }
    else {
        cout << "Could not open file: " << filename << endl;
        if  (!file_exists(filename)) {
            cout << "File does not exist" << endl;

            char buffer[FILENAME_MAX];
            getcwd(buffer, FILENAME_MAX);
            cout << "current dir: " << buffer << endl;
        }
    }
    return result;
}

/**
 * Compares two floating point numbers, only checking
 * for accuracy in some of the leading digits specified by
 * the precision parameter.
 *
 * @return
 */
bool isEqual(double x, double y, int precision=7) {
    return fabs(x - y) < (max(x,y) * pow(10, -1 * precision));
}

/**
 * @brief Compare two lists of fragments and return true if they are the same.
 * This is written to allow for showing a diff later.
 *
 * @param a
 * @param b
 * @return
 */
bool compareFragments(const vector<FragmentIon> &a, const vector<FragmentIon> &b) {
    if (a.size() != b.size()) {
        return false;
    }
    for (unsigned int i = 0; i < a.size(); ++i) {
        if (!isEqual(a[i].mz, b[i].mz)
            || !isEqual(a[i].intensity, b[i].intensity)) {
            return false;
        }
    }
    return true;
}

void printDifference(vector<FragmentIon> a, vector<FragmentIon> b) {
    int i = 0;
    int j = 0;
    while (true) {
        if (i < a.size()) {
            cout << a.at(i).mz << "\t" << a.at(i).intensity << "\t";
            ++i;
        }
        else {
            cout << "\t\t";
        }

        if (j < b.size()) {
            cout << b.at(j).mz << "\t" << b.at(j).intensity << "\t";
            ++j;
        }
        else {
            cout << "\t\t";
        }

        cout << endl;
        if (i >= a.size() && j >= b.size()) {
            break;
        }
    }
}
