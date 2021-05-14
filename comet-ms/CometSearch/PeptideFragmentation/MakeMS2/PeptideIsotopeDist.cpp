
#include "PeptideIsotopeDist.h"

#include "dist.h"
#include "mathlib.h"

#include <iostream>
#include <math.h>

#define MIN_PEPTIDE_DIST_P 0.0000001

using namespace std;

void PeptideIsotopeDist::initialize() {
    AADensity0 = AANaturalIsotopicDists();
    AADist = AANaturalIsotopicDists();
    hasNewDist = false;

    HDist = distAA({"H"});
    HODist = distAA({"H", "O"});
    H2ODist = distAA({"H", "H", "O"});

    H3O2Dist = distAA({"H", "H", "H", "O", "O"});

    removeDistDensity(AADensity0, H3O2Dist);
    removeDistDensity(AADist, H3O2Dist);

    bDeficit.resize(100);
    vector<string> comp = {"H", "H", "O"};
    for (int i = 0; i < 100; ++i) {
        bDeficit[i] = distAA(comp);
        comp.push_back("H");
        comp.push_back("O");
    }

    yDeficit.resize(100);
    comp.resize(0);
    comp = {"H", "H", "O", "O"};
    for (int i = 0; i < 100; ++i) {
        yDeficit[i] = distAA(comp);
        comp.push_back("H");
        comp.push_back("O");
    }
}

bool replace(std::string& str, const std::string& from, const std::string& to) {
    size_t start_pos = str.find(from);
    if(start_pos == std::string::npos)
        return false;
    str.replace(start_pos, from.length(), to);
    return true;
}

void PeptideIsotopeDist::setOption(string key, string value) {
    string aa = key;
    replace(aa, "dist-", "");
    if(AADist.find(aa) == AADist.end()) {
        throw "Invalid key in parameters" + key;
    }

    AADist[aa] = readDist(value);
    hasNewDist = true;
}

void PeptideIsotopeDist::printOptions() {
    cout << "\n\nUsing Distribution:" << endl;
    for (const auto& option : AADist)
    {
        cout << option.first << " : ";
        for (const auto& value : option.second) {
            cout << value << ",";
        }
        cout << endl;
    }
}

vector< vector<float> >
PeptideIsotopeDist::averageDist(
        const vector< vector<float> > aOld,
        const vector< vector<float> > aNew,
        const vector< vector<float> > bOld,
        const vector< vector<float> > bNew,
        float p,
        float q,
        int nHeavy) {

    int length = aOld.size();
    //int outputRows = min(length, nHeavy + 1);
    int outputRows = nHeavy+1;

    vector< vector<float> > output(outputRows, vector<float>(length));
    vector<float> pB(length);
    vector<float> qB(length);
    vector<float> pCondB(length);
    vector<float> qCondB(length);

   //    for (auto i(0); i <= nHeavy && i < length; ++i) {
    for (auto i(0); i <= nHeavy; ++i) {
        for (auto j(0u); j < length; ++j) {
            pB[j] = aOld.at(j).at(i);
        }

        for (auto j(0u); j < length; ++j) {
            qB[j] = aNew.at(j).at(i);
        }

        for (auto j(2u); j <= length; ++j) {
            pCondB[j - 2] = bOld.at(length - j).at(nHeavy - i);
            qCondB[j - 2] = bNew.at(length - j).at(nHeavy - i);
        }

        if (i == nHeavy) {
            pCondB[length - 1] = 1;
            qCondB[length - 1] = 1;
        } else {
            pCondB[length - 1] = 0;
            qCondB[length - 1] = 0;
        }

        // .5 * ((p_cb * p_b / p_) + (q_cb * q_b / q_))
        vectorMultiply(pCondB, pB);
        vectorDivide(pCondB, p);

        vectorMultiply(qCondB, qB);
        vectorDivide(qCondB, q);

        vectorAdd(pCondB, qCondB);
        vectorMultiply(pCondB, 0.5);
        output[i] = pCondB;
    }
    return output;
}

ionDistOutput PeptideIsotopeDist::makeFragmentDist(string peptide, int nHeavy, bool isUseSqrtNormalization, bool debug) {
    auto peptideLength = peptide.length();
    auto minSize = nHeavy;
    auto oldIons = ionDist(peptide, AADensity0, bDeficit, yDeficit, minSize);
    ionDistOutput newIons;
    if (hasNewDist) {
        newIons = ionDist(peptide, AADist, bDeficit, yDeficit, minSize);
    }
    else {
        newIons = oldIons;
    }

    auto oldPep = oldIons.peptide;
    auto newPep = newIons.peptide;

    float p = oldPep.at(nHeavy);
    float q = newPep.at(nHeavy);

    auto bOld = oldIons.bIons;
    auto bNew = newIons.bIons;

    auto yOld = oldIons.yIons;
    auto yNew = newIons.yIons;

    ionDistOutput output;

    //if p or q is very close to zero, peptide is not very interesting, return all zeros.
    if (static_cast<double>(p) < MIN_PEPTIDE_DIST_P || static_cast<double>(q) < MIN_PEPTIDE_DIST_P) {
        for (int i = 0; i <= nHeavy; i++) {
            vector<float> zeroVector = vector<float>(peptideLength, 0);
            output.bIons.push_back(zeroVector);
            output.yIons.push_back(zeroVector);
        }
        return (output);
    }

    output.bIons = averageDist(bOld, bNew, yOld, yNew, p, q, nHeavy);
    output.yIons = averageDist(yOld, yNew, bOld, bNew, p, q, nHeavy);

    if (isUseSqrtNormalization) {
        ionDistOutput normalizedOutput;

        normalizedOutput.bIons = normalizeIons(output.bIons, false);
        normalizedOutput.yIons = normalizeIons(output.yIons, false);

       if (debug) {
            cout << "ORIGINAL b- IONS:" << endl;
            printIonsData(output.bIons);
            cout << endl;
            cout << "NORMALIZED b- IONS:" << endl;
            printIonsData(normalizedOutput.bIons);
            cout << endl;
            cout << "ORIGINAL y- IONS:" << endl;
            printIonsData(output.yIons);
            cout << endl;
            cout << "NORMALIZED y- IONS:" << endl;
            printIonsData(normalizedOutput.yIons);
            cout << endl;
        }

        return normalizedOutput;
    } else {
        if (debug) {
            cout << "b- IONS:" << endl;
            printIonsData(output.bIons);
            cout << endl;

            cout << "y- IONS:" << endl;
            printIonsData(output.yIons);
            cout << endl;
        }

        return output;
    }

}

vector<vector<float>> PeptideIsotopeDist::normalizeIons(vector<vector<float>> unnormalizedIons, bool debug) {

    map<unsigned int, vector<float>> ionToIntensityValues = {};

    if (debug) cout << "CHARGE STATES=" << unnormalizedIons.size() << endl;

    //[1] collect all of the intensities together.

    // i = charge state number
    // j = bIon number
    for (unsigned int i = 0; i < unnormalizedIons.size(); i++) {

        vector<float> crossIonsSingleIncorporationNumber = unnormalizedIons.at(i);

        if (debug && i == 0) cout << "NUM IONS=" << crossIonsSingleIncorporationNumber.size() << endl;

        for (unsigned int j = 0; j < crossIonsSingleIncorporationNumber.size(); j ++) {

            float intensity = crossIonsSingleIncorporationNumber.at(j);

            map<unsigned int, vector<float> >::iterator it = ionToIntensityValues.find(j);

            unsigned int ion = it->first;
            if (it != ionToIntensityValues.end()) {
                if (debug) cout << "(ion, chg) = (" << j << ", " << i << ")=" << to_string(intensity) << endl;
                ionToIntensityValues[ion].push_back(intensity);
            } else {
                if (debug) cout << "(ion, chg) = (" << j << ", " << i << ")=" << to_string(intensity) << endl;
                vector<float> ionIntensityGroup;
                ionIntensityGroup.push_back(intensity);
                ionToIntensityValues.insert(make_pair(j, ionIntensityGroup));
            }
        }
    }

    if (debug) cout << "re-structured data as map." << endl;
    if (debug) cout << "map size:" << ionToIntensityValues.size() << endl;

    //[2] Now, normalize
    map<unsigned int, vector<float> > ionToIntensityValuesNormalized = {};

    map<unsigned int, vector<float > >::iterator it = ionToIntensityValues.begin();
    while (it != ionToIntensityValues.end()){

        unsigned int ion(it->first);
        vector<float> isotopomers = it->second;

        float sumSqrts = 0;

        vector<float> sqrt_isotopomers;
        for (auto isotopomer : isotopomers) {
            float sqrtIsotopomer = sqrt(isotopomer);
            sumSqrts += sqrtIsotopomer;
            sqrt_isotopomers.push_back(sqrtIsotopomer);
        }

        vector<float> sqrt_isotopomer_normalized;
        for (auto sqrt_isotopomer : sqrt_isotopomers) {
            float normalized_sqrt_isotopomer = sumSqrts > static_cast<float>(1e-10) ? sqrt_isotopomer/sumSqrts : 0;
            sqrt_isotopomer_normalized.push_back(normalized_sqrt_isotopomer);
        }

        ionToIntensityValuesNormalized.insert(make_pair(ion, sqrt_isotopomer_normalized));

        ++it;
    }

    if (debug) cout << "normalized data." << endl;

    //[3] restructure output as vector of vectors
    vector<vector<float>> ionsNormalized = unnormalizedIons;

    if (debug) cout << "ITERATION THROUGH NORMALIZED MAP." << endl;
    map<unsigned int, vector<float > >::iterator normalizedIt = ionToIntensityValuesNormalized.begin();
    while (normalizedIt != ionToIntensityValuesNormalized.end()) {

        unsigned int ion = normalizedIt->first;
        vector<float> normalizedIsotopomers = normalizedIt->second;

        if (debug) {
            cout << "IN MAP: " << ion << endl;
            for (auto f : normalizedIsotopomers){
                cout << f << " ";
            }
            cout << endl;
        }

        for (unsigned int chargeNum = 0; chargeNum < normalizedIsotopomers.size(); chargeNum++) {
            float intensity = normalizedIsotopomers.at(chargeNum);
            ionsNormalized.at(chargeNum).at(ion) = intensity;
        }

        ++normalizedIt;
    }

    if (debug) cout << "updated normalized values in vector<vector<>>." << endl;

    return ionsNormalized;
}



void PeptideIsotopeDist::printIonsData(vector<vector<float>> ions) {
    map<unsigned int, vector<float>> ionToIntensityValues = {};

    //[1] collect all of the intensities together.

    // i = charge state number
    // j = bIon number
    for (unsigned int i = 0; i < ions.size(); i++) {

        vector<float> singleIncorporationNumber = ions.at(i);

        for (unsigned int j = 0; j < singleIncorporationNumber.size(); j ++) {

            float intensity = singleIncorporationNumber.at(j);

            map<unsigned int, vector<float> >::iterator it = ionToIntensityValues.find(j);

            unsigned int ion = it->first;
            if (it != ionToIntensityValues.end()) {
                ionToIntensityValues[ion].push_back(intensity);
            } else {
                vector<float> ionIntensityGroup;
                ionIntensityGroup.push_back(intensity);
                ionToIntensityValues.insert(make_pair(j, ionIntensityGroup));
            }
        }
    }

    map<unsigned int, vector<float > >::iterator it = ionToIntensityValues.begin();
    while (it != ionToIntensityValues.end()){

        unsigned int ion = it->first;

        string ionString = to_string(ion+1);
        ionString.append(": ");

        vector<float> intensities = it->second;
        float intensitySum = 0;

        for (auto f : intensities) {
            intensitySum += f;
            ionString.append(to_string(f));
            ionString.append(" ");
        }

        cout << ionString << "(sum= " << to_string(intensitySum) << ")" << endl;

        ++it;
    }
}

