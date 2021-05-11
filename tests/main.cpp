
#include "SpectrumGenerator.h"
#include "peptidemass.h"
#include "testing.h"
#include "hyper.h"
#include "benchmark.h"
#include "FragmentModels/MobileProtonModel.h"
#include "FragmentModels/MobileProtonIsotopeModel.h"

#include <cmath>
#include <iostream>
#include <fstream>
#include <string.h>

using namespace std;

void testFragments(string testFile) {
    SpectrumGenerator generator;
    generator.initialize();

    string peptide = "EQGLR";
    auto fwd = generateFwd(peptide);
    auto rev = generateRev(peptide);

    FragmentModelData data;
    data.nHeavy = 2;
    data.obsCharge = 2;
    data.maxCharge = 1;
    data.minMz = 0;
    data.maxMz = 2000;
    data.aaFwdMasses = fwd.data();
    data.aaRevMasses = rev.data();
    auto fragments = generator.getChargeFragments(peptide, data);
    vector<FragmentIon> result = flattenFragments(fragments);

    // For generating the test data:
    // serializeFragments(result, "data/expected.tsv");
    auto expected = unserializeFragments(testFile);
    if(!compareFragments(result, expected)) {
        printDifference(result, expected);
        cout << "error: result does not match expected" << endl;
        throw 1;
    }
    cout << "results match" << endl;
}

void showFragments(string peptide, ofstream& outputFileStream) {

    SpectrumGenerator generator;
    generator.initialize();

    auto fwd = generateFwd(peptide);
    auto rev = generateRev(peptide);

    FragmentModelData data;
    data.nHeavy = 2;
    data.obsCharge = 2;
    data.maxCharge = 1;
    data.minMz = 0;
    data.maxMz = 2000;
    data.aaFwdMasses = fwd.data();
    data.aaRevMasses = rev.data();
    auto fragments = generator.getChargeFragments(peptide, data);
    vector<FragmentIon> result = flattenFragments(fragments);

    outputFileStream << "mz\tintensity\n";
    for (auto f : result) {
        outputFileStream << to_string(f.mz) << "\t" << to_string(f.intensity) << "\n";
    }

}

void runBenchmark() {
    SpectrumGenerator generator;
    generator.initialize();
    string peptide = "EQGLR";
    auto fwd = generateFwd(peptide);
    auto rev = generateRev(peptide);

    FragmentModelData data;
    data.nHeavy = 2;
    data.obsCharge = 2;
    data.maxCharge = 1;
    data.minMz = 0;
    data.maxMz = 2000;
    data.aaFwdMasses = fwd.data();
    data.aaRevMasses = rev.data();

    BENCH("5aa2z10k", 10000) {
        auto fragments = generator.getChargeFragments(peptide, data);
    }
}

void testHyper() {
    auto p1 = dhyper(1, 2, 2, 3);
    if(p1 != 0.5) {
        cout << "error: dhyper1 is " << p1 << endl;
        throw 1;
    }
    auto p2 = dhyper(14, 70, 30, 20);
    if(fabs(p2 - 0.2140911) > 0.000001) {
        cout << "error: dhyper2 is " << p2 << endl;
        throw 1;
    }
    auto p3 = dhyper(1, 3, 3, 3);
    if(p3 != 0.45) {
        cout << "error: dhyper3 is " << p3 << endl;
        throw 1;
    }
}

void testMobileProton(string peptide, int obsCharge, string testDir, string testFilename) {
    FragmentModelOptions options;
    options.isotopes.initialize();
    MobileProtonModel model(options);

    auto fwd = generateFwd(peptide);
    auto rev = generateRev(peptide);

    FragmentModelData data;
    data.obsCharge = obsCharge;
    data.maxCharge = 6;
    data.aaFwdMasses = fwd.data();
    data.aaRevMasses = rev.data();
    data.nHeavy = 0;

    auto fragments = model.run(peptide, data);
    auto result = flattenFragments(fragments);

    // writeFragments(result, testDir + testFilename);
    vector<FragmentIon> expected;
    if (!testFilename.empty()) {
        expected = unserializeFragments(testDir + testFilename);
    }
    if(!compareFragments(result, expected)) {
        printDifference(result, expected);
        cout << "error: testMobileProton: result does not match expected for " << peptide << " at +" << data.obsCharge << endl;
        throw 1;
    }
}

void testMobileProtonIsotope(string peptide, int obsCharge, string testDir, string testFilename) {
    FragmentModelOptions options;
    options.isotopes.initialize();
    MobileProtonIsotopeModel model(options);

    auto fwd = generateFwd(peptide);
    auto rev = generateRev(peptide);

    FragmentModelData data;
    data.obsCharge = obsCharge;
    data.maxCharge = 6;
    data.aaFwdMasses = fwd.data();
    data.aaRevMasses = rev.data();
    data.nHeavy = 0;

    auto fragments = model.run(peptide, data);
    auto result = flattenFragments(fragments);

    // writeFragments(result, testDir + testFilename);
    vector<FragmentIon> expected;
    if (!testFilename.empty()) {
        expected = unserializeFragments(testDir + testFilename);
    }
    if(!compareFragments(result, expected)) {
        printDifference(result, expected);
        cout << "error: testMobileProtonIsotope: result does not match expected for " << peptide << " at +" << data.obsCharge << endl;
        throw 1;
    }
}

int main(int argc, char *argv[]) {

    cout << "start" << endl;

    string testFile = "data/expected.tsv";
    string testDir = "data/";

    if (argc > 2) {
        if (strcmp(argv[1], "-d") == 0){
            testFile = string(argv[2]);
            cout << "Test file located at: \"" << testFile << "\"" << endl;
        } else if (strcmp(argv[1], "-peptide") == 0) {

            string peptide = argv[2];
            string outputFile = argv[3];

            cout << "Output file located at: \"" << outputFile << "\"" << endl;

            ofstream outputFileStream;
            outputFileStream.open(outputFile);

            showFragments(peptide, outputFileStream);
            outputFileStream.close();
            cout << "All Processes Completed Successfully!" << endl;
            cout << "done" << endl;
            return 0;
        }
    }

    cout << "test file is " << testFile << endl;

    testFragments(testFile);
    testMobileProton("EQGLR", 2, testDir, "mobile_proton1.tsv");
    testMobileProton("EQGHGLR", 3, testDir, "mobile_proton2.tsv");
    testMobileProton("GGKGGR", 2, testDir, "mobile_proton3.tsv");
    testMobileProton("GGR", 100, testDir, "");
    testMobileProtonIsotope("PEPTIDER", 2, testDir, "mobile_proton_isotope1.tsv");
    testMobileProtonIsotope("GGR", 100, testDir, "");
    runBenchmark();
    testHyper();

    string peptide = "EQHHGLR";
    auto fwd = generateFwd(peptide);
    auto rev = generateRev(peptide);

    FragmentModelData data;
    data.nHeavy = 0;
    data.obsCharge = 3;
    data.maxCharge = 3;
    data.minMz = 0;
    data.maxMz = 2000;
    data.aaFwdMasses = fwd.data();
    data.aaRevMasses = rev.data();

    SpectrumGenerator generator2;
    generator2.initialize(FragmentModelOptions::CHARGE);
    auto result3 = generator2.getChargeFragments(peptide, data);

    cout << "done" << endl;
    return 0;
}
