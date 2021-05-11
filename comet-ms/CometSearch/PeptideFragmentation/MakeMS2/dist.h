#ifndef DIST_H
#define DIST_H

#include <map>
#include <string>
#include <vector>

struct ionDistOutput {
    std::vector<std::vector<float> > bIons;
    std::vector<std::vector<float> > yIons;
    std::vector<float> peptide;
};

std::vector<float> distAA(const std::vector<std::string>& x);

/**
 * Calculate Isotopic Distributions for each amino acid
 * using the natural abundance isotope distributions.
 */
std::map<std::string, std::vector<float> > AANaturalIsotopicDists();

/**
 * Calculate Isotopic Distributions for each amino acid
 * With a chance of incorporating a deuterium.
 */
std::map<std::string, std::vector<float> > AADeuteratedIsotopicDists(float p);

ionDistOutput ionDist(std::string peptide,
        const std::map<std::string, std::vector<float> >& AAdensity,
        const std::vector<std::vector<float> > &bDeficit,
        const std::vector<std::vector<float> > &yDeficit, int minSize);

void removeDistDensity(
      std::map<std::string, std::vector<float> >& dist,
      std::vector<float> remove);

std::vector<float> readDist(std::string in);

ionDistOutput makeFragmentDist(std::string peptide, int nHeavy);

#endif // DIST_H
