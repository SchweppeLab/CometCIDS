#ifndef PEPTIDE_MASS_H
#define PEPTIDE_MASS_H

#include <string>
#include <vector>

// Fragment Types
#define MASSES_B 10
#define MASSES_Y 20

#define DEUTERIUM_MASS_DIFF 1.0057281805
/*
 * Note that 6*DEUTERIUM_MASS_DIFF = 6.034369083
 *
 * 80Se - 74Se (also 6 neutrons) = 5.9940449
 *
 * delta m/z = 0.040324183
 *
 * ppm gap of delta m/z:
 * at m/z = 300: 134 ppm
 * at m/z = 500: 80 ppm
 * at m/z = 700: 57 ppm
 *
 * For reference,
 * delta 13C - 15N = 0.00632010556
 *
 * A peptide fragment with an 80Se (instead of 74Se) and no deuteria should be distinguishable
 * from a peptide fragment with 6 deuteria. With a 20 ppm fragment match tolerance, will probably miss the fragments.
 */
std::vector<double> generateFwd(const std::string peptide);

std::vector<double> generateRev(const std::string peptide);

void generateFragments(const std::string peptide, int fragmentType, std::vector<double>& output);

float masstomz(float mass, int charge);

#endif // PEPTIDE_MASS_H
