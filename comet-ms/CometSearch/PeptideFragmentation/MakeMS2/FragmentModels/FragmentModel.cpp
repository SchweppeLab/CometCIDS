
#include "FragmentModel.h"

#include <vector>
#include <math.h>

using namespace std;

void FragmentModel::normalize(vector< vector<FragmentIon> > &fragments) {

    if (options.doSqrt) {
        // Take sqrt of fragment intensities.
        for (auto& charge : fragments) {
            for (auto& fragment : charge) {
                fragment.intensity = sqrt(fragment.intensity);
            }
        }
    }

    // Find max
    double maxMz = 0;
    double maxIntensity = 0;
    for (const auto charge : fragments) {
        for (const auto fragment : charge) {
            if(fragment.mz > maxMz) {
                maxMz = fragment.mz;
            }
            if(fragment.intensity > maxIntensity) {
                maxIntensity = fragment.intensity;
            }
        }
    }

    // filter
    double minIntensity = maxIntensity * options.intensityCutoff;
    vector< vector<FragmentIon> > filtered(fragments.size());
    for (int i = 0; i < fragments.size(); ++i) {
        for (const auto fragment : fragments.at(i)) {
            if(fragment.intensity > minIntensity) {
                filtered[i].push_back(fragment);
            }
        }
    }

    if(options.doWindowedScaling) {
        // find window maxs
        vector<double> binMax(11);

        double windowSize = maxMz / 10;
        for (const auto charge : filtered) {
            for (const auto fragment : charge) {
                int bin = fragment.mz / windowSize;
                if (binMax.at(bin) < fragment.intensity) {
                    binMax[bin] = fragment.intensity;
                }
            }
        }

        double maxOfAllBins = -1;

        for (auto binVal: binMax){
            if (binVal > maxOfAllBins){
                maxOfAllBins = binVal;
            }
        }

        // apply scaling
        for (int i = 0; i < filtered.size(); ++i) {
            for (auto& fragment : filtered[i]) {
                int bin = fragment.mz / windowSize;

                /**
                * Changed by Jon + Phil on 2019-07-11
                */
                //fragment.intensity /= binMax[bin];
                fragment.intensity /= maxOfAllBins;
            }
        }
    }

    fragments = filtered;
}
