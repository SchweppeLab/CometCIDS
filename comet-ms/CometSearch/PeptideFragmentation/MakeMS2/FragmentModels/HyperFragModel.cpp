#include "HyperFragModel.h"

#include "../../../../../extern/Rmath/dhyper.c"
#include "../../../../../extern/Rmath/dbinom.c"
#include "../../../../../extern/Rmath/stirlerr.c"
#include "../../../../../extern/Rmath/bd0.c"

using namespace std;

int R_finite(double x){
    return isfinite(x);
}

vector<double> HyperFragModel::hyperFrag(int ion, int peptideLength, int numDeuteria){

    vector<double> intensities(numDeuteria+1);

    for (unsigned int i = 0; i <= intensities.size(); i++) {
        auto intensity = dhyper(i, numDeuteria, (peptideLength-numDeuteria), ion, false);
        if (isfinite(intensity)) {
            intensities[i] = intensity;
        }
    }

    return intensities;
}

vector< vector<FragmentIon> > HyperFragModel::run(string peptide, const FragmentModelData inputData) {
    //TODO
}
