#ifndef FRAGMENTMODELDATA_H
#define FRAGMENTMODELDATA_H


class FragmentModelData
{
public:
    /**
     * The observed offset from the number of heavy isotopes.
     */
    int nHeavy;

    /**
     * The observed precursor charge.
     */
    int obsCharge;

    /**
     * The max charge to use when creating fragment ions.
     */
    int maxCharge;

    /**
     * Minimum m/z range to consider.
     */
    double minMz;

    /**
     * Maximum m/z range to consider.
     */
    double maxMz;

    /**
     * Masses at each position along the peptide going forward.
     * This is precomputed by comet.
     */
    double* aaFwdMasses;

    /**
     * Masses at each position along the peptide going reverse.
     * This is precomputed by comet.
     */
    double* aaRevMasses;
};

#endif // FRAGMENTMODELDATA_H
