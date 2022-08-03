#ifndef FRAGMENTMODELOPTIONS_H
#define FRAGMENTMODELOPTIONS_H

#include "PeptideFragmentation/MakeMS2/PeptideIsotopeDist.h"

class FragmentModelOptions
{
public:

    enum Model
    {
        PLAIN,
        ISOTOPE,
        CHARGE,
        CHARGE_ISOTOPE,
        MOBILE_PROTON,
        MOBILE_PROTON_ISOTOPE,
        MOBILE_PROTON2,
        MOBILE_PROTON_ISOTOPE2,
        MOBILE_PROTON3,
        MOBILE_PROTON_ISOTOPE3,
        HYPER_FRAG
    };

    FragmentModelOptions() {
        intensityCutoff = 0.05;
        obsChargeOffset = 0;
        normalize = true;
        doSqrt = true;
        doWindowedScaling = true;
        isUseSqrtNormalization = false;
        isDebug = false;
        useMzRange = false;
    }

    Model type;

    PeptideIsotopeDist isotopes;

    /**
     * Min intensity relative to base peak.
     *
     * intensity Cutoff is multiplied by base peak, and lower
     * intensity peaks are removed from theoretical fragments.
     */
    float intensityCutoff;

    /**
     * value added to max charge for fragment ions.
     */
    int obsChargeOffset;

    /**
     * Apply binned scaling to fragments.
     */
    bool normalize;

    /**
     * Use square root of expected abundance
     */
    bool doSqrt;

    /**
     * Scale fragments to 1 in each window.
     */
    bool doWindowedScaling;

    /**
     * all charge state ions for a particular b-/y- ion are normalized
     * according to sqrt(ion)/(sum(sqrt(ions))
     *
     * so, b3 with charge states 0, 1, and 2 (aka, [M+0], [M+1], and [M+2]) has new intensity values
     * for charge states 0, 1 and 2 of
     *
     * 0 = sqrt(0)/(sqrt(0)+(sqrt(1)+sqrt(2))
     * 1 = sqrt(1)/(sqrt(0)+(sqrt(1)+sqrt(2))
     * 2 - sqrt(2)/(sqrt(0)+(sqrt(1)+sqrt(2))
     */
    bool isUseSqrtNormalization;

    /**
      * special parameter for --osxtm debugging.
      * turned off by default.
      */
    bool isDebug;

    /**
     * Enable to consider only fragment ions within the observed mz range of the spectrum.
     */
    bool useMzRange;

};

#endif // FRAGMENTMODELOPTIONS_H
