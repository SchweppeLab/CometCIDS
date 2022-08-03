
#include "FragmentModelFactory.h"

#include "ChargeIsotopeModel.h"
#include "ChargeModel.h"
#include "IsotopeModel.h"
#include "MobileProtonIsotopeModel.h"
#include "MobileProtonModel.h"
#include "ChargeConditionalFragmentModel.h"
#include "PlainModel.h"
#include "HyperFragModel.h"

using namespace std;

unique_ptr<FragmentModel> FragmentModelFactory::get(FragmentModelOptions& options)
{
    switch (options.type) {
    case FragmentModelOptions::CHARGE_ISOTOPE:
        return unique_ptr<FragmentModel>(new ChargeIsotopeModel(options));
        break;
    case FragmentModelOptions::CHARGE:
        return unique_ptr<FragmentModel>(new ChargeModel(options));
        break;
    case FragmentModelOptions::ISOTOPE:
        return unique_ptr<FragmentModel>(new IsotopeModel(options));
        break;
    case FragmentModelOptions::MOBILE_PROTON:
        return unique_ptr<FragmentModel>(new MobileProtonModel(options));
        break;
    case FragmentModelOptions::MOBILE_PROTON_ISOTOPE:
        return unique_ptr<FragmentModel>(new MobileProtonIsotopeModel(options));
        break;
    case FragmentModelOptions::MOBILE_PROTON2:
        return unique_ptr<FragmentModel>(
                    new ChargeConditionalFragmentModel(options,
                                                       unique_ptr<FragmentModel>(new MobileProtonModel(options)),
                                                       unique_ptr<FragmentModel>(new PlainModel(options)),
                                                       2)
                    );
        break;
    case FragmentModelOptions::MOBILE_PROTON_ISOTOPE2:
        return unique_ptr<FragmentModel>(
                    new ChargeConditionalFragmentModel(options,
                                                       unique_ptr<FragmentModel>(new MobileProtonIsotopeModel(options)),
                                                       unique_ptr<FragmentModel>(new IsotopeModel(options)),
                                                       2)
                    );
        break;
    case FragmentModelOptions::MOBILE_PROTON3:
        return unique_ptr<FragmentModel>(
                    new ChargeConditionalFragmentModel(options,
                                                       unique_ptr<FragmentModel>(new MobileProtonModel(options)),
                                                       unique_ptr<FragmentModel>(new PlainModel(options)),
                                                       3)
                    );
        break;
    case FragmentModelOptions::MOBILE_PROTON_ISOTOPE3:
        return unique_ptr<FragmentModel>(
                    new ChargeConditionalFragmentModel(options,
                                                       unique_ptr<FragmentModel>(new MobileProtonIsotopeModel(options)),
                                                       unique_ptr<FragmentModel>(new IsotopeModel(options)),
                                                       3)
                    );
        break;
    case FragmentModelOptions::HYPER_FRAG:
        return unique_ptr<FragmentModel>(new HyperFragModel(options));
        break;
    default:
        return unique_ptr<FragmentModel>(new PlainModel(options));
        break;
    }
}
