#ifndef FRAGMENT_MODEL_FACTORY_H
#define FRAGMENT_MODEL_FACTORY_H

#include <memory>

#include "FragmentModel.h"

class FragmentModelFactory {
public:
   static std::unique_ptr<FragmentModel> get(FragmentModelOptions &options);
};

#endif // FRAGMENT_MODEL_FACTORY_H
