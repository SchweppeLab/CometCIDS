#pragma once
#include <numbers>
class HyperGeometric
{
public:
     /**
     * Get instance of singleton.
     *
     * @return FragmentCalculator
     */
    static HyperGeometric& instance();

    float binomialCoeff(float n, float k);

    float pmf(float n, float k, float N, float K);

    vector<float> hyperFrag(int i, int n, int h);


private:
    /**
     * Make the constructor private since this is a singleton.
     */
    HyperGeometric();

    /**
     * Disable the copy constructor to make sure
     * the singleton isnt accidentally copied.
     */
    HyperGeometric(const HyperGeometric&) = delete;
};

