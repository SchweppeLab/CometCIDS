#include "HyperGeometric.h"
#include "CometDataInternal.h"
#include <vector>
using namespace std;

HyperGeometric::HyperGeometric() {

}

HyperGeometric& HyperGeometric::instance()
{
    static HyperGeometric hgDist;
    return hgDist;
}

/* https://www.geeksforgeeks.org/binomial-coefficient-dp-9/
 * Recursive bionmial calculation
*/
float HyperGeometric::binomialCoeff(float n, float k)
{
    // Base Cases
    if (k > n)
        return 0;
    if (k == 0 || k == n)
        return 1;

    // Recur
    return binomialCoeff(n - 1, k - 1)
        + binomialCoeff(n - 1, k);
}

/*
 * Generate the probability mass function based on recursive binomial
 */
float HyperGeometric::pmf(float n, float k, float N, float K)
{
    return binomialCoeff(K, k) * binomialCoeff(N - K, n - k) / binomialCoeff(N, n);
}

/*
 * Generate the HG probability mass function based on recursive binomial
 * i is the b or y ion subscript.  
 * n is number of amino acids
 * h is the number of modifications/isotopes observed
 */
vector<float> HyperGeometric::hyperFrag(int i, int n, int h)
{
    int j;
    vector<float> isotopeDistribution;
    for (j = 0; j <= h; j++) {
        isotopeDistribution[i] = pmf(j, h, n - j, i);
    }
    return isotopeDistribution;
}
