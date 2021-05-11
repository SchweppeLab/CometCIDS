
#include <boost/math/distributions/hypergeometric.hpp>
#include <algorithm>

using namespace std;

double dhyper(int x, int m, int n, int k) {
    int low = max(0, k - n);
    int high = min(k, m);
    if (x < low || x > high) {
        return 0;
    }
    // param 1: "r" the number of possible fail outcomes
    // param 2: "n" sample size
    // param 3: "N" total population size
    boost::math::hypergeometric_distribution<> dist(m, k, m + n);
    return boost::math::pdf(dist, x);
}
