// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2006-2016 Tiago de Paula Peixoto <tiago@skewed.de>
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 3
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#ifndef UTIL_HH
#define UTIL_HH

#include "config.h"

#include <cmath>
#include <iostream>
#include <boost/math/special_functions/zeta.hpp>

#include "cache.hh"

namespace graph_tool
{
using namespace boost;

// Warning: lgamma(x) is not thread-safe! However, since in the context of this
// program the outputs should _always_ be positive, this can be overlooked.

inline double lbinom(double N, double k)
{
    if (N == 0 || k == 0 || k > N)
        return 0;
    return (lgamma(N + 1) - lgamma(k + 1)) - lgamma(N - k + 1);
}

inline double lbinom_fast(int N, int k)
{
    if (N == 0 || k == 0 || k > N)
        return 0;
    return lgamma_fast(N + 1) - lgamma_fast(N - k + 1) - lgamma_fast(k + 1);
}

inline double lbinom_careful(double N, double k)
{
    if (N == 0 || k == 0 || k >= N)
        return 0;
    double lgN = lgamma(N + 1);
    double lgk = lgamma(k + 1);
    if (lgN - lgk > 1e8)
    {
        // We have N >> k. Use Stirling's approximation: ln N! ~ N ln N - N
        // and reorder
        return - N * log1p(-k / N) - k * log1p(-k / N) - k - lgk + k * log(N);
    }
    else
    {
        return lgN - lgamma(N - k + 1) - lgk;
    }
}

} // namespace graph_tool


// polylogarithm and degree-distribution description length (xi)

double spence(double x);

template <class Type>
Type polylog(int n, Type z, Type epsilon=1e-6)
{
    if (n == 2)
        return spence(1 - z);

    int k = 1;
    Type S = 0;
    Type delta = epsilon + 1;
    Type zk = z;
    while (delta > epsilon)
    {
        Type dS = zk / pow(k, n);
        k++;
        zk *= z;
        S += dS;
        delta = dS;
    }
    return S;
}

template <class NType, class EType>
void get_mu_l(NType N, EType E, double& mu, double& l,
              double epsilon=1e-8)
{
    mu = sqrt(polylog<double>(2, 1.) / double(E));
    l = 1. - exp(-double(N) * mu);

    double delta = epsilon + 1;
    while (delta > epsilon)
    {
        double nmu = sqrt(polylog<double>(2, l) / double(E));
        double nl = 1. - exp(-double(N) * mu);

        delta = std::abs(nmu - mu) + std::abs(nl - l);
        mu = nmu;
        l = nl;
    }

    l = -log(l);
}

template <class NType, class EType>
double get_xi_fast(NType N, EType E)
{
    if (E == 0 || N == 0)
        return 0;
    const double z2 = boost::math::zeta(2.);
    return 2 * sqrt(z2 * E);
}

template <class NType, class EType>
double get_xi(NType N, EType E, double epsilon=1e-8)
{
    if (E == 0 || N == 0)
        return 0;
    double mu = 0, l = 0;
    get_mu_l(N, E, mu, l, epsilon);
    double S = double(N) * l + 2 * double(E) * mu;
    return S;
}


#endif // UTIL_HH
