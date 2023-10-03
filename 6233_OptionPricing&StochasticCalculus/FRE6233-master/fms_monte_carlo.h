// fms_monte_carlo.h - Monte Carlo simulation
#pragma once
#include <cmath>
#include <type_traits>

namespace fms::monte_carlo {

	// Numerically stable average.
	// s_n = (x_1 + ... + x_n)/n
	// n s_n - (n-1) s_{n-1} = x_n
	// s_n = s_{n-1} + (x_n - s_{n-1})/n
	template<class X, class S = std::invoke_result<X>::type>
	inline S average(size_t n, X& x)
	{
		S s = 0;

		for (size_t m = 1; m <= n; ++m) {
			s += (x() - s) / m;
		}

		return s;
	}

	//Welford's online algo, compute the var in one pass
	template<class X, class S = std::invoke_result<X>::type>
	inline S stddev(size_t n, X& x) 
	{
		S s = 0;
		S avglast = 0;
		S avg = x();
		for (unsigned int m = 2; m <= n; ++m) {
			S temp = x();
			avglast = avg;
			avg += (temp - s) / m;
			s += ((temp - avglast) * (temp - avg) - s) / m;
		}

		return std::sqrt(s);
	}
}
