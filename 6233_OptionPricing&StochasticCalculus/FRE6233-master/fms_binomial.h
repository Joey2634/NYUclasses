// fms_binomial.h - binomial model
#pragma once
#include <cmath>
#include <functional>

// indicate error
#define ensure(e) if (!(e)) { return std::numeric_limits<double>::quiet_NaN(); }

namespace fms::binomial {

	// F_j = F e^{s W_j/sqrt(n)}/cosh(s/sqrt(n))^n
	// v_j(i) = E_j[nu(F_n) | F_j(i) = S, tau >= t_j]
	inline double value(int i, int j, int n, double f, double s, 
		const std::function<double(double)>& nu, bool american = false)
	{
		ensure(n > 0);

		double sn = s / sqrt(n);
		double F_j = f * exp(sn * (j - 2. * i)) / pow(cosh(sn), j);

		if (j == n) {
			return nu(F_j);
		}

		double v = (value(i, j + 1, n, f, s, nu, american)
			      + value(i + 1, j + 1, n, f, s, nu, american)) / 2;

		if (american) {
			//!!! optimal exercise code goes here !!!
		}

		return v;
	}

	// American put (k < 0) or call (p > 0) value at time j given W_j = i
	inline double value(int i, int j, int n, double f, double s, double k, bool american = false)
	{
		std::function<double(double)> nu;

		if (k < 0) { // put
			nu = [k](double F) { return std::max(-k - F, 0.); };
		}
		else if (k > 0) { // call
			nu = [k](double F) { return std::max(F - k, 0.); };
		}
		else {
			return std::numeric_limits<double>::quiet_NaN();
		}

		return value(i, j, n, f, s, nu, american);
	}

} // namespace fms

#undef ensure