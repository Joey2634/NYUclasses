// fms_variate_normal.h - Normally distributed random variate
#pragma once
#include <cmath>
#include "fms_variate.h"

namespace fms::variate {

	// sqrt(2 pi)
	constexpr double M_SQRT2PI = 2.50662827463100050240;
#ifndef M_SQRT2
	// sqrt(2)
	constexpr double M_SQRT2 = 1.41421356237309504880;
#endif

	// standard normal random variate
	struct normal : public fms::variate::base {

		// Hermite polynomials
		// H_{n+1}(x) = x H_n(x) - n H_{n-1}(x), H_0(x) = 1, H_1(x) = x
		static constexpr double H(unsigned n, double x) noexcept
		{
			if (n == 0) {
				return 1;
			}
			if (n == 1) {
				return x;
			}

			return x * H(n - 1, x) - (n - 1) * H(n - 2, x);
		}

		// P(X <= x) and derivatives
		static double N(double x, unsigned n = 0)
		{
			if (n == 0) {
				return  (1 + erf(x / M_SQRT2)) / 2;
			}

			double phi = exp(-x * x / 2) / M_SQRT2PI;

			if (n == 1) {
				return phi;
			}

			// (d/dx)^n phi(x) = (-1)^n phi(x) H_n(x)
			return phi * H(n - 1, x) * ((n & 1) ? 1 : -1);
		}

		// P^s(X <= x) = P(X <= x - s) and derivatives
		double _cdf(double x, double s, unsigned nx = 0, unsigned ns = 0) const override
		{
			return N(x - s, nx + ns) * (ns & 1 ? -1 : 1);
		}

		// kappa(s) = log E[e^{s X}] = s^2/2 and derivativs
		double _cumulant(double s, unsigned n = 0) const override
		{
			if (n == 0) {
				return s * s / 2;
			}
			if (n == 1) {
				return s;
			}
			if (n == 2) {
				return 1;
			}

			return 0;
		}
	};

} // namespace fms::variate
