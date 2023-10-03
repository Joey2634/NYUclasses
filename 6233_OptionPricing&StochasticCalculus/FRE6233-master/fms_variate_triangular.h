// fms_variate_triangular.h - Triangular variate
// A triangular distribution is defined by three values:
// l, m, h. Its density is a linear function connecting
// (l, 0), (m, 2/(h - l)), and (h, 0).
//
// f(x) = 0, x < l;
// f(x) = 2(x - l)/(m - l)(h - l), l <= x <= m;
// f(x) = 2(h - x)/(h - m)(h - l), m <= x <= h;
// f(x) = 0, x > h;
#pragma once
#include <algorithm>
#include <limits>
#include "fms_variate.h"

namespace fms::variate {

	// triangular distribution
	struct triangular : public fms::variate::base {
		double l, m, h;
		double a, b;
		
		static double I(double x, double s) {
			return exp(s * x) / s;
		};
		static double Ix(double x, double s) {
			return exp(s * x) * (x / s - 1 / (s * s));
		};

		triangular(double l, double m, double h)
			: l(l), m(m), h(h), a(2 / ((m - l) * (h - l))), b(2 / ((h - m) * (h - l)))
		{
			// ensure(l < m)
			// ensure(m < h);
		}
		// P^s(X <= x) = E[e^{s X - kappa(s)} 1(X <= x)] and derivatives
		// cdf(x, s, nx, ns) = int_{-infty^x} e^{s y - kappa(y)} f(y) dy.
		double _cdf(double x, double s, unsigned nx = 0, unsigned ns = 0) const override
		{
			double mgfs = mgf(s); // e^{kappa(s)}

			if (nx == 0 && ns == 0) {
				double Ps = 0;
				if (l <= x) {
					double xm = std::min(x, m);
					Ps = (a * Ix(xm, s) - a * l * I(xm, s)) - (a * Ix(l, s) - a * l * I(l, s));
					if (m <= x && x <= h) {
						Ps += (b * h * I(x, s) - b * Ix(x, s)) - (b * h * I(m, s) - b * Ix(m, s));
					}
				}

				return Ps / mgfs;
			}
			if (nx == 1 && ns == 0) {
				double ps = 0;
				if (l <= x && x <= m) {
					ps = a * (x - l);
				}
				else if (m <= x && x <= h) {
					ps = b * (h - x);
				}
				ps *= exp(s * x) / mgfs;

				return ps;
			}

			//!!! implement for nx = 0, ns = 0
			//!!! implement for nx = 1, ns = 0
			//!!! implement for nx = 1, ns = 1

			return std::numeric_limits<double>::quiet_NaN();
		}

		// kappa(s) = log E[e^{s X}] and derivatives
		// E[e^{s X}] =
		//     int_l^m e^{sx} a(x - l) dx
		//   + int_m^h e^{sx} b(h - x) dx
		//
		// int e^{sx} dx = e^{sx}/s
		// int x e^{sx} dx = e^{sx}(x/s - 1/s^2)
		//
		// E[e^{s X}] =
		//   [a e^{sx}(x/s - 1/s^2) - al e^{sx}/s]_l^m
		// + [bh e^{sx}/s - b e^{sx}(x/s - 1/s^2)]_m^h
		//
		double mgf(double s) const
		{
			if (s == 0) {
				return 1;
			}


			double Esx = (a * Ix(m, s) - a * l * I(m, s)) - (a * Ix(l, s) - a * l * I(l, s));
			Esx += (b * h * I(h, s) - b * Ix(h, s)) - (b * h * I(m, s) - b * Ix(m, s));

			return Esx;
		}
		double _cumulant(double s, unsigned n = 0) const override
		{
			if (n != 0) {
				return std::numeric_limits<double>::quiet_NaN();
			}

			return log(mgf(s));
		}
	};

} // namespace fms::variate
