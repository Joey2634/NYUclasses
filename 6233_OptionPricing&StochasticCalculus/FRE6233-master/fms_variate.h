// fms_variate.h - NVI base class for all variate bases
#pragma once

namespace fms::variate {

	struct base {
		virtual ~base()
		{ }

		// P^s(X <= x) = E[e^{s X - kappa(s)} 1(X <= x)] and derivatives wrt x and s
		double cdf(double x, double s = 0, unsigned nx = 0, unsigned ns = 0) const
		{
			return _cdf(x, s, nx, ns);
		}
		// kappa(s) = log E[exp(s X)] 
		double cumulant(double s, unsigned n = 0) const
		{
			return _cumulant(s, n);
		}
	private:
		// overridden in derived class
		virtual double _cdf(double x, double s, unsigned nx, unsigned ns) const = 0;
		virtual double _cumulant(double s, unsigned n) const = 0;
	};

} // namespace fms