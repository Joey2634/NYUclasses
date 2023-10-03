// fsm_option.h - Option value and greeks
#pragma once
#include <cmath>
#include <limits>
#include <tuple>
#include "fms_variate.h"

using namespace fms::variate;

namespace fms {

	// Return NaN to indicate error.
	constexpr double NaN = std::numeric_limits<double>::quiet_NaN();

	namespace option {
		enum contract {
			PUT = 'P',
			DIGITAL_PUT = 'Q',
			CALL = 'C',
			DIGITAL_CALL = 'D',
		};

		//  moneyness
		inline double moneyness(const variate::base& v, double f, double s, double k)
		{
			if (f <= 0 || s <= 0 || k <= 0) {
				return NaN;
			}

			return (log(k / f) + v.cumulant(s)) / s;
		}

		// E[(F/f)^n 1(F <= k)] = e^{kappa(ns) - n kappa(s)} P_{ns}(X <= x)
		// E[(F/f)^n 1(F > k)] = e^{kappa(ns) - n kappa(s)} P_{ns}(X > x)
		inline double partial_moment(const variate::base& v, double f, double s, double k, int n)
		{
			double x = moneyness(v, f, s, fabs(k));

			if (n == 0) {
				return k < 0 ? v.cdf(x) : 1 - v.cdf(x);
			}

			if (n == 1) {
				return k < 0 ? v.cdf(x, s) : 1 - v.cdf(x,s);
			}

			double Kn = v.cumulant(n * s) - n * v.cumulant(s);

			return exp(Kn) * (k < 0 ? v.cdf(x, n * s) : 1 - v.cdf(x, n*s));
		}

		// Use 0 rate and forward values
		namespace black {
			// put (k < 0) or call (k > 0) option value
			inline double value(const variate::base& v, double f, double s, double k)
			{
				if (k < 0) { // put
					double x = moneyness(v, f, s, -k);

					return (-k) * v.cdf(x, 0) - f * v.cdf(x, s);
				}
				else if (k > 0) { // call
					// c = p + f - k
					return value(v, f, s, -k) + f - k;
				}

				// k = -/+ 0
				return signbit(k) ? 0 : f;
			}

			// put (k < 0) or call (k > 0) option delta, dv/df
			inline double delta(const variate::base& v, double f, double s, double k)
			{
				if (k < 0) { // put
					double x = moneyness(v, f, s, -k);

					return -v.cdf(x, s, 0, 0);
				}
				else if (k > 0) { // call
					// dc/df = dp/df + 1
					return delta(v, f, s, -k) + 1;
				}

				return signbit(k) ? 0 : 1;
			}

			// put (k < 0) or call (k > 0) option gamma, d^2v/df^2
			inline double gamma(const variate::base& v, double f, double s, double k)
			{
				double x = moneyness(v, f, s, std::fabs(k));

				return v.cdf(x, s, 1, 0) / (f * s);
			}

			// n-th derivative with respect to f
			inline double value(const variate::base& v, double f, double s, double k, unsigned n)
			{
				if (n == 0) {
					return value(v, f, s, k);
				}
				if (n == 1) {
					return delta(v, f, s, k);
				}

				double x = moneyness(v, f, s, std::fabs(k));

				return v.cdf(x, s, n - 1, 0) / pow(f * s, n - 1);
			}

			// put (k < 0) or call (k > 0) option vega, dv/ds
			inline double vega(const variate::base& v, double f, double s, double k)
			{
				double x = moneyness(v, f, s, fabs(k));

				return -f * v.cdf(x, s, 0, 1);
			}

			// put (k < 0) or call (k > 0) option theta, -dv/dt
			inline double theta(const variate::base& v, double f, double sigma, double k, double t, double dt = 1. / 250)
			{
				double s = sigma * sqrt(t);
				double v0 = value(v, f, s, k);
				s = sigma * sqrt(t - dt);
				double v_ = value(v, f, s, k);

				return (v_ - v0) / dt;
			}

			// implied volatility using initial guess, max number of iterations, and tolerance
			inline double implied(const variate::base& v, double f, double v0, double k,
				double s = 0, unsigned n = 0, double tol = 0)
			{
				// max(k - f,0) >= k - f
				// max(k - f,0) <= k
				if (k < 0) {
					if (v0 <= std::max(-k - f, 0.) || v0 >= k) {
						return NaN;
					}
				}
				// max(f - k,0) >= f - k
				// max(f - k,0) <= f
				else if (k > 0) {
					if (v0 <= std::max(f - k, 0.) || v0 >= f) {
						return NaN;
					}
				}

				if (s == 0) {
					s = 0.1; // initial vol guess
				}
				if (n == 0) {
					n = 100; // maximum number of iterations
				}
				if (tol == 0) {
					tol = sqrt(std::numeric_limits<double>::epsilon()); // absolute tolerance
				}

				double v_ = value(v, f, s, k);
				double dv_ = vega(v, f, s, k); // dv/ds
				double s_ = s - (v_ - v0) / dv_; // Newton-Raphson
				if (s_ < 0) {
					s_ = s / 2;
				}
				while (fabs(s_ - s) > tol) {
					v_ = value(v, f, s_, k);
					dv_ = vega(v, f, s_, k);
					s = s_ - (v_ - v0) / dv_;
					if (s < 0) {
						s = s_ / 2;
					}
					std::swap(s_, s);
					if (n == 0) {
						return NaN;
					}
					--n;
				}

				return s_;
			}

			// Var((k - F)^+) = E[(k - F)^2 1(F <= k)] - E[(k - F) 1(F <= k)]^2
			inline double variance(const variate::base& v, double f, double s, double k)
			{
				double P0 = partial_moment(v, f, s, k, 0);
				double P1 = partial_moment(v, f, s, k, 1);
				double P2 = partial_moment(v, f, s, k, 2);

				if (k < 0) {
					double o = -k * P0 - f * P1;
					double o2 = k * k * P0 + 2 * k * f * P1 + f * f * P2;

					return o2 - o * o;
				}
				else if (k > 0) {
					double o = f * P1 - k * P0;
					double o2 = f * f * P2 - 2 * k * f * P1 + k * k * P0;

					return o2 - o * o;
				}

				return signbit(k) ? 0 : f * f * (exp(v.cumulant(2 * s) - 2 * v.cumulant(s)) - 1);
			}

			inline double moment4(const variate::base& v, double f, double s, double k)
			{
				double P0 = partial_moment(v, f, s, k, 0);
				double P1 = partial_moment(v, f, s, k, 1);
				double P2 = partial_moment(v, f, s, k, 2);
				double P3 = partial_moment(v, f, s, k, 3);
				double P4 = partial_moment(v, f, s, k, 4);
				if (k < 0) {
					double o = -k * P0 - f * P1;
					double o2 = k * k * P0 + 2 * k * f * P1 + f * f * P2;
					double o3 = -k * k * k * P0 - 3 * k * k * f * P1 - 3 * k * f * f * P2 - f * f * f * P3;
					double o4 = k * k * k * k * P0 + 4 * k * k * k * f * P1 + 6 * k * k * f * f * P2 + 4 * k * f * f * f * P3 + f * f * f * f * P4;
					return o4 - 4 * o3 * o + 6 * o2 * o * o - 3 * o * o * o * o;
				}
				else if (k > 0) {
					double o = f * P1 - k * P0;
					double o2 = f * f * P2 - 2 * k * f * P1 + k * k * P0;
					double o3 = f * f * f * P3 - 3 * f * f * k * P2 + 3 * f * k * k * P1 - k * k * k * P0;
					double o4 = f * f * f * f * P4 - 4 * f * f * f * k * P3 + 6 * f * f * k * k * P2 - 4 * f * k * k * k * P1 + k * k * k * k * P0;
					return o4 - 4 * o3 * o + 6 * o2 * o * o - 3 * o * o * o * o;
				}

				return signbit(k) ? 0 : f * f * (exp(v.cumulant(2 * s) - 2 * v.cumulant(s)) - 1);
			}
		}

		namespace digital {

			// q = P(F <= -k), k < 0, or d = P(F > k), k > 0
			inline double value(const variate::base& v, double f, double s, double k)
			{
				double x = moneyness(v, f, s, fabs(k));
				double v0 = v.cdf(x, 0);

				if (k < 0) {
					return v0;
				}
				else if (k > 0) {
					return 1 - v0;
				}

				return signbit(k) ? 0 : 1;
			}
			// dq/df or dd/df
			inline double delta(const variate::base& v, double f, double s, double k)
			{
				double x = moneyness(v, f, s, fabs(k));
				double v0 = -v.cdf(x, 0, 1) / (f * s);

				if (k < 0) {
					return v0;
				}
				else if (k > 0) {
					return -v0;
				}

				return 0;
			}
			// d^2q/df^2 or d^d/df^2
			inline double gamma(const variate::base& v, double f, double s, double k)
			{
				double x = moneyness(v, f, s, fabs(k));
				double v0 = (s * v.cdf(x, 0, 1) + v.cdf(x, 0, 2)) / (f * f * s * s);

				if (k < 0) {
					return v0;
				}
				else if (k > 0) {
					return -v0;
				}

				return 0;
			}
			// dq/ds or dd/ds
			inline double vega(const variate::base& v, double f, double s, double k)
			{
				double x = moneyness(v, f, s, fabs(k));
				double v0 = v.cdf(x, 0, 1) * (v.cumulant(s, 1) - x) / s;

				if (k < 0) {
					return v0;
				}
				else if (k > 0) {
					return -v0;
				}

				return 0;
			}
		}

		namespace bsm {

			// Convert B-S/M parameters to Black forward parameters.
			inline auto Dfs(double r, double S, double sigma, double t)
			{
				if (t == 0) {
					t = 1;
				}
				double D = exp(-r * t);
				double f = S / D;
				double s = sigma * sqrt(t);

				return std::tuple(D, f, s);
			}

			inline double moneyness(const variate::base& v, double r, double S, double sigma, double k, double t)
			{
				auto [D, f, s] = Dfs(r, S, sigma, t);

				return option::moneyness(v, f, s, fabs(k));
			}

			inline double value(const variate::base& v, double r, double S, double sigma, int c, double k, double t)
			{
				auto [D, f, s] = Dfs(r, S, sigma, t);

				switch (c) {
				case option::contract::PUT:
					return D * black::value(v, f, s, -k);
				case option::contract::CALL:
					return D * black::value(v, f, s, k);
				case option::contract::DIGITAL_PUT:
					return D * digital::value(v, f, s, -k);
				case option::contract::DIGITAL_CALL:
					return D * digital::value(v, f, s, k);
				}

				return NaN;
			}

			// delta
			inline double delta(const variate::base& v, double r, double S, double sigma, int c, double k, double t)
			{
				auto [D, f, s] = Dfs(r, S, sigma, t);

				switch (c) {
				case option::contract::PUT:
					return black::delta(v, f, s, -k);
				case option::contract::CALL:
					return black::delta(v, f, s, k);
				case option::contract::DIGITAL_PUT:
					return digital::delta(v, f, s, -k);
				case option::contract::DIGITAL_CALL:
					return digital::delta(v, f, s, k);
				}

				return NaN;
			}
			// gamma
			inline double gamma(const variate::base& v, double r, double S, double sigma, int c, double k, double t)
			{
				auto [D, f, s] = Dfs(r, S, sigma, t);

				switch (c) {
				case option::contract::PUT:
					return black::gamma(v, f, s, -k) / D;
				case option::contract::CALL:
					return black::gamma(v, f, s, k) / D;
				case option::contract::DIGITAL_PUT:
					return digital::gamma(v, f, s, -k) / D;
				case option::contract::DIGITAL_CALL:
					return digital::gamma(v, f, s, k) / D;
				}

				return NaN;
			}
			// vega
			inline double vega(const variate::base& v, double r, double S, double sigma, int c, double k, double t)
			{
				auto [D, f, s] = Dfs(r, S, sigma, t);

				switch (c) {
				case option::contract::PUT:
					return black::vega(v, f, s, -k) * D * sqrt(t);
				case option::contract::CALL:
					return black::vega(v, f, s, k) * D * sqrt(t);
				case option::contract::DIGITAL_PUT:
					return digital::vega(v, f, s, -k) * D * sqrt(t);
				case option::contract::DIGITAL_CALL:
					return digital::vega(v, f, s, k) * D * sqrt(t);
				}

				return NaN;
			}

			// put (k < 0) or call (k > 0) option theta, -dv/dt
			inline double theta(const variate::base& v, double r, double S, double sigma, int c, double k, double t, double dt = 1. / 250)
			{
				double v0 = value(v, r, S, sigma, c, k, t);
				double v_ = value(v, r, S, sigma, c, k, t - dt);

				return (v_ - v0) / dt;
			}
			// variance
			inline double variance(const variate::base& v, double r, double S, double sigma, int c, double k, double t)
			{
				auto [D, f, s] = Dfs(r, S, sigma, t);

				switch (c) {
				case option::contract::PUT:
					return D * D * black::variance(v, f, s, -k);
				case option::contract::CALL:
					return D * D * black::variance(v, f, s, k);
				}

				return NaN;
			}

		} // namespace bsm
	}
}
