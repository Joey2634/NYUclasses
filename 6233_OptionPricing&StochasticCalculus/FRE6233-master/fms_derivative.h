// fms_derivative.h - Test derivatives
#pragma once
#include <functional>
#include <limits>
#include <valarray>
namespace fms {

	inline constexpr double epsilon = std::numeric_limits<double>::epsilon();

	inline double relative_difference(double a, double b)
	{
		if (a == 0 && b == 0) {
			return 0;
		}

		if (std::signbit(a) != std::signbit(b)) {
			return std::numeric_limits<double>::max();
		}

		if (a == 0) {
			a = epsilon;
		}
		if (b == 0) {
			b = epsilon;
		}

		return std::fabs(std::fabs(a - b)/std::min(a,b));
	}
	inline double epsilon_difference(double a, double b)
	{
		return relative_difference(a,b)/epsilon;
	}

	// sequence from b to e with increment dx
	inline std::valarray<double> sequence(double b, double e, double dx)
	{
		size_t n = static_cast<unsigned>(1 + (e - b) / dx);
		std::valarray<double> a(n);

		for (size_t i = 0; i < n; ++i) {
			a[i] = b + i * dx;
		}

		return a;
	}

	// symmetric difference quotient
	template<class X, class Y>
	inline auto difference_quotient(const std::function<Y(X)>& f, X dx)
	{
		return [dx, &f](X x) { 
			return (f(x + dx) - f(x - dx)) / (2 * dx); 
		};
	}

	// (f(x + h) - f(x - h))/2h = f'(x) + f'''(x) h^2/3! + O(h^3)
	template<class X, class Y>
	inline bool derivative_test(const std::function<Y(X)>& f, X x, X h, X df, X dddf, X O = 1)
	{
		dddf = std::max(fabs(dddf), 1.);
		auto Df = difference_quotient(f, h);

		double lhs = Df(x) - df;
		double rhs = dddf * h * h / 6;
		
		bool b = fabs(lhs) <= O * rhs;
		if (!b) {
			O = fabs(lhs)/fabs(rhs); // set breakpoint to check tolerance
		}

		return b;
	}
}
