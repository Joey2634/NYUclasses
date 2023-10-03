// fms_variate_normal.t.cpp - Test fms::variate::normal
#ifdef _DEBUG
// Only test in debug mode
#include <cassert>
#include <algorithm>
#include "fms_variate_normal.h"
#include "fms_derivative.h"

using namespace fms;
using namespace fms::variate;

int normal_H_test()
{
	double xs[] = { -1, 0, 1, 2 };
	{
		// H_0(x) = 1
		for (double x : xs) {
			assert(1 == normal::H(0, x));
		}

		// H_1(x) = x
		for (double x : xs) {
			assert(x == normal::H(1, x));
		}

		// H_2(x) = x^2 - 1
		for (double x : xs) {
			assert(x * x - 1 == normal::H(2, x));
		}
	}

	return 0;
}
// cause the test to be run when the dll is loaded
int normal_H_test_ = normal_H_test();

// test N^{(n)}(x)
template<class X = double, class Y = double>
inline bool normal_derivative_test(int n, X x, X h)
{
	Y df = normal::N(x, n + 1);
	Y dddf = normal::N(x, n + 3);
	auto f = [n](double x) { return normal::N(x, n); };

	return derivative_test<X, Y>(f, x, h, df, dddf);
}

int normal_test()
{
	{
		// sanity checks
		assert(0.5 == normal::N(0));
		assert(1 / M_SQRT2PI == normal::N(0, 1));
		assert(0 == normal::N(0, 2));
	}
	{
		double xs[] = { -1, 0, 1, 2 };
		double hs[] = { 0.1, 0.01, 0.001, 0.0001 };
		for (int n : { 0, 1, 2 }) {
			for (double x : xs) {
				for (double h : hs) {
					assert(normal_derivative_test(n, x, h));
				}
			}
		}
	}

	return 0;
}
int normal_test_ = normal_test();

// test d^nx/dx^nx d^ns/ds^ns cdf(x)
template<class X = double, class Y = double>
inline bool normal_cdf_derivative_test(int nx, int ns, X x, X h)
{
	normal N;

	auto f = [&N](double x) { return N.cdf(x); };
	Y df = N.cdf(x, h, nx + 1, ns);
	Y dddf = N.cdf(x, h, nx + 3, ns);

	return derivative_test<X, Y>(f, x, h, df, dddf);
}

int normal_cdf_test()
{
	{
		double xs[] = { -1, 0, 1, 2 };
		double hs[] = { 0.1, 0.01, 0.001, 0.0001 };
		for (int nx : { 0, 1, 2 }) {
			for (int ns : {0, 1, 2}) {
				for (double x : xs) {
					for (double h : hs) {
						assert(normal_cdf_derivative_test(nx, ns, x, h));
					}
				}
			}
		}
	}

	return 0;
}
int normal_cdf_test_ = normal_cdf_test();

#endif // _DEBUG