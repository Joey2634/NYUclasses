// fms_binomial.t.cpp - get binomial
#include <cassert>
#include "fms_binomial.h"

using namespace fms;

int binomial_test(int n, double f, double s, double k, double tol = 1e-13)
{
	{
		double c = binomial::value(0, 0, n, f, s, k);
		double p = binomial::value(0, 0, n, f, s, -k);
		double err = (c - p) - (f - k);
		assert(fabs(err) < tol);

		double ap = binomial::value(0, 0, n, f, s, -k, true);
		assert(ap >= p);
		double ac = binomial::value(0, 0, n, f, s, k, true);
		assert(fabs(ac - c) < tol);
	}

	return 0;
}

int binomial_tests()
{
	binomial_test(10, 100, .2, 100);
	binomial_test(10, 100, .1, 90);
	binomial_test(10, 110, .1, 100);

	return 0;
}
int binomial_tests_ = binomial_tests();