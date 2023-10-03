// xll_option.cpp - Black-Scholes/Merton option value and greeks.
#include "fms_option.h"
#include "fms_binomial.h"
#include "fms_variate_normal.h"
#include "xll_FRE6233.h"

using namespace xll;
using namespace fms;
using namespace fms::option;

static variate::normal N;

// return normal variate if v == 0
inline variate::base* pv(HANDLEX v)
{
	if (!v) {
		return &N;
	}

	handle<base> v_(v);
	ensure(v_);

	return v_.ptr();
}

AddIn xai_option_moneyness(
	Function(XLL_DOUBLE, "xll_option_moneyness", "OPTION.MONEYNESS")
	.Arguments({
		Arg(XLL_HANDLEX, "v", "is a handle to a variate."),
		Arg(XLL_DOUBLE, "f", "is the forward."),
		Arg(XLL_DOUBLE, "sigma", "is the volatility."),
		Arg(XLL_DOUBLE, "k", "is the strike."),
		Arg(XLL_DOUBLE, "t", "is the time in years to expiration."),
		})
	.FunctionHelp("Return the option moneyness.")
	.Category(CATEGORY)
	.Documentation(R"(
Moneyness is \((\log(k/f) + \kappa(s))/s\).
)")
);
double WINAPI xll_option_moneyness(HANDLEX v, double f, double sigma, double k, double t)
{
#pragma XLLEXPORT
	double result = XLL_NAN;

	try {
		result = moneyness(*pv(v), f, sigma * sqrt(t), fabs(k));
	}
	catch (const std::exception& ex) {
		XLL_ERROR(ex.what());
	}
	catch (...) {
		XLL_ERROR(__FUNCTION__ ": unknown exception");
	}

	return result;
}

#define OPTION_URL "https://en.wikipedia.org/wiki/Option_(finance)"

XLL_CONST(WORD, OPTION_PUT, contract::PUT, "European put option", CATEGORY, OPTION_URL);
XLL_CONST(WORD, OPTION_CALL, contract::CALL, "European call option", CATEGORY, OPTION_URL);
XLL_CONST(WORD, OPTION_DIGITAL_PUT, contract::DIGITAL_PUT, "European digital put option", CATEGORY, OPTION_URL);
XLL_CONST(WORD, OPTION_DIGITAL_CALL, contract::DIGITAL_CALL, "European digital call option", CATEGORY, OPTION_URL);

AddIn xai_option_value(
	Function(XLL_DOUBLE, "xll_option_value", "OPTION.VALUE")
	.Arguments({
		Arg(XLL_HANDLEX, "v", "is a handle to a variate."),
		Arg(XLL_DOUBLE, "S", "is the spot."),
		Arg(XLL_DOUBLE, "sigma", "is the volatility."),
		Arg(XLL_WORD, "option", "is the contract type from OPTION_*."),
		Arg(XLL_DOUBLE, "k", "is the strike."),
		Arg(XLL_DOUBLE, "t", "is the time in years to expiration."),
		Arg(XLL_DOUBLE, "_r", "is the optional continuously compouned interest rate. Default is 0."),
		Arg(XLL_LONG, "_n", "is the number of steps for binomial pricing. Default is 0.")
		})
	.FunctionHelp("Return the option value.")
	.Category(CATEGORY)
	.Documentation(R"(
The underlying at \(t\) is \(S_t = S e^{rt + \sigma\sqrt{t}X - \kappa(\sigma\sqrt{t})}\)
and \(X\) is a normalized variate determined by then handle <code>v</code>. 
Option value is 
\(e^{-rt}E[\phi(S_t)]\)
where \(\phi\) is the option payoff.
Note if \(r = 0\) this gives the Black value where
\(S\) is the forward.
If <code>n</code> is not 0 then a binomial model with <code>n</code> steps is used.
If <code>n &gt; 0</code> then the option is American.
If <code>n &lt; 0</code> then the option is European.)")
);
double WINAPI xll_option_value(HANDLEX v, double S, double sigma, contract flag, double k, double t, double r)
{
#pragma XLLEXPORT
	double result = XLL_NAN;

	try {
		result = bsm::value(*pv(v), r, S, sigma, flag, k, t);
	}
	catch (const std::exception& ex) {
		XLL_ERROR(ex.what());
	}
	catch (...) {
		XLL_ERROR(__FUNCTION__ ": unknown exception");
	}

	return result;
}

AddIn xai_option_delta(
	Function(XLL_DOUBLE, "xll_option_delta", "OPTION.DELTA")
	.Arguments({
		Arg(XLL_HANDLEX, "v", "is a handle to a variate."),
		Arg(XLL_DOUBLE, "S", "is the spot."),
		Arg(XLL_DOUBLE, "sigma", "is the volatility."),
		Arg(XLL_WORD, "option", "is the contract type from OPTION_*."),
		Arg(XLL_DOUBLE, "k", "is the strike."),
		Arg(XLL_DOUBLE, "t", "is the time in years to expiration."),
		Arg(XLL_DOUBLE, "r", "is the continuously compouned interest rate. Default is 0."),
		})
	.FunctionHelp("Return the option call or put delta.")
	.Category(CATEGORY)
	.Documentation(R"(
Option delta is the derivative of option value with respect to forward.
)")
);
double WINAPI xll_option_delta(HANDLEX v, double S, double sigma, contract flag, double k, double t, double r)
{
#pragma XLLEXPORT
	double result = XLL_NAN;

	try {
	    result = bsm::delta(*pv(v), r, S, sigma, flag, k, t);
	}
	catch (const std::exception& ex) {
		XLL_ERROR(ex.what());
	}
	catch (...) {
		XLL_ERROR(__FUNCTION__ ": unknown exception");
	}

	return result;
}

AddIn xai_option_gamma(
	Function(XLL_DOUBLE, "xll_option_gamma", "OPTION.GAMMA")
	.Arguments({
		Arg(XLL_HANDLEX, "v", "is a handle to a variate."),
		Arg(XLL_DOUBLE, "S", "is the spot."),
		Arg(XLL_DOUBLE, "sigma", "is the volatility."),
		Arg(XLL_WORD, "option", "is the contract type from OPTION_*."),
		Arg(XLL_DOUBLE, "k", "is the strike."),
		Arg(XLL_DOUBLE, "t", "is the time in years to expiration."),
		Arg(XLL_DOUBLE, "r", "is the continuously compouned interest rate. Default is 0."),
		})
		.FunctionHelp("Return the option call (k > 0) or put (k < 0) gamma.")
	.Category(CATEGORY)
	.Documentation(R"(
Option gamma is the second derivative of option value with respect to forward.
)")
);
double WINAPI xll_option_gamma(HANDLEX v, double S, double sigma, contract flag, double k, double t, double r)
{
#pragma XLLEXPORT
	double result = XLL_NAN;

	try {
	    result = bsm::gamma(*pv(v), r, S, sigma, flag, k, t);
	}
	catch (const std::exception& ex) {
		XLL_ERROR(ex.what());
	}
	catch (...) {
		XLL_ERROR(__FUNCTION__ ": unknown exception");
	}

	return result;
}

AddIn xai_option_vega(
	Function(XLL_DOUBLE, "xll_option_vega", "OPTION.VEGA")
	.Arguments({
		Arg(XLL_HANDLEX, "v", "is a handle to a variate."),
		Arg(XLL_DOUBLE, "S", "is the spot."),
		Arg(XLL_DOUBLE, "sigma", "is the volatility."),
		Arg(XLL_WORD, "option", "is the contract type from OPTION_*."),
		Arg(XLL_DOUBLE, "k", "is the strike."),
		Arg(XLL_DOUBLE, "t", "is the time in years to expiration."),
		Arg(XLL_DOUBLE, "r", "is the continuously compouned interest rate. Default is 0."),
		})
	.FunctionHelp("Return the option call (k > 0) or put (k < 0) vega.")
	.Category(CATEGORY)
	.Documentation(R"(
Option vega is the derivative of option value with respect to vol.
)")
);
double WINAPI xll_option_vega(HANDLEX v, double S, double sigma, contract flag, double k, double t, double r)
{
#pragma XLLEXPORT
	double result = XLL_NAN;

	try {
	    result = bsm::vega(*pv(v), r, S, sigma, flag, k, t);
	}
	catch (const std::exception& ex) {
		XLL_ERROR(ex.what());
	}
	catch (...) {
		XLL_ERROR(__FUNCTION__ ": unknown exception");
	}

	return result;
}

AddIn xai_option_theta(
	Function(XLL_DOUBLE, "xll_option_theta", "OPTION.THETA")
	.Arguments({
		Arg(XLL_HANDLEX, "v", "is a handle to a variate."),
		Arg(XLL_DOUBLE, "S", "is the spot."),
		Arg(XLL_DOUBLE, "sigma", "is the volatility."),
		Arg(XLL_WORD, "option", "is the contract type from OPTION_*."),
		Arg(XLL_DOUBLE, "k", "is the strike."),
		Arg(XLL_DOUBLE, "t", "is the time in years to expiration"),
		Arg(XLL_DOUBLE, "r", "is the continuously compouned interest rate. Default is 0."),
		Arg(XLL_DOUBLE, "_dt", "is an optional time increment. Default is 1/250."),
		})
	.FunctionHelp("Return the option call (k > 0) or put (k < 0) theta.")
	.Category(CATEGORY)
	.Documentation(R"(
Option theta is the derivative of option value with respect to time.
)")
);
double WINAPI xll_option_theta(HANDLEX v, double S, double sigma, int flag, double k, double t, double r, double dt)
{
#pragma XLLEXPORT
	double result = XLL_NAN;

	try {
		if (dt == 0) {
			dt = 1. / 250;
		}

	    result = bsm::theta(*pv(v), r, S, sigma, flag, k, t, dt);
	}
	catch (const std::exception& ex) {
		XLL_ERROR(ex.what());
	}
	catch (...) {
		XLL_ERROR(__FUNCTION__ ": unknown exception");
	}

	return result;
}

AddIn xai_option_implied(
	Function(XLL_DOUBLE, "xll_option_implied", "OPTION.IMPLIED")
	.Arguments({
		Arg(XLL_HANDLEX, "v", "is a handle to a variate."),
		Arg(XLL_DOUBLE, "f", "is the forward."),
		Arg(XLL_DOUBLE, "v", "is the option value."),
		Arg(XLL_DOUBLE, "k", "is the strike."),
		Arg(XLL_DOUBLE, "t", "is the time in years to expiration."),
		Arg(XLL_DOUBLE, "_sigma", "is an optional initial guess. Default is 0.1."),
		Arg(XLL_WORD, "_n", "is an optional maximum number of iterations. Default is 100."),
		Arg(XLL_DOUBLE, "_tol", "is an optional absolute tolerance. Default is square root of machine epsilon."),
		})
	.FunctionHelp("Return the option call (k > 0) or put (k < 0) implied vol.")
	.Category(CATEGORY)
	.Documentation(R"(
Option implied volatility is the inverse of value.
)")
);
double WINAPI xll_option_implied(HANDLEX v, double f, double v0, double k, double t, double sigma, unsigned n, double tol)
{
#pragma XLLEXPORT
	double result = XLL_NAN;

	try {
	    result = black::implied(*pv(v), f, v0, k, sigma * sqrt(t), n, tol);
	}
	catch (const std::exception& ex) {
		XLL_ERROR(ex.what());
	}
	catch (...) {
		XLL_ERROR(__FUNCTION__ ": unknown exception");
	}

	return result;
}

AddIn xai_option_variance(
	Function(XLL_DOUBLE, "xll_option_variance", "OPTION.VARIANCE")
	.Arguments({
		Arg(XLL_HANDLEX, "v", "is a handle to a variate."),
		Arg(XLL_DOUBLE, "S", "is the spot."),
		Arg(XLL_DOUBLE, "sigma", "is the volatility."),
		Arg(XLL_WORD, "option", "is the contract type from OPTION_*."),
		Arg(XLL_DOUBLE, "k", "is the strike."),
		Arg(XLL_DOUBLE, "t", "is the time in years to expiration."),
		Arg(XLL_DOUBLE, "r", "is the continuously compouned interest rate. Default is 0."),
		})
		.FunctionHelp("Return the option variance.")
	.Category(CATEGORY)
	.Documentation(R"(
Option variance.
)")
);
double WINAPI xll_option_variance(HANDLEX v, double S, double sigma, contract flag, double k, double t, double r)
{
#pragma XLLEXPORT
	double result = XLL_NAN;

	try {
		result = bsm::variance(*pv(v), r, S, sigma, flag, k, t);
	}
	catch (const std::exception& ex) {
		XLL_ERROR(ex.what());
	}
	catch (...) {
		XLL_ERROR(__FUNCTION__ ": unknown exception");
	}

	return result;
}
