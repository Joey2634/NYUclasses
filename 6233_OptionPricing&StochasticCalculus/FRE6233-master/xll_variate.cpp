// xll_variate.cpp - Cumulative distribution, derivatives, and cumulant
#include "fms_variate.h"
#include "xll_FRE6233.h"

using namespace xll;
using namespace fms::variate;

AddIn xai_variate_cdf(
	Function(XLL_DOUBLE, "xll_variate_cdf", "VARIATE.CDF")
	.Arguments({
		Arg(XLL_HANDLEX, "model", "is a handle to a variate model."),
		Arg(XLL_DOUBLE, "x", "is the value at which to compute the CDF."),
		Arg(XLL_DOUBLE, "s", "is the Esscher parameter."),
		Arg(XLL_WORD, "nx", "is the number of derivatives with respect to x. Default is 0."),
		Arg(XLL_WORD, "ns", "is the number of derivatives with respect to s. Default is 0."),
		})
		.Category(CATEGORY)
	.FunctionHelp("Compute the Esscher transform of the cumulative distribution function.")
	.Documentation(R"(
Return the Esscher transform of the cumulative distribution function and its derivatives.
)")
);
double WINAPI xll_variate_cdf(HANDLEX v, double x, double s, WORD nx, WORD ns)
{
#pragma XLLEXPORT
	double result = XLL_NAN;

	try {
		handle<base> v_(v);
		ensure(v_);

		result = v_->cdf(x, s, nx, ns);
	}
	catch (const std::exception& ex) {
		XLL_ERROR(ex.what());
	}
	catch (...) {
		XLL_ERROR(__FUNCTION__ ": unknown exception");
	}

	return result;
}


AddIn xai_variate_cumulant(
	Function(XLL_DOUBLE, "xll_variate_cumulant", "VARIATE.CUMULANT")
	.Arguments({
		Arg(XLL_HANDLEX, "model", "is a handle to a variate model."),
		Arg(XLL_DOUBLE, "s", "is the Esscher parameter."),
		Arg(XLL_WORD, "ns", "is the number of derivatives with respect to s. Default is 0."),
		})
		.Category(CATEGORY)
	.FunctionHelp("Compute the cumulant of a variate.")
	.Documentation(R"(
Return the cumulant and derivatives of a variate.
)")
);
double WINAPI xll_variate_cumulant(HANDLEX v, double s, WORD ns)
{
#pragma XLLEXPORT
	double result = XLL_NAN;

	try {
		handle<base> v_(v);
		ensure(v_);

		result = v_->cumulant(s, ns);
	}
	catch (const std::exception& ex) {
		XLL_ERROR(ex.what());
	}
	catch (...) {
		XLL_ERROR(__FUNCTION__ ": unknown exception");
	}

	return result;
}