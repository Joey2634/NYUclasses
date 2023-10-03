// xll_pwflat.cpp - Piecewise constant curves
#include "fms_pwflat.h"
#include "xll_FRE6233.h"

using namespace fms;
using namespace xll;

#ifdef _DEBUG
int fms_pwflat_value_test = pwflat::value_test();
int fms_pwflat_integral_test = pwflat::integral_test();
int fms_pwflat_discount_test = pwflat::discount_test();
int fms_pwflat_spot_test = pwflat::spot_test();
#endif // _DEBUG

AddIn xai_pwflat_curve_(
	Function(XLL_HANDLEX, "xll_pwflat_curve_", "\\PWF.CURVE")
	.Arguments({
		Arg(XLL_FP, "t", "is an array of times."),
		Arg(XLL_FP, "f", "is an array of values."),
		Arg(XLL_DOUBLE, "_f", "is an optional extrapolated value.")
		})
	.Uncalced()
	.Category(CATEGORY)
	.FunctionHelp("Return the value of a piecewise constant function.")
);
HANDLEX WINAPI xll_pwflat_curve_(const _FPX* pt, const _FPX* pf, double _f)
{
#pragma XLLEXPORT
	handle<pwflat::curve<>> h_(new pwflat::curve(size(*pt), pt->array, pf->array, _f));

	return h_.get();
}


AddIn xai_pwflat_curve_value(
	Function(XLL_DOUBLE, "xll_pwflat_curve_value", "PWF.CURVE.FORWARD")
	.Arguments({
		Arg(XLL_HANDLEX, "curve", "is a handle to a curve."),
		Arg(XLL_DOUBLE, "u", "is the value at which to calculate the forward."),
		})
	.Category(CATEGORY)
	.FunctionHelp("Return the value of a piecewise constant function.")
);
double WINAPI xll_pwflat_curve_value(HANDLEX h, double u)
{
#pragma XLLEXPORT
	handle<pwflat::curve<>> h_(h);

	return h_->forward(u);
}

AddIn xai_pwflat_value(
	Function(XLL_DOUBLE, "xll_pwflat_value", "PWF.VALUE")
	.Arguments({
		Arg(XLL_DOUBLE, "u", "is the value at which to calculate the function."),
		Arg(XLL_FP, "t", "is an array of times."),
		Arg(XLL_FP, "f", "is an array of values."),
		Arg(XLL_DOUBLE, "_f", "is an optional extrapolated value.")
		})
	.Category(CATEGORY)
	.FunctionHelp("Return the value of a piecewise constant function.")
);
double WINAPI xll_pwflat_value(double u, const _FPX* pt, const _FPX* pf, double _f)
{
#pragma XLLEXPORT
	return pwflat::value(u, size(*pt), pt->array, pf->array, _f);
}
