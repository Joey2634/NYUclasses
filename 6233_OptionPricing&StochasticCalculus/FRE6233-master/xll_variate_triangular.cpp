// xll_triangular.cpp - Standard triangular distribution
#include "fms_variate_triangular.h"
#include "xll_FRE6233.h"

using namespace xll;
using namespace fms::variate;

AddIn xai_variate_triangular(
	Function(XLL_HANDLEX, "xll_variate_triangular", "\\VARIATE.TRIANGULAR")
	.Arguments({
		Arg(XLL_DOUBLE, "l", "is the low."),
		Arg(XLL_DOUBLE, "m", "is the middle."),
		Arg(XLL_DOUBLE, "h", "is the high."),
		})
	.Uncalced()
	.Category(CATEGORY)
	.FunctionHelp("Return a handle to a standard triangular model.")
	.Documentation(R"(
The standard triangular random variate has density 
!!!.
)")
);
HANDLEX WINAPI xll_variate_triangular(double l, double m, double h)
{
#pragma XLLEXPORT
	handle<base> t(new triangular(l, m, h));

	return t.get();
}