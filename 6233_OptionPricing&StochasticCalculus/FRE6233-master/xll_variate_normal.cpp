// xll_normal.cpp - Standard normal distribution
#include "fms_variate_normal.h"
#include "xll_FRE6233.h"

using namespace xll;
using namespace fms::variate;

AddIn xai_variate_normal(
	Function(XLL_HANDLEX, "xll_variate_normal", "\\VARIATE.NORMAL")
	.Uncalced()
	.Category(CATEGORY)
	.FunctionHelp("Return a handle to a standard normal model.")
	.Documentation(R"(
The standard normal random variate has density \(e^{-x^2/2}/\sqrt{2\pi}\), \(-\infty < x < \infty\).
)")
);
HANDLEX WINAPI xll_variate_normal()
{
#pragma XLLEXPORT
	handle<base> h(new normal{});

	return h.get();
}