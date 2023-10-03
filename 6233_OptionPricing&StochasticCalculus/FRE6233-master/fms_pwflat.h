// fms_pwflat.h - Piecewise flat left-continuous curve
#pragma once
#ifdef _DEBUG
#include <cassert>
#endif
#include <algorithm>
#include <limits>
#include <vector>
/*
		   { f[i] if t[i-1] < t <= t[i];
	f(t) = { _f   if t > t[n-1];
		   { NaN  if t < 0
	|                                   _f
	|        f[1]             f[n-1] o--------
	| f[0] o------          o--------x
	x------x      ... ------x
	|
	0-----t[0]--- ... ---t[n-2]---t[n-1]
*/
namespace fms::pwflat {

	template<class X>
	inline static constexpr X NaN = std::numeric_limits<X>::quiet_NaN();

	// f(u) = f[i], t[i-1] < u <= t[i], 0 <= i < n
	// f(u) = _f, u > t[n-1]
	// f(u) = NaN, u < 0
	template<class T = double, class F = double>
	inline F value(T u, size_t n, const T* t, const F* f,
		F _f = NaN<F>)
	{
		if (u < 0) {
			return NaN<F>;
		}

		// first element in t >= u
		auto i = std::lower_bound(t, t + n, u);

		return i == t + n ? _f : f[i - t];
	}

#ifdef _DEBUG
	inline int value_test()
	{
		double t[] = { 1,2,3 };
		double f[] = { .1,.2,.3 };
		assert(_isnan(value(-1., 3, t, f)));
		assert(f[0] == value(0., 3, t, f));
		double t0 = 0;
		for (int i : {0, 1, 2}) {
			assert(f[i] == value(t[i], 3, t, f));
			double ti = (t[i] + t0) / 2;
			assert(f[i] == value(ti, 3, t, f));
			t0 = t[i];
		}
		assert(4 == value(t[2] + 1, 3, t, f, 4.));

		return 0;
	}
#endif // _DEBUG

	// int_0^u f(t) dt
	template<class T, class F>
	inline F integral(T u, size_t n, const T* t, const F* f, F _f = NaN<F>)
	{
		if (u < 0)
			return NaN<F>;

		F I = 0;
		T t_ = 0;

		size_t i;
		for (i = 0; i < n && t[i] <= u; ++i) {
			I += f[i] * (t[i] - t_);
			t_ = t[i];
		}
		I += (n == 0 || u > t[n - 1] ? _f : f[i]) * (u - t_);

		return I;
	}

#ifdef _DEBUG
	inline int integral_test()
	{
		double t[] = { 1,2,3 };
		double f[] = { .1,.2,.3 };

		// test integral(u, n, t, f) for u = -0.5, 0, 0.5, 1, ..., 3.5
		// by comparing hand computation.
		// E.g., assert(integral(0.5, ...) == .1*0.5);
		assert(_isnan(integral(-0.5, 3, t, f)));
		assert(integral(0.0, 3, t, f) == f[0] * 0.0);
		assert(integral(0.5, 3, t, f) == f[0] * 0.5);
		assert(integral(1.0, 3, t, f) == f[0] * 1.0);
		assert(integral(1.5, 3, t, f) == f[0] * 1.0 + f[1] * 0.5);
		assert(integral(2.0, 3, t, f) == f[0] * 1.0 + f[1] * 1.0);
		assert(integral(2.5, 3, t, f) == f[0] * 1.0 + f[1] * 1.0 + f[2] * 0.5);
		assert(integral(3.0, 3, t, f) == f[0] * 1.0 + f[1] * 1.0 + f[2] * 1.0);
		assert(integral(3.5, 3, t, f, 4.0) == f[0] * 1.0 + f[1] * 1.0 + f[2] * 1.0 + 4.0 * 0.5);
		
		return 0;
	}
#endif // _DEBUG

	// D(u) = exp(-int_0^u f(t) dt)
	template<class T, class F>
	inline F discount(T u, size_t n, const T* t, const F* f, F _f = NaN<F>)
	{
		return exp(-integral(u, n, t, f, _f));
	}

#ifdef _DEBUG

	inline int discount_test() {
		// !!! add tests
		return 0;
	}

#endif // _DEBUG

	// r(u) = (1/u) int_0^u f(t) dt
	template<class T, class F>
	inline F spot(T u, size_t n, const T* t, const F* f, F _f = NaN<F>)
	{
		if (u < 0) {
			return NaN<F>;
		}
		if (n == 0) {
			return _f;
		}

		return u <= t[0] ? f[0] : integral(u, n, t, f, _f) / u;
	}
#ifdef _DEBUG

	inline int spot_test()
	{
		// !!! add tests
		return 0;
	}

#endif // _DEBUG

	// pwflat forward value type
	template<class T = double, class F = double>
	class curve {
		std::vector<T> t;
		std::vector<F> f;
		F _f;
	public:
		// default constructable
		curve(F _f = NaN<F>)
			: _f(_f)
		{ }
		curve(size_t n, const T* t, const F* f, F _f = NaN<F>)
			: t(t, t + n), f(f, f + n), _f(_f)
		{ }
		curve(const curve&) = default;
		curve& operator=(const curve&) = default;
		~curve()
		{ }

		curve& extend(T t_, F f_)
		{
			// ensure(t_ > t.back());
			t.push_back(t_);
			f.push_back(f_);

			return *this;
		}

		F forward(T u) const
		{
			return value(u, t.size(), t.data(), f.data(), _f);
		}
		F integral(T u) const
		{
			return integral(u, t.size(), t.data(), f.data(), _f);
		}
		F spot(T u) const
		{
			return spot(u, t.size(), t.data(), f.data(), _f);
		}
		F discount(T u) const
		{
			return discount(u, t.size(), t.data(), f.data(), _f);
		}

#ifdef _DEBUG
		int test()
		{
			{
				curve c;
			}

			return 0;
		}
#endif // _DEBUG

	};

} // namespace fms
