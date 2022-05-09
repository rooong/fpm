#ifndef FPM_FIXED_HPP
#define FPM_FIXED_HPP

#include <cassert>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <type_traits>
#include <string>

namespace fpm
{

//! Fixed-point number type
//! \tparam BaseType         the base integer type used to store the fixed-point number. This can be a signed or unsigned type.
//! \tparam IntermediateType the integer type used to store intermediate results during calculations.
//! \tparam FractionBits     the number of bits of the BaseType used to store the fraction
template <typename BaseType, typename IntermediateType, unsigned int FractionBits>
class fixed
{
    static_assert(std::is_integral<BaseType>::value, "BaseType must be an integral type");
    static_assert(FractionBits > 0, "FractionBits must be greater than zero");
    static_assert(FractionBits <= sizeof(BaseType) * 8, "BaseType must at least be able to contain entire fraction");
    static_assert(FractionBits <= 62, "Fraction may be no more than 62 bits");
    static_assert(sizeof(IntermediateType) > sizeof(BaseType), "IntermediateType must be larger than BaseType");
    static_assert(std::is_signed<IntermediateType>::value == std::is_signed<BaseType>::value, "IntermediateType must have same signedness as BaseType");

    static const BaseType FRACTION_MULT = BaseType(1) << FractionBits;

    struct raw_construct_tag {};
    inline fixed(BaseType val, raw_construct_tag) : m_value(val) {}

public:
    inline fixed() = default;

    // Converts an integral number to the fixed-point type.
    // Like static_cast, this truncates bits that don't fit.
    template <typename T, typename std::enable_if<std::is_integral<T>::value>::type* = nullptr>
    inline fixed(T val)
        : m_value(static_cast<BaseType>(val * FRACTION_MULT))
    {}

    // Converts an floating-point number to the fixed-point type.
    // Like static_cast, this truncates bits that don't fit.
    template <typename T, typename std::enable_if<std::is_floating_point<T>::value>::type* = nullptr>
    inline fixed(T val)
        : m_value(static_cast<BaseType>((val >= 0.0) ? (val * FRACTION_MULT + T{0.5}) : (val * FRACTION_MULT - T{0.5})))
    {}

    // Constructs from another fixed-point type with possibly different underlying representation.
    // Like static_cast, this truncates bits that don't fit.
    template <typename B, typename I, unsigned int F>
    inline explicit fixed(fixed<B,I,F> val)
        : m_value(from_fixed_point<F>(val.raw_value()).raw_value())
    {}

    // Explicit conversion to a floating-point type
    template <typename T, typename std::enable_if<std::is_floating_point<T>::value>::type* = nullptr>
    inline explicit operator T()
    {
        return static_cast<T>(m_value) / FRACTION_MULT;
    }

    // Explicit conversion to an integral type
    template <typename T, typename std::enable_if<std::is_integral<T>::value>::type* = nullptr>
    inline explicit operator T()
    {
        return static_cast<T>(m_value / FRACTION_MULT);
    }

    // Returns the raw underlying value of this type.
    // Do not use this unless you know what you're doing.
    const inline BaseType raw_value() const
    {
        return m_value;
    }

    const inline float toFloat() const
    {
        return static_cast<float>(m_value) / FRACTION_MULT;
    }

    const inline int toInt() const
    {
        return static_cast<int>(m_value / FRACTION_MULT);
    }

    //! Constructs a fixed-point number from another fixed-point number.
    //! \tparam NumFractionBits the number of bits used by the fraction in \a value.
    //! \param value the integer fixed-point number
    template <unsigned int NumFractionBits, typename T, typename std::enable_if<(NumFractionBits > FractionBits)>::type* = nullptr>
    static const inline fixed from_fixed_point(T value)
    {
        // To correctly round the last bit in the result, we need one more bit of information.
        // We do this by multiplying by two before dividing and adding the LSB to the real result.
        return fixed(static_cast<BaseType>(
             value / (T(1) << (NumFractionBits - FractionBits)) +
            (value / (T(1) << (NumFractionBits - FractionBits - 1)) % 2)),
            raw_construct_tag{});
    }

    template <unsigned int NumFractionBits, typename T, typename std::enable_if<(NumFractionBits <= FractionBits)>::type* = nullptr>
    static const inline fixed from_fixed_point(T value)
    {
        return fixed(static_cast<BaseType>(
            value * (T(1) << (FractionBits - NumFractionBits))),
            raw_construct_tag{});
    }

    // Constructs a fixed-point number from its raw underlying value.
    // Do not use this unless you know what you're doing.
    static const inline fixed from_raw_value(BaseType value)
    {
        return fixed(value, raw_construct_tag{});
    }

    //
    // Constants
    //
    static const fixed e() { return from_fixed_point<61>(6267931151224907085ll); }
    static const fixed pi() { return from_fixed_point<61>(7244019458077122842ll); }
    static const fixed half_pi() { return from_fixed_point<62>(7244019458077122842ll); }
    static const fixed two_pi() { return from_fixed_point<60>(7244019458077122842ll); }

    //
    // Arithmetic member operators
    //

    const inline fixed operator-() const
    {
        return fixed::from_raw_value(-m_value);
    }

    inline fixed& operator+=(const fixed& y)
    {
        m_value += y.m_value;
        return *this;
    }

	template <typename I, typename std::enable_if<std::is_floating_point<I>::value>::type* = nullptr>
	inline fixed& operator+=(I y)
	{
		m_value += fixed(y);
		return *this;
	}

    template <typename I, typename std::enable_if<std::is_integral<I>::value>::type* = nullptr>
    inline fixed& operator+=(I y)
    {
        m_value += y * FRACTION_MULT;
        return *this;
    }

    inline fixed& operator-=(const fixed& y)
    {
        m_value -= y.m_value;
        return *this;
    }

	template <typename I, typename std::enable_if<std::is_floating_point<I>::value>::type* = nullptr>
	inline fixed& operator-=(I y)
	{
		m_value -= fixed(y);
		return *this;
	}

    template <typename I, typename std::enable_if<std::is_integral<I>::value>::type* = nullptr>
    inline fixed& operator-=(I y)
    {
        m_value -= y * FRACTION_MULT;
        return *this;
    }

    inline fixed& operator*=(const fixed& y)
    {
        // Normal fixed-point multiplication is: x * y / 2**FractionBits.
        // To correctly round the last bit in the result, we need one more bit of information.
        // We do this by multiplying by two before dividing and adding the LSB to the real result.
        auto value = (static_cast<IntermediateType>(m_value) * y.m_value) / (FRACTION_MULT / 2);
        m_value = static_cast<BaseType>((value / 2) + (value % 2));
        return *this;
    }

	template <typename I, typename std::enable_if<std::is_floating_point<I>::value>::type* = nullptr>
	inline fixed& operator*=(I y)
	{
		m_value *= fixed(y);
		return *this;
	}

    template <typename I, typename std::enable_if<std::is_integral<I>::value>::type* = nullptr>
    inline fixed& operator*=(I y)
    {
        m_value *= y;
        return *this;
    }

    inline fixed& operator/=(const fixed& y)
    {
        assert(y.m_value != 0);
        // Normal fixed-point division is: x * 2**FractionBits / y.
        // To correctly round the last bit in the result, we need one more bit of information.
        // We do this by multiplying by two before dividing and adding the LSB to the real result.
        auto value = (static_cast<IntermediateType>(m_value) * FRACTION_MULT * 2) / y.m_value;
        m_value = static_cast<BaseType>((value / 2) + (value % 2));
        return *this;
    }

	template <typename I, typename std::enable_if<std::is_floating_point<I>::value>::type* = nullptr>
	inline fixed& operator/=(I y)
	{
		m_value /= fixed(y);
		return *this;
	}

    template <typename I, typename std::enable_if<std::is_integral<I>::value>::type* = nullptr>
    inline fixed& operator/=(I y)
    {
        m_value /= y;
        return *this;
    }

    std::string toString()
    {
        std::string str(std::to_string(m_value));
        str += "(" + std::to_string(static_cast<float>(m_value) / FRACTION_MULT) +")";
        return str;
    }

private:
    BaseType m_value;
};

//
// Convenience typedefs
//

using fixed_16_16 = fixed<std::int32_t, std::int64_t, 16>;
using fixed_24_8 = fixed<std::int32_t, std::int64_t, 8>;
using fixed_8_24 = fixed<std::int32_t, std::int64_t, 24>;

//
// Addition
//

template <typename B, typename I, unsigned int F>
const inline fixed<B, I, F> operator+(const fixed<B, I, F>& x, const fixed<B, I, F>& y)
{
    return fixed<B, I, F>(x) += y;
}

template <typename B, typename I, unsigned int F, typename T, typename std::enable_if<std::is_floating_point<T>::value>::type* = nullptr>
const inline fixed<B, I, F> operator+(const fixed<B, I, F>& x, T y)
{
	return fixed<B, I, F>(x) += fixed<B, I, F>(y);
}

template <typename B, typename I, unsigned int F, typename T, typename std::enable_if<std::is_integral<T>::value>::type* = nullptr>
const inline fixed<B, I, F> operator+(const fixed<B, I, F>& x, T y)
{
    return fixed<B, I, F>(x) += y;
}

template <typename B, typename I, unsigned int F, typename T, typename std::enable_if<std::is_floating_point<T>::value>::type* = nullptr>
const inline fixed<B, I, F> operator+(T x, const fixed<B, I, F>& y)
{
	return fixed<B, I, F>(y) += fixed<B, I, F>(x);
}

template <typename B, typename I, unsigned int F, typename T, typename std::enable_if<std::is_integral<T>::value>::type* = nullptr>
const inline fixed<B, I, F> operator+(T x, const fixed<B, I, F>& y)
{
    return fixed<B, I, F>(y) += x;
}

//
// Subtraction
//

template <typename B, typename I, unsigned int F>
const inline fixed<B, I, F> operator-(const fixed<B, I, F>& x, const fixed<B, I, F>& y)
{
    return fixed<B, I, F>(x) -= y;
}

template <typename B, typename I, unsigned int F, typename T, typename std::enable_if<std::is_floating_point<T>::value>::type* = nullptr>
const inline fixed<B, I, F> operator-(const fixed<B, I, F>& x, T y)
{
	return fixed<B, I, F>(x) -= fixed<B, I, F>(y);
}

template <typename B, typename I, unsigned int F, typename T, typename std::enable_if<std::is_integral<T>::value>::type* = nullptr>
const inline fixed<B, I, F> operator-(const fixed<B, I, F>& x, T y)
{
    return fixed<B, I, F>(x) -= y;
}

template <typename B, typename I, unsigned int F, typename T, typename std::enable_if<std::is_floating_point<T>::value>::type* = nullptr>
const inline fixed<B, I, F> operator-(T x, const fixed<B, I, F>& y)
{
	return fixed<B, I, F>(x) -= fixed<B, I, F>(y);
}

template <typename B, typename I, unsigned int F, typename T, typename std::enable_if<std::is_integral<T>::value>::type* = nullptr>
const inline fixed<B, I, F> operator-(T x, const fixed<B, I, F>& y)
{
    return fixed<B, I, F>(x) -= y;
}

//
// Multiplication
//

template <typename B, typename I, unsigned int F>
const inline fixed<B, I, F> operator*(const fixed<B, I, F>& x, const fixed<B, I, F>& y)
{
    return fixed<B, I, F>(x) *= y;
}

template <typename B, typename I, unsigned int F, typename T, typename std::enable_if<std::is_floating_point<T>::value>::type* = nullptr>
const inline fixed<B, I, F> operator*(const fixed<B, I, F>& x, T y)
{
	return fixed<B, I, F>(x) *= fixed<B, I, F>(y);
}

template <typename B, typename I, unsigned int F, typename T, typename std::enable_if<std::is_integral<T>::value>::type* = nullptr>
const inline fixed<B, I, F> operator*(const fixed<B, I, F>& x, T y)
{
    return fixed<B, I, F>(x) *= y;
}

template <typename B, typename I, unsigned int F, typename T, typename std::enable_if<std::is_floating_point<T>::value>::type* = nullptr>
const inline fixed<B, I, F> operator*(T x, const fixed<B, I, F>& y)
{
	return fixed<B, I, F>(y) *= fixed<B, I, F>(x);
}

template <typename B, typename I, unsigned int F, typename T, typename std::enable_if<std::is_integral<T>::value>::type* = nullptr>
const inline fixed<B, I, F> operator*(T x, const fixed<B, I, F>& y)
{
    return fixed<B, I, F>(y) *= x;
}

//
// Division
//

template <typename B, typename I, unsigned int F>
const inline fixed<B, I, F> operator/(const fixed<B, I, F>& x, const fixed<B, I, F>& y)
{
    return fixed<B, I, F>(x) /= y;
}

template <typename B, typename I, unsigned int F, typename T, typename std::enable_if<std::is_floating_point<T>::value>::type* = nullptr>
const inline fixed<B, I, F> operator/(const fixed<B, I, F>& x, T y)
{
	return fixed<B, I, F>(x) /= fixed<B, I, F>(y);
}

template <typename B, typename I, unsigned int F, typename T, typename std::enable_if<std::is_integral<T>::value>::type* = nullptr>
const inline fixed<B, I, F> operator/(const fixed<B, I, F>& x, T y)
{
    return fixed<B, I, F>(x) /= y;
}

template <typename B, typename I, unsigned int F, typename T, typename std::enable_if<std::is_floating_point<T>::value>::type* = nullptr>
const inline fixed<B, I, F> operator/(T x, const fixed<B, I, F>& y)
{
	return fixed<B, I, F>(x) /= fixed<B, I, F>(y);
}

template <typename B, typename I, unsigned int F, typename T, typename std::enable_if<std::is_integral<T>::value>::type* = nullptr>
const inline fixed<B, I, F> operator/(T x, const fixed<B, I, F>& y)
{
    return fixed<B, I, F>(x) /= y;
}

//
// Comparison operators
//

template <typename B, typename I, unsigned int F>
const inline bool operator==(const fixed<B, I, F>& x, const fixed<B, I, F>& y)
{
    return x.raw_value() == y.raw_value();
}

template <typename B, typename I, unsigned int F>
const inline bool operator!=(const fixed<B, I, F>& x, const fixed<B, I, F>& y)
{
    return x.raw_value() != y.raw_value();
}

template <typename B, typename I, unsigned int F>
const inline bool operator<(const fixed<B, I, F>& x, const fixed<B, I, F>& y)
{
    return x.raw_value() < y.raw_value();
}

template <typename B, typename I, unsigned int F>
const inline bool operator>(const fixed<B, I, F>& x, const fixed<B, I, F>& y)
{
    return x.raw_value() > y.raw_value();
}

template <typename B, typename I, unsigned int F>
const inline bool operator<=(const fixed<B, I, F>& x, const fixed<B, I, F>& y)
{
    return x.raw_value() <= y.raw_value();
}

template <typename B, typename I, unsigned int F>
const inline bool operator>=(const fixed<B, I, F>& x, const fixed<B, I, F>& y)
{
    return x.raw_value() >= y.raw_value();
}

namespace detail
{
// Number of base-10 digits required to fully represent a number of bits
static const int max_digits10(int bits)
{
    // 8.24 fixed-point equivalent of (int)ceil(bits * std::log10(2));
    using T = long long;
    return static_cast<int>((T{bits} * 5050445 + (T{1} << 24) - 1) >> 24);
}

// Number of base-10 digits that can be fully represented by a number of bits
static const int digits10(int bits)
{
    // 8.24 fixed-point equivalent of (int)(bits * std::log10(2));
    using T = long long;
    return static_cast<int>((T{bits} * 5050445) >> 24);
}

} // namespace detail
} // namespace fpm

// Specializations for customization points
namespace std
{

template <typename B, typename I, unsigned int F>
struct hash<fpm::fixed<B,I,F>>
{
    using argument_type = fpm::fixed<B, I, F>;
    using result_type = std::size_t;

    result_type operator()(argument_type arg) const{
        return m_hash(arg.raw_value());
    }

private:
    std::hash<B> m_hash;
};

template <typename B, typename I, unsigned int F>
struct numeric_limits<fpm::fixed<B,I,F>>
{
    static const bool is_specialized = true;
    static const bool is_signed = std::numeric_limits<B>::is_signed;
    static const bool is_integer = false;
    static const bool is_exact = true;
    static const bool has_infinity = false;
    static const bool has_quiet_NaN = false;
    static const bool has_signaling_NaN = false;
    static const bool has_denorm = std::denorm_absent;
    static const bool has_denorm_loss = false;
    static const std::float_round_style round_style = std::round_to_nearest;
    static const bool is_iec_559 = false;
    static const bool is_bounded = true;
    static const bool is_modulo = std::numeric_limits<B>::is_modulo;
    static const int digits = std::numeric_limits<B>::digits;

    // Any number with `digits10` significant base-10 digits (that fits in
    // the range of the type) is guaranteed to be convertible from text and
    // back without change. Worst case, this is 0.000...001, so we can only
    // guarantee this case. Nothing more.
    static const int digits10 = 1;

    // This is equal to max_digits10 for the integer and fractional part together.
    static const int max_digits10 =
        fpm::detail::max_digits10(std::numeric_limits<B>::digits - F) + fpm::detail::max_digits10(F);

    static const int radix = 2;
    static const int min_exponent = 1 - F;
    static const int min_exponent10 = -fpm::detail::digits10(F);
    static const int max_exponent = std::numeric_limits<B>::digits - F;
    static const int max_exponent10 = fpm::detail::digits10(std::numeric_limits<B>::digits - F);
    static const bool traps = true;
    static const bool tinyness_before = false;

    static const fpm::fixed<B,I,F> lowest() {
        return fpm::fixed<B,I,F>::from_raw_value(std::numeric_limits<B>::lowest());
    };

    static const fpm::fixed<B,I,F> min() {
        return lowest();
    }

    static const fpm::fixed<B,I,F> max() {
        return fpm::fixed<B,I,F>::from_raw_value(std::numeric_limits<B>::max());
    };

    static const fpm::fixed<B,I,F> epsilon() {
        return fpm::fixed<B,I,F>::from_raw_value(1);
    };

    static const fpm::fixed<B,I,F> round_error() {
        return fpm::fixed<B,I,F>(1) / 2;
    };

    static const fpm::fixed<B,I,F> denorm_min() {
        return min();
    }
};

template <typename B, typename I, unsigned int F>
const bool numeric_limits<fpm::fixed<B,I,F>>::is_specialized;
template <typename B, typename I, unsigned int F>
const bool numeric_limits<fpm::fixed<B,I,F>>::is_signed;
template <typename B, typename I, unsigned int F>
const bool numeric_limits<fpm::fixed<B,I,F>>::is_integer;
template <typename B, typename I, unsigned int F>
const bool numeric_limits<fpm::fixed<B,I,F>>::is_exact;
template <typename B, typename I, unsigned int F>
const bool numeric_limits<fpm::fixed<B,I,F>>::has_infinity;
template <typename B, typename I, unsigned int F>
const bool numeric_limits<fpm::fixed<B,I,F>>::has_quiet_NaN;
template <typename B, typename I, unsigned int F>
const bool numeric_limits<fpm::fixed<B,I,F>>::has_signaling_NaN;
template <typename B, typename I, unsigned int F>
const bool numeric_limits<fpm::fixed<B,I,F>>::has_denorm;
template <typename B, typename I, unsigned int F>
const bool numeric_limits<fpm::fixed<B,I,F>>::has_denorm_loss;
template <typename B, typename I, unsigned int F>
const std::float_round_style numeric_limits<fpm::fixed<B,I,F>>::round_style;
template <typename B, typename I, unsigned int F>
const bool numeric_limits<fpm::fixed<B,I,F>>::is_iec_559;
template <typename B, typename I, unsigned int F>
const bool numeric_limits<fpm::fixed<B,I,F>>::is_bounded;
template <typename B, typename I, unsigned int F>
const bool numeric_limits<fpm::fixed<B,I,F>>::is_modulo;
template <typename B, typename I, unsigned int F>
const int numeric_limits<fpm::fixed<B,I,F>>::digits;
template <typename B, typename I, unsigned int F>
const int numeric_limits<fpm::fixed<B,I,F>>::digits10;
template <typename B, typename I, unsigned int F>
const int numeric_limits<fpm::fixed<B,I,F>>::max_digits10;
template <typename B, typename I, unsigned int F>
const int numeric_limits<fpm::fixed<B,I,F>>::radix;
template <typename B, typename I, unsigned int F>
const int numeric_limits<fpm::fixed<B,I,F>>::min_exponent;
template <typename B, typename I, unsigned int F>
const int numeric_limits<fpm::fixed<B,I,F>>::min_exponent10;
template <typename B, typename I, unsigned int F>
const int numeric_limits<fpm::fixed<B,I,F>>::max_exponent;
template <typename B, typename I, unsigned int F>
const int numeric_limits<fpm::fixed<B,I,F>>::max_exponent10;
template <typename B, typename I, unsigned int F>
const bool numeric_limits<fpm::fixed<B,I,F>>::traps;
template <typename B, typename I, unsigned int F>
const bool numeric_limits<fpm::fixed<B,I,F>>::tinyness_before;

}

#endif
