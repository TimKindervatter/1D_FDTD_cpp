#pragma once

#include <vector>
#include <type_traits>
#include <cmath>
#include <algorithm>
#include <complex>
#include <execution>

// Array Creation

template<typename T>
std::vector<double> linspace(T start_in, T end_in, int num_in)
{
    std::vector<double> linspaced;

    double start = static_cast<double>(start_in);
    double end = static_cast<double>(end_in);
    double num = static_cast<double>(num_in);

    if (num == 0) { return linspaced; }
    if (num == 1)
    {
        linspaced.push_back(start);
        return linspaced;
    }

    double delta = (end - start) / (num - 1);

    for (int i = 0; i < num - 1; ++i)
    {
        linspaced.push_back(start + delta * i);
    }
    linspaced.push_back(end); // I want to ensure that start and end
                              // are exactly the same as the input
    return linspaced;
}

template<typename T>
std::vector<T> arange(T start, T end, T step = 1)
{
    size_t N = (int)floor(end - start);

    std::vector<T> vec;
    vec.reserve(N);
    for (T value = start; value < end; value += step)
        vec.push_back(value);

    return vec;
}

// Unary Math Broadcasts

template<typename T>
auto abs(const std::vector<T>& vec1)
{
    using vec_type = decltype(std::abs(std::declval<T>()));
    std::vector<vec_type> return_vector;
    return_vector.reserve(vec1.size());
    std::transform(std::begin(vec1),
        std::end(vec1),
        std::back_inserter(return_vector),
        [](auto elem) {return std::abs(elem); });

    return return_vector;
}
template<typename T>
auto squared_magnitude(const std::vector<T>& vec1)
{
    using vec_type = decltype(std::norm(std::declval<T>()));
    std::vector<vec_type> return_vector;
    return_vector.reserve(vec1.size());
    std::transform(std::begin(vec1),
        std::end(vec1),
        std::back_inserter(return_vector),
        [](auto elem) {return std::norm(elem); }); // std::real(elem) * std::real(elem) + std::imag(elem) * std::imag(elem)

    return return_vector;
}

template<typename T>
auto sqrt(const std::vector<T>& vec1)
{
    using vec_type = decltype(std::sqrt(std::declval<T>()));
    std::vector<vec_type> return_vector;
    return_vector.reserve(vec1.size());
    std::transform(std::begin(vec1),
        std::end(vec1),
        std::back_inserter(return_vector),
        [](auto elem) {return std::sqrt(elem); });

    return return_vector;
}

template<typename T>
auto exp(const std::vector<T>& vec1)
{
    using vec_type = decltype(std::exp(std::declval<T>()));
    std::vector<vec_type> return_vector;
    return_vector.reserve(vec1.size());
    std::transform(std::begin(vec1),
        std::end(vec1),
        std::back_inserter(return_vector),
        [](auto elem) {return std::exp(elem); });

    return return_vector;
}

// Max and Min

template<typename T>
auto max(const std::vector<T>& vec, T default_value = T{})
{
    auto max_element_iter = std::max_element(std::begin(vec), std::end(vec));

    T max_element{ default_value };
    if (max_element_iter != std::end(vec))
    {
        max_element = *max_element_iter;
    }

    return max_element;
}

template<typename T>
auto min(const std::vector<T>& vec, T default_value = T{})
{
    auto min_element_iter = std::min_element(std::begin(vec), std::end(vec));

    T min_element{ default_value };
    if (min_element_iter != std::end(vec))
    {
        min_element = *min_element_iter;
    }

    return min_element;
}

// Binary Math Broadcasts

template<typename T1, typename T2>
auto pow(const std::vector<T1>& vec1, T2 power)
{
    using vec_type = decltype(std::pow(std::declval<T1>(), std::declval<T2>()));
    std::vector<vec_type> return_vector;
    return_vector.reserve(vec1.size());
    std::transform(std::begin(vec1),
        std::end(vec1),
        std::back_inserter(return_vector),
        [power](auto elem) {return std::pow(elem, power); });

    return return_vector;
}


// Arithmetic Operations

// +

template<typename T1, typename T2>
auto operator+(const std::vector<T1>& vec1, T2 scalar)
{
    using vec_type = decltype(std::declval<T1>() + std::declval<T2>());
    std::vector<vec_type> return_vector;
    return_vector.reserve(vec1.size());
    std::transform(std::begin(vec1),
        std::end(vec1),
        std::back_inserter(return_vector),
        [scalar](auto elem) {return elem + scalar; });

    return return_vector;
}

template<typename T1, typename T2>
auto operator+(T1 scalar, const std::vector<T2>& vec1)
{
    return operator+(vec1, scalar);
}

template<typename T1, typename T2>
auto operator+(const std::vector<T1>& vec1, const std::vector<T2>& vec2)
{
    using vec_type = decltype(std::declval<T1>() + std::declval<T2>());
    std::vector<vec_type> return_vector;
    return_vector.reserve(vec1.size());
    std::transform(std::begin(vec1),
        std::end(vec1),
        std::begin(vec2),
        std::back_inserter(return_vector),
        [](auto elem1, auto elem2) {return elem1 + elem2; });

    return return_vector;
}

// Binary -

template<typename T1, typename T2>
auto operator-(const std::vector<T1>& vec1, T2 scalar)
{
    using vec_type = decltype(std::declval<T1>() - std::declval<T2>());
    std::vector<vec_type> return_vector;
    return_vector.reserve(vec1.size());
    std::transform(std::begin(vec1),
        std::end(vec1),
        std::back_inserter(return_vector),
        [scalar](auto elem) {return elem - scalar; });

    return return_vector;
}

template<typename T1, typename T2>
auto operator-(T1 scalar, const std::vector<T2>& vec1)
{
    using vec_type = decltype(std::declval<T1>() - std::declval<T2>());
    std::vector<vec_type> return_vector;
    return_vector.reserve(vec1.size());
    std::transform(std::begin(vec1),
        std::end(vec1),
        std::back_inserter(return_vector),
        [scalar](auto elem) {return scalar - elem; });

    return return_vector;
}

template<typename T1, typename T2>
auto operator-(const std::vector<T1>& vec1, const std::vector<T2>& vec2)
{
    using vec_type = decltype(std::declval<T1>() - std::declval<T2>());
    std::vector<vec_type> return_vector;
    return_vector.reserve(vec1.size());
    std::transform(std::begin(vec1),
        std::end(vec1),
        std::begin(vec2),
        std::back_inserter(return_vector),
        [](auto elem1, auto elem2) {return elem1 - elem2; });

    return return_vector;
}

// *

template<typename T1, typename T2>
auto operator*(const std::vector<T1>& vec1, T2 scalar)
{
    using vec_type = decltype(std::declval<T1>() * std::declval<T2>());
    std::vector<vec_type> return_vector;
    return_vector.reserve(vec1.size());
    std::transform(std::begin(vec1),
        std::end(vec1),
        std::back_inserter(return_vector),
        [scalar](auto elem) {return elem * scalar; });

    return return_vector;
}

template<> inline
auto operator*<std::complex<double>, double>(const std::vector<std::complex<double>>& vec1, double scalar)
{
    using vec_type = decltype(std::declval<std::complex<double>>() * std::declval<double>());
    std::vector<vec_type> return_vector;
    return_vector.reserve(vec1.size());
    std::transform(std::begin(vec1),
        std::end(vec1),
        std::back_inserter(return_vector),
        [scalar](auto elem) {return std::complex{ elem.real() * scalar, elem.imag() * scalar }; });

    return return_vector;
}

template<typename T1, typename T2>
auto operator*(T1 scalar, const std::vector<T2>& vec1)
{
    return operator*(vec1, scalar);
}

template<typename T1, typename T2>
auto operator*(const std::vector<T1>& vec1, const std::vector<T2>& vec2)
{
    using vec_type = decltype(std::declval<T1>() * std::declval<T2>());
    std::vector<vec_type> return_vector;
    return_vector.reserve(vec1.size());
    std::transform(std::begin(vec1),
        std::end(vec1),
        std::begin(vec2),
        std::back_inserter(return_vector),
        [](auto elem1, auto elem2) {return elem1 * elem2; });

    return return_vector;
}

// /

template<typename T1, typename T2>
auto operator/(const std::vector<T1>& vec1, T2 scalar)
{
    using vec_type = decltype(std::declval<T1>() / std::declval<T2>());
    std::vector<vec_type> return_vector;
    return_vector.reserve(vec1.size());
    std::transform(std::begin(vec1),
        std::end(vec1),
        std::back_inserter(return_vector),
        [scalar](auto elem) {return elem / scalar; });

    return return_vector;
}

template<typename T1, typename T2>
auto operator/(T1 scalar, const std::vector<T2>& vec1)
{
    using vec_type = decltype(std::declval<T1>() / std::declval<T2>());
    std::vector<vec_type> return_vector;
    return_vector.reserve(vec1.size());
    std::transform(std::begin(vec1),
        std::end(vec1),
        std::back_inserter(return_vector),
        [scalar](auto elem) {return scalar / elem; });

    return return_vector;
}

template<typename T1, typename T2>
auto operator/(const std::vector<T1>& vec1, const std::vector<T2>& vec2)
{
    using vec_type = decltype(std::declval<T1>() / std::declval<T2>());
    std::vector<vec_type> return_vector;
    return_vector.reserve(vec1.size());
    std::transform(std::begin(vec1),
        std::end(vec1),
        std::begin(vec2),
        std::back_inserter(return_vector),
        [](auto elem1, auto elem2) {return elem1 / elem2; });

    return return_vector;
}

// Unary -

template<typename T1>
auto operator-(const std::vector<T1>& vec1)
{
    return operator*(-1, vec1);
}