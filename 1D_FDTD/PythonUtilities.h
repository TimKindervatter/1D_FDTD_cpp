//#pragma once
//
//#include <vector>
//#include <type_traits>
//#include <cmath>
//#include <algorithm>
//#include <complex>
//#include <execution>
//
//// Array Creation
//
//template<typename T1, typename T2>
//std::vector<floating_point_t> linspace(T1 start_in, T2 end_in, int num_in)
//{
//    std::vector<floating_point_t> linspaced;
//
//    floating_point_t start = static_cast<floating_point_t>(start_in);
//    floating_point_t end = static_cast<floating_point_t>(end_in);
//    floating_point_t num = static_cast<floating_point_t>(num_in);
//
//    if (num == 0) { return linspaced; }
//    if (num == 1)
//    {
//        linspaced.push_back(start);
//        return linspaced;
//    }
//
//    floating_point_t delta = (end - start) / (num - 1);
//
//    for (int i = 0; i < num - 1; ++i)
//    {
//        linspaced.push_back(start + delta * i);
//    }
//    linspaced.push_back(end); // I want to ensure that start and end
//                              // are exactly the same as the input
//    return linspaced;
//}
//
//template<typename T>
//std::vector<T> arange(T start, T end, T step = 1)
//{
//    size_t N = (int)floor(end - start);
//
//    std::vector<T> vec;
//    vec.reserve(N);
//    for (T value = start; value < end; value += step)
//        vec.push_back(value);
//
//    return vec;
//}
//
//// Unary Math Broadcasts
//
//template<typename T>
//auto abs(const std::vector<T>& vec1)
//{
//    using vec_type = decltype(std::abs(std::declval<T>()));
//    std::vector<vec_type> return_vector;
//    return_vector.reserve(vec1.size());
//    std::transform(std::begin(vec1),
//        std::end(vec1),
//        std::back_inserter(return_vector),
//        [](auto elem) {return std::abs(elem); });
//
//    return return_vector;
//}
//template<typename T>
//auto squared_magnitude(const std::vector<T>& vec1)
//{
//    using vec_type = decltype(std::norm(std::declval<T>()));
//    std::vector<vec_type> return_vector;
//    return_vector.reserve(vec1.size());
//    std::transform(std::begin(vec1),
//        std::end(vec1),
//        std::back_inserter(return_vector),
//        [](auto elem) {return std::norm(elem); }); // std::real(elem) * std::real(elem) + std::imag(elem) * std::imag(elem)
//
//    return return_vector;
//}
//
//template<typename T>
//auto sqrt(const std::vector<T>& vec1)
//{
//    using vec_type = decltype(std::sqrt(std::declval<T>()));
//    std::vector<vec_type> return_vector;
//    return_vector.reserve(vec1.size());
//    std::transform(std::begin(vec1),
//        std::end(vec1),
//        std::back_inserter(return_vector),
//        [](auto elem) {return std::sqrt(elem); });
//
//    return return_vector;
//}
//
//template<typename T>
//auto exp(const std::vector<T>& vec1)
//{
//    using vec_type = decltype(std::exp(std::declval<T>()));
//    std::vector<vec_type> return_vector;
//    return_vector.reserve(vec1.size());
//    std::transform(std::begin(vec1),
//        std::end(vec1),
//        std::back_inserter(return_vector),
//        [](auto elem) {return std::exp(elem); });
//
//    return return_vector;
//}
//
//// Max and Min
//
//template<typename T1, typename T2>
//auto max(const std::vector<T1>& vec, T2 default_value = T{})
//{
//    auto max_element_iter = std::max_element(std::begin(vec), std::end(vec));
//
//    using promoted_type = decltype(std::declval<T1>() + std::declval<T2>());
//    promoted_type max_element{ default_value };
//    if (max_element_iter != std::end(vec))
//    {
//        max_element = *max_element_iter;
//    }
//
//    return max_element;
//}
//
//template<typename T1, typename T2>
//auto min(const std::vector<T1>& vec, T2 default_value = T{})
//{
//    auto min_element_iter = std::min_element(std::begin(vec), std::end(vec));
//
//    using promoted_type = decltype(std::declval<T1>() + std::declval<T2>());
//    promoted_type min_element{ default_value };
//    if (min_element_iter != std::end(vec))
//    {
//        min_element = *min_element_iter;
//    }
//
//    return min_element;
//}
//
//// Binary Math Broadcasts
//
//template<typename T1, typename T2>
//auto pow(const std::vector<T1>& vec1, T2 power)
//{
//    using vec_type = decltype(std::pow(std::declval<T1>(), std::declval<T2>()));
//    std::vector<vec_type> return_vector;
//    return_vector.reserve(vec1.size());
//    std::transform(std::begin(vec1),
//        std::end(vec1),
//        std::back_inserter(return_vector),
//        [power](auto elem) {return std::pow(elem, power); });
//
//    return return_vector;
//}
//
//
//// Arithmetic Operations
//
//// +
//
//template<typename T1, typename T2>
//auto operator+(const std::vector<T1>& vec1, T2 scalar)
//{
//    using vec_type = decltype(std::declval<T1>() + std::declval<T2>());
//    std::vector<vec_type> return_vector;
//    return_vector.reserve(vec1.size());
//    std::transform(std::begin(vec1),
//        std::end(vec1),
//        std::back_inserter(return_vector),
//        [scalar](auto elem) {return elem + scalar; });
//
//    return return_vector;
//}
//
//template<typename T1, typename T2>
//auto operator+(T1 scalar, const std::vector<T2>& vec1)
//{
//    return operator+(vec1, scalar);
//}
//
//template<typename T1, typename T2>
//auto operator+(const std::vector<T1>& vec1, const std::vector<T2>& vec2)
//{
//    using vec_type = decltype(std::declval<T1>() + std::declval<T2>());
//    std::vector<vec_type> return_vector;
//    return_vector.reserve(vec1.size());
//    std::transform(std::begin(vec1),
//        std::end(vec1),
//        std::begin(vec2),
//        std::back_inserter(return_vector),
//        [](auto elem1, auto elem2) {return elem1 + elem2; });
//
//    return return_vector;
//}
//
//// Binary -
//
//template<typename T1, typename T2>
//auto operator-(const std::vector<T1>& vec1, T2 scalar)
//{
//    using vec_type = decltype(std::declval<T1>() - std::declval<T2>());
//    std::vector<vec_type> return_vector;
//    return_vector.reserve(vec1.size());
//    std::transform(std::begin(vec1),
//        std::end(vec1),
//        std::back_inserter(return_vector),
//        [scalar](auto elem) {return elem - scalar; });
//
//    return return_vector;
//}
//
//template<typename T1, typename T2>
//auto operator-(T1 scalar, const std::vector<T2>& vec1)
//{
//    using vec_type = decltype(std::declval<T1>() - std::declval<T2>());
//    std::vector<vec_type> return_vector;
//    return_vector.reserve(vec1.size());
//    std::transform(std::begin(vec1),
//        std::end(vec1),
//        std::back_inserter(return_vector),
//        [scalar](auto elem) {return scalar - elem; });
//
//    return return_vector;
//}
//
//template<typename T1, typename T2>
//auto operator-(const std::vector<T1>& vec1, const std::vector<T2>& vec2)
//{
//    using vec_type = decltype(std::declval<T1>() - std::declval<T2>());
//    std::vector<vec_type> return_vector;
//    return_vector.reserve(vec1.size());
//    std::transform(std::begin(vec1),
//        std::end(vec1),
//        std::begin(vec2),
//        std::back_inserter(return_vector),
//        [](auto elem1, auto elem2) {return elem1 - elem2; });
//
//    return return_vector;
//}
//
//// *
//
//template<typename T1, typename T2>
//auto operator*(const std::vector<T1>& vec1, T2 scalar)
//{
//    using vec_type = decltype(std::declval<T1>() * std::declval<T2>());
//    std::vector<vec_type> return_vector;
//    return_vector.reserve(vec1.size());
//    std::transform(std::begin(vec1),
//        std::end(vec1),
//        std::back_inserter(return_vector),
//        [scalar](auto elem) {return elem * scalar; });
//
//    return return_vector;
//}
//
//template<> inline
//auto operator*<std::complex<floating_point_t>, floating_point_t>(const std::vector<std::complex<floating_point_t>>& vec1, floating_point_t scalar)
//{
//    using vec_type = decltype(std::declval<std::complex<floating_point_t>>() * std::declval<floating_point_t>());
//    std::vector<vec_type> return_vector;
//    return_vector.reserve(vec1.size());
//    std::transform(std::begin(vec1),
//        std::end(vec1),
//        std::back_inserter(return_vector),
//        [scalar](auto elem) {return std::complex{ elem.real() * scalar, elem.imag() * scalar }; });
//
//    return return_vector;
//}
//
//template<typename T1, typename T2>
//auto operator*(T1 scalar, const std::vector<T2>& vec1)
//{
//    return operator*(vec1, scalar);
//}
//
//template<typename T1, typename T2>
//auto operator*(const std::vector<T1>& vec1, const std::vector<T2>& vec2)
//{
//    using vec_type = decltype(std::declval<T1>() * std::declval<T2>());
//    std::vector<vec_type> return_vector;
//    return_vector.reserve(vec1.size());
//    std::transform(std::begin(vec1),
//        std::end(vec1),
//        std::begin(vec2),
//        std::back_inserter(return_vector),
//        [](auto elem1, auto elem2) {return elem1 * elem2; });
//
//    return return_vector;
//}
//
//inline std::complex<floating_point_t> operator*(std::complex<floating_point_t> z, std::complex<floating_point_t> w)
//{
//    auto x = z.real();
//    auto y = z.imag();
//    auto u = w.real();
//    auto v = w.imag();
//
//    return std::complex{ x*u - y*v, x*v + y*u };
//}
//
//// /
//
//template<typename T1, typename T2>
//auto operator/(const std::vector<T1>& vec1, T2 scalar)
//{
//    using vec_type = decltype(std::declval<T1>() / std::declval<T2>());
//    std::vector<vec_type> return_vector;
//    return_vector.reserve(vec1.size());
//    std::transform(std::begin(vec1),
//        std::end(vec1),
//        std::back_inserter(return_vector),
//        [scalar](auto elem) {return elem / scalar; });
//
//    return return_vector;
//}
//
//template<typename T1, typename T2>
//auto operator/(T1 scalar, const std::vector<T2>& vec1)
//{
//    using vec_type = decltype(std::declval<T1>() / std::declval<T2>());
//    std::vector<vec_type> return_vector;
//    return_vector.reserve(vec1.size());
//    std::transform(std::begin(vec1),
//        std::end(vec1),
//        std::back_inserter(return_vector),
//        [scalar](auto elem) {return scalar / elem; });
//
//    return return_vector;
//}
//
//template<typename T1, typename T2>
//auto operator/(const std::vector<T1>& vec1, const std::vector<T2>& vec2)
//{
//    using vec_type = decltype(std::declval<T1>() / std::declval<T2>());
//    std::vector<vec_type> return_vector;
//    return_vector.reserve(vec1.size());
//    std::transform(std::begin(vec1),
//        std::end(vec1),
//        std::begin(vec2),
//        std::back_inserter(return_vector),
//        [](auto elem1, auto elem2) {return elem1 / elem2; });
//
//    return return_vector;
//}
//
//
//inline std::complex<floating_point_t> operator/(std::complex<floating_point_t> z, std::complex<floating_point_t> w)
//{
//    auto x = z.real();
//    auto y = z.imag();
//    auto u = w.real();
//    auto v = w.imag();
//
//    auto denominator = u * u + v * v;
//    auto real_numerator = x * u + y * v;
//    auto imag_numerator = -x * v + y * u;
//
//    auto real = real_numerator / denominator;
//    auto imag = imag_numerator / denominator;
//
//    return std::complex{ real, imag };
//}
//
//// Unary -
//
//template<typename T1>
//auto operator-(const std::vector<T1>& vec1)
//{
//    return operator*(-1, vec1);
//}