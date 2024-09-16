#ifndef MY_TYPE_INFO_H_
#define MY_TYPE_INFO_H_

#include <tuple> // std::tuple
#include <array> // std::array
#include <vector> // std::vector
#include <string> // std::string
#include <concepts> // std::same_as
#include <type_traits> // std::{true,false}_type, std::is_array

template <class T>
struct is_tuple : std::false_type {};

template <class ...Ts>
struct is_tuple<std::tuple<Ts...>> : std::true_type {};

template <class T>
inline constexpr bool is_tuple_v = is_tuple<T>::value;

template <class T>
struct is_std_array : std::false_type {};

template <class T, size_t N>
struct is_std_array<std::array<T, N>> : std::true_type {};

template <class T>
inline constexpr bool is_array_type = is_std_array<T>::value || std::is_array<T>::value;

template <class T>
concept array_type = is_array_type<T>;

template <array_type T>
struct array_info;

template <class T, size_t N>
struct array_info<T[N]>
{
    using value_type = T;
    static constexpr size_t size = N;
};

template <class T, size_t N>
struct array_info<std::array<T, N>>
{
    using value_type = T;
    static constexpr size_t size = N;
};

template <class T>
inline constexpr size_t array_size = array_info<T>::size;

template <class T>
concept string_type = std::same_as<T, std::string> || std::same_as<T, std::vector<char>>;

template <class T>
concept real_type = std::same_as<T, float> || std::same_as<T, double>;

template <class T>
concept index_type = std::is_integral_v<T>;

#endif
