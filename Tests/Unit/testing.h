/*
 * Copyright (c) 2017 Reed A. Cartwright
 * Authors:  Reed A. Cartwright <reed@cartwrig.ht>
 *
 * This file is part of DeNovoGear.
 *
 * DeNovoGear is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once
#ifndef TESTS_UNIT_TESTING_H
#define TESTS_UNIT_TESTING_H

#include <cstddef>
#include <type_traits>
#include <boost/test/unit_test.hpp>
#include <boost/preprocessor/seq/for_each.hpp>
#include <boost/version.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/range/size.hpp>

#define GETTER1( C, V ) \
    static auto V(const C& x) -> const decltype(x.V##_)& { return x.V##_; } \
    static auto V(C& x) -> decltype(x.V##_)& { return x.V##_; } \
/**/

#define GETTER2( C, V, S ) \
    static auto S(C& x) -> decltype(x.V##_.S)& { return x.V##_.S; } \
/**/

#define GETTERF( C, F ) \
    template<typename ... Args> \
    static auto F(C &x, Args... args) -> decltype(x.F(args...)) { return x.F(args...); } \
/**/

#define GETTER1_ELEM( R, C, V ) \
    GETTER1( C, V ) \
/**/

#define GETTERS_FOR_MEMBER_VARIABLES( C, S ) \
    BOOST_PP_SEQ_FOR_EACH(GETTER1_ELEM, C, S) \
/**/


#if defined(BOOST_VERSION) && BOOST_VERSION < 105900
// Version for Boost.Test 2

#define CHECK_EQUAL_RANGES( L, R ) \
    BOOST_CHECK_EQUAL_COLLECTIONS( ::std::begin( L ), ::std::end( L ), \
                                   ::std::begin( R ), ::std::end( R )) \
/**/

#define CHECK_CLOSE_RANGES( L, R, T )    do { \
    BOOST_REQUIRE_EQUAL( ::boost::size( L ), ::boost::size( R ) ); \
    auto CHECK_CLOSE_RANGES_A = ::std::begin( L ); \
    auto CHECK_CLOSE_RANGES_B = ::std::begin( R ); \
    ::std::size_t CHECK_CLOSE_RANGES_N = 0; \
    while(CHECK_CLOSE_RANGES_A < ::std::end( L )) { \
        if(*CHECK_CLOSE_RANGES_A == 0) { \
            BOOST_CHECK_SMALL(*CHECK_CLOSE_RANGES_B, T); \
        } else if(*CHECK_CLOSE_RANGES_B == 0 ) { \
            BOOST_CHECK_SMALL(*CHECK_CLOSE_RANGES_A, T); \
        } else { \
            BOOST_CHECK_CLOSE_FRACTION(*CHECK_CLOSE_RANGES_A, *CHECK_CLOSE_RANGES_B, T); \
        } \
        ++CHECK_CLOSE_RANGES_A; \
        ++CHECK_CLOSE_RANGES_B; \
    } \
    } while(false) \
/**/

#define CHECK_NE_RANGES( L, R )    do { \
    BOOST_REQUIRE_EQUAL( ::boost::size( L ), ::boost::size( R ) ); \
    auto CHECK_CLOSE_RANGES_A = ::std::begin( L ); \
    auto CHECK_CLOSE_RANGES_B = ::std::begin( R ); \
    ::std::size_t CHECK_CLOSE_RANGES_N = 0; \
    while(CHECK_CLOSE_RANGES_A < ::std::end( L )) { \
        BOOST_CHECK_NE(*CHECK_CLOSE_RANGES_A, *CHECK_CLOSE_RANGES_B); \
        ++CHECK_CLOSE_RANGES_A; \
        ++CHECK_CLOSE_RANGES_B; \
    } \
    } while(false) \
/**/

#ifndef BOOST_TEST_CONTEXT
#define BOOST_TEST_CONTEXT(D) \
    if(true) \
/**/
#endif

#ifndef BOOST_TEST_INFO
#define BOOST_TEST_INFO(D) \
/**/
#endif

#else
// Version for modern Boost.Test
#define CHECK_EQUAL_RANGES( L, R ) \
    BOOST_TEST(L == R, ::boost::test_tools::per_element() ) \
/**/

#define CHECK_NE_RANGES( L, R ) \
    BOOST_TEST(L != R, ::boost::test_tools::per_element() ) \
/**/

#define CHECK_CLOSE_RANGES( L, R, T )    do { \
    ::boost::test_tools::local_fpc_tolerance<double> t_o_l( T ); \
    BOOST_TEST( L == R, ::boost::test_tools::per_element() ); \
    } while(false) \
/**/
#endif

namespace dng {
namespace detail {

// RAII class that creates and opens a temporary file for reading and writing
// and deletes it when done.
struct AutoTempFile {
    boost::filesystem::path path;
    boost::filesystem::fstream file;

    AutoTempFile() {
        using namespace boost::filesystem;
        path = temp_directory_path();
        path /= unique_path();
        file.open(path, std::ios::out); //create file
        file.close();
        file.open(path); // open rw
    } 
    ~AutoTempFile() {
        using namespace boost::filesystem;
        file.close();
        remove(path);
    }
};

// function to create a 

namespace {
template<typename T>
struct test_range {
    using iterator = T;
    using const_iterator = T;
    using value_type = typename std::iterator_traits<T>::value_type;
    using reference = typename std::iterator_traits<T>::reference;
    using difference_type = typename std::iterator_traits<T>::difference_type;
    using pointer = typename std::iterator_traits<T>::pointer;
    using iterator_category = typename std::iterator_traits<T>::iterator_category;
    
    test_range(T b, T e) : begin_{b}, end_{e} { }

    const_iterator begin() const { return begin_; }
    const_iterator end() const { return end_; }
    difference_type size() const { return std::distance(begin(),end()); }

    const_iterator begin_;
    const_iterator end_;
};


template<typename M>
auto make_test_range(const M& m) -> test_range<decltype(m.data())>
{
    return {m.data(), m.data()+m.size()};
}

template<typename M>
auto make_test_range(M& m) -> test_range<decltype(m.data())>
{
    return {m.data(), m.data()+m.size()};
}

template<typename M>
auto make_test_range(const M& m) -> test_range<decltype(m.begin())>
{
    return {m.begin(), m.end()};
}

template<typename M>
auto make_test_range(M& m) -> test_range<decltype(m.begin())>
{
    return {m.begin(), m.end()};
}


template<typename T>
auto make_test_range(const T* b, const T* e) -> test_range<const T*>
{
    return {b, e};
}

template<typename T, size_t N>
auto make_test_range(const T (&a)[N]) -> test_range<const T*>
{
    return {&a[0], &a[N]};
}

} // anonymous namespace

} // namespace dng::detail
} // namespace dng

namespace dng {
namespace detail {
namespace io {
// allow enum's to be printed
template<typename CharType, typename CharTrait, typename E>
inline
typename std::enable_if<std::is_enum<E>::value, std::basic_ostream<CharType, CharTrait> >::type&
operator<<(std::basic_ostream<CharType, CharTrait>& o, E m) {
    // This will probably not do what we want if enum's int type is char.
    o << static_cast<typename std::underlying_type<E>::type>(m);
    return o;
}
}}} // namespace dng::detail::io

// allow pair of strings to be printed
namespace std {
template<typename CharType, typename CharTrait>
inline
typename std::basic_ostream<CharType, CharTrait>&
operator<<(std::basic_ostream<CharType, CharTrait>& o, const std::pair<std::string,std::string>& m) {
    o << "{'" << m.first << "', '" << m.second << "'}";
    return o;
}

}

#endif
