//
// Created by steven on 2/4/16.
//

#ifndef DENOVOGEAR_BOOST_UTILS_H
#define DENOVOGEAR_BOOST_UTILS_H

namespace utils {


// Helper function that mimics boost::istream_range
template<class Elem, class Traits> inline
boost::iterator_range <std::istreambuf_iterator<Elem, Traits>> istreambuf_range(
        std::basic_istream <Elem, Traits> &in) {
    return boost::iterator_range < std::istreambuf_iterator < Elem, Traits >> (
            std::istreambuf_iterator<Elem, Traits>(in),
                    std::istreambuf_iterator<Elem, Traits>());
}


} // namespace utils

#endif //DENOVOGEAR_BOOST_UTILS_H
