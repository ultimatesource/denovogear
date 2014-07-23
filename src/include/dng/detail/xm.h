/******************************************************************************
 * Copyright (c) 2007-2009 Reed A. Cartwright, PhD <reed@scit.us>             *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER *
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING    *
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER        *
 * DEALINGS IN THE SOFTWARE.                                                  *
 ******************************************************************************/

#include <boost/preprocessor.hpp>

#ifndef XMACROS_HELPERS
#define XMACROS_HELPERS

/******************************************************************************
 *    X-Helpers List                                                          *
 ******************************************************************************/

// The JS macro cats a seq 'seq' with separator 'sep'

#define JS_OP(s, data, elem) BOOST_PP_CAT(data, elem)

#define JS(sep, seq) BOOST_PP_IF( BOOST_PP_EQUAL(BOOST_PP_SEQ_SIZE(seq),1), \
	JS_1, JS_2)(sep,seq)

#define JS_1(sep, seq) BOOST_PP_SEQ_HEAD(seq)

#define JS_2(sep, seq) BOOST_PP_SEQ_CAT(( BOOST_PP_SEQ_HEAD(seq) ) \
	BOOST_PP_SEQ_TRANSFORM(JS_OP, sep, BOOST_PP_SEQ_TAIL(seq)) \
)

// The SS macro is similiar to _JS except that it stringizes everything

#define SS_OP(r, data, elem) data BOOST_PP_STRINGIZE(elem)

#define SS(sep, seq) BOOST_PP_STRINGIZE(BOOST_PP_SEQ_HEAD(seq)) \
	BOOST_PP_SEQ_FOR_EACH(SS_OP, sep, BOOST_PP_SEQ_TAIL(seq))

// Output result if seq is defined

#define IFD(seq,res) BOOST_PP_EXPR_IF(BOOST_PP_GREATER(BOOST_PP_SEQ_SIZE((_)seq),1), res)

#define XV(lname) JS(_, lname)
#define XS(lname) SS("-", lname)
#define XP(lname) SS(".", lname)
#define DL(a,b) a,b

#else

/******************************************************************************
 *    Cleanup                                                                 *
 ******************************************************************************/

#undef XMACROS_HELPERS
#undef JS_OP
#undef JS
#undef JS_1
#undef JS_2
#undef SS_OP
#undef SS
#undef IFD
#undef XV
#undef XS
#undef XP
#undef DL

#endif

