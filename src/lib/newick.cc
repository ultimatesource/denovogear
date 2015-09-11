/*
 * Copyright (c) 2014-2015 Reed A. Cartwright
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
#define BOOST_SPIRIT_USE_PHOENIX_V3 1

#include <iostream>

#include <dng/newick.h>

#include <boost/fusion/adapted/std_pair.hpp>
#include <boost/fusion/include/std_pair.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/phoenix/core.hpp>
#include <boost/phoenix/operator.hpp>
#include <boost/phoenix/bind/bind_function.hpp>

#include <boost/algorithm/string/trim_all.hpp>

#define DNG_SM_PREFIX "SM-" // define also in pedigree.cc

using namespace dng::newick;

namespace qi = boost::spirit::qi;
namespace standard = boost::spirit::standard;
namespace phoenix = boost::phoenix;

struct node_t {
    std::string label;
    double length;
    std::size_t parent;
};

typedef std::vector<node_t> tree_t;
typedef std::size_t index_t;

struct make_tip_impl {
    typedef index_t result_type;

    index_t operator()(std::string label, double length, tree_t &g) const {
        index_t id = g.size();
        g.push_back({std::move(label), length, static_cast<std::size_t>(-1)});
        return id;
    }
};
const phoenix::function<make_tip_impl> make_tip;

struct make_inode_impl {
    typedef index_t result_type;

    index_t operator()(std::vector<index_t> v,
                       std::string label, double length, tree_t &g) const {
        using namespace boost;
        index_t id = g.size();
        g.push_back({std::move(label), length, static_cast<std::size_t>(-1)});
        for(auto a : v) {
            g[a].parent = id;
        }
        return id;
    }
};
const phoenix::function<make_inode_impl> make_inode;

template <typename Iterator>
struct newick_grammar :
qi::grammar<Iterator, void(tree_t &), standard::space_type> {
    // http://evolution.genetics.washington.edu/phylip/newick_doc.html
    newick_grammar() : newick_grammar::base_type(start) {
        using standard::char_; using qi::eps; using qi::attr;
        using qi::double_; using qi::lexeme; using qi::raw; using qi::omit;
        using qi::as_string;
        using qi::_1; using qi::_2; using qi::_3; using qi::_val; using qi::_r1;
        using standard::space;
        using phoenix::bind; using phoenix::val;
        using boost::vertex_distance; using boost::vertex_name;

        start    = -(node(_r1) || ';');
        node     = tip(_r1) | inode(_r1);

        length = (':' >> double_) | attr(1.0);

        tip      = (label >> length)[_val = make_tip(_1, _2, _r1)];
        inode    = (('(' >> (node(_r1) % ',') >> ')') >> ilabel >> length)[_val =
                       make_inode(_1,
                                  _2, _3, _r1)];

        ilabel = label | attr("");
        label    = unquoted | squoted | dquoted;
        unquoted = as_string[+(char_ - (char_(":,)(;'\"[]|") | space))];
        squoted   = as_string[lexeme['\'' >> *(char_ - '\'') >> '\'']];
        dquoted   = as_string[lexeme['"' >> *(char_ - '"') >> '"']];
    }

    qi::rule<Iterator, void(tree_t &), standard::space_type> start;
    qi::rule<Iterator, index_t(tree_t &), standard::space_type> node, tip, inode;
    qi::rule<Iterator, std::string(), standard::space_type> label, unquoted,
    squoted, dquoted, ilabel;
    qi::rule<Iterator, double(), standard::space_type> length;
};

int dng::newick::parse(const std::string &text, vertex_t root, Graph &graph) {
    using standard::space; using phoenix::ref;
    newick_grammar<std::string::const_iterator> newick_parser;
    std::string::const_iterator first = text.begin();
    qi::parse(first, text.end(), *space);
    if(first == text.end()) {
        return 0;
    }
    tree_t tree;
    bool r = qi::phrase_parse(first, text.end(), newick_parser(phoenix::ref(tree)),
                              space);
    if(first != text.end() || !r) {
        return -1;
    }
    std::size_t offset = num_vertices(graph) + tree.size() - 1;
    for(std::size_t i = tree.size(); i > 0; --i) {
        auto &&a = tree[i - 1];
        boost::trim_fill(a.label, "_");
        if(!a.label.empty()) {
            a.label = DNG_SM_PREFIX + a.label;
        }
        vertex_t v = add_vertex(a.label, graph);
        add_edge((a.parent == -1) ? root : offset - a.parent,
                 v, {dng::EdgeType::Mitotic, a.length}, graph);
    }

    return 1;
}
