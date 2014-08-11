/*
 * Copyright (c) 2014 Reed A. Cartwright
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

#include <dng/newick.h>

#include <boost/spirit/include/qi.hpp>

#include <boost/phoenix/core.hpp>
#include <boost/phoenix/operator.hpp>
#include <boost/phoenix/bind/bind_function.hpp>

using namespace dng::newick;

namespace qi = boost::spirit::qi;
namespace standard = boost::spirit::standard;
namespace phoenix = boost::phoenix;

typedef dng::Graph tree_t;
typedef boost::graph_traits<tree_t>::vertex_descriptor node_t;
typedef boost::graph_traits<tree_t>::edge_descriptor branch_t;

struct make_tip_impl {
	typedef node_t result_type;

	node_t operator()(std::string label, tree_t &g) const {
		using namespace boost;
		auto id = add_vertex(g);
		put(vertex_name,g,id,label);
		put(vertex_distance,g,id,0.0);
		return id;
	}
};
const phoenix::function<make_tip_impl> make_tip;

struct make_inode_impl {
	typedef node_t result_type;

	node_t operator()(std::vector<node_t> v, tree_t &g) const {
		using namespace boost;
		auto id = add_vertex(g);
		for(auto a : v)
			add_edge(id,a,g);
		put(vertex_name,g,id,"");
		put(vertex_distance,g,id,0.0);
		return id;
	}
};
const phoenix::function<make_inode_impl> make_inode;

struct set_prop_impl {
	typedef void result_type;

	template<typename Prop, typename Val>
	void operator()(Prop p, tree_t &g, node_t id, Val v) const {
		using namespace boost;
		put(p,g,id,v);
	}
};
const phoenix::function<set_prop_impl> set_prop;

template <typename Iterator>
struct newick_grammar : qi::grammar<Iterator, tree_t(), standard::space_type> {
	// http://evolution.genetics.washington.edu/phylip/newick_doc.html
	newick_grammar() : newick_grammar::base_type(start) {
		using standard::char_; using qi::eps;
		using qi::float_; using qi::lexeme; using qi::raw; using qi::omit;
		using qi::_1; using qi::_2; using qi::_val; using qi::_r1;
		using standard::space;
		using phoenix::bind; using phoenix::val;
		using boost::vertex_distance; using boost::vertex_name;
		
		start    = omit[node(_val)] >> ';';
		node     = tip(_r1) | inode(_r1);
		
		tip      = label[_val=make_tip(_1,_r1)]
		           >> -(':' >> float_[set_prop(vertex_distance,_r1,_val,_1)]);
		inode    = ('(' >> (node(_r1) % ',') >> ')')[_val=make_inode(_1,_r1)]
					    >> -(label[set_prop(vertex_name,_r1,_val,_1)]
					    || (':' >> float_[set_prop(vertex_distance,_r1,_val,_1)]));
		label    = unquoted | squoted | dquoted;
		unquoted = lexeme[+(char_ - (char_(":,)(;'[]|")|space))];
		squoted   = raw[lexeme['\'' >> *(char_ - '\'') >> '\'']];
		dquoted   = raw[lexeme['"' >> *(char_ - '"') >> '"']];
	}
	
	qi::rule<Iterator, tree_t(), standard::space_type> start;
	qi::rule<Iterator, node_t(tree_t&), standard::space_type> node,tip,inode;
	qi::rule<Iterator, std::string(), standard::space_type> label,unquoted,squoted,dquoted;
};


bool dng::newick::parse(const std::string& text, Graph& graph) {
	using standard::space;
	newick_grammar<std::string::const_iterator> newick_parser;
	std::string::const_iterator first = text.begin();
	bool r = qi::phrase_parse(first, text.end(), newick_parser, space, graph);
	if( first != text.end() || !r ) {
		return false;
	}
	return true;
}
