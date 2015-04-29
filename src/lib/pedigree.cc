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

#include <dng/pedigree.h>
#include <dng/graph.h>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/biconnected_components.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/range/algorithm/find.hpp>

// Calculate an Felsenstein (1981) Substitution Matrix based on rate and nucleotide frequency
// If bNoMut is true, mask out the off diagonal.
Eigen::Matrix4d f81_matrix(double mu, std::array<double,4> nuc_freq, bool bNoMut=false) {
	double beta = 1.0;
	for(auto d : nuc_freq)
		beta -= d*d;
	beta = 1.0/beta;
	beta = exp(-beta*mu);

	Eigen::Matrix4d m;
	for(int i : {0,1,2,3}) {
		for(int j : {0,1,2,3}) {
			m(i,j) = nuc_freq[i]*(1.0-beta);
		}
		m(i,i) += beta;
	}
	return (bNoMut) ? (m.diagonal()).asDiagonal() : m;
}

dng::TransitionMatrix mitosis_matrix(double mu, std::array<double,4> nuc_freq, bool bNoMut=false) {
	using dng::nucleotides;

	// construct haploid mutation matrix
	Eigen::Matrix4d m = f81_matrix(mu,nuc_freq,bNoMut);

	// Matrix to hold the return value
	dng::TransitionMatrix ret{10,10};

	// Populate the matrix by combining two haploid mutation matrices
	for(int p = 0; p < 10; ++p) {
		for(int k = 0; k < 10; ++k) {
			ret(p,k) = 0.5*( m(nucleotides[p][0],nucleotides[k][0])*m(nucleotides[p][1],nucleotides[k][1])
				            +m(nucleotides[p][0],nucleotides[k][1])*m(nucleotides[p][1],nucleotides[k][0]));
			if(nucleotides[p][0] == nucleotides[p][1])
				continue;
			ret(p,k) += 0.5*( m(nucleotides[p][1],nucleotides[k][0])*m(nucleotides[p][0],nucleotides[k][1])
				             +m(nucleotides[p][1],nucleotides[k][1])*m(nucleotides[p][0],nucleotides[k][0]));
		}
	}
	return ret;
}

// TODO: separate mother/father sub matrices
dng::TransitionMatrix meiosis_matrix(double mu, std::array<double,4> nuc_freq, bool bNoMut=false) {
	// Construct Mutation Process
	using dng::nucleotides;

	Eigen::Matrix4d m = f81_matrix(mu,nuc_freq,bNoMut);

	dng::TransitionMatrix ret{100,10};
	for(int i = 0; i < 10; ++i) {
		for(int j = 0; j < 10; ++j) {
			int p = 10*i+j;
			for(int k = 0; k < 10; ++k) {
				ret(p,k) =
					  ( m(nucleotides[i][0],nucleotides[k][0])
					  + m(nucleotides[i][1],nucleotides[k][0]))
				    * ( m(nucleotides[j][0],nucleotides[k][1])
				      + m(nucleotides[j][1],nucleotides[k][1]))/4.0 ;
				if(nucleotides[k][0] == nucleotides[k][1])
					continue;
				ret(p,k) +=
					  ( m(nucleotides[i][0],nucleotides[k][1])
					  + m(nucleotides[i][1],nucleotides[k][1]))
				    * ( m(nucleotides[j][0],nucleotides[k][0])
				      + m(nucleotides[j][1],nucleotides[k][0]))/4.0 ;
			}
		}
	}
	return ret;
}

dng::GenotypeArray population_prior(double theta, std::array<double,4> nuc_freq, std::array<double,4> prior) {
	dng::GenotypeArray ret{10};
	double alpha[4] = {
		theta*nuc_freq[0]+prior[0], theta*nuc_freq[1]+prior[1],
		theta*nuc_freq[2]+prior[2], theta*nuc_freq[3]+prior[3]
	};
	double alpha_sum = alpha[0]+alpha[1]+alpha[2]+alpha[3];
	ret <<
	        alpha[0]*(1.0+alpha[0])/alpha_sum/(1.0+alpha_sum), // AA
	    2.0*alpha[0]*(    alpha[1])/alpha_sum/(1.0+alpha_sum), // AC
	    2.0*alpha[0]*(    alpha[2])/alpha_sum/(1.0+alpha_sum), // AG
	    2.0*alpha[0]*(    alpha[3])/alpha_sum/(1.0+alpha_sum), // AT
	        alpha[1]*(1.0+alpha[1])/alpha_sum/(1.0+alpha_sum), // CC
	    2.0*alpha[1]*(    alpha[2])/alpha_sum/(1.0+alpha_sum), // CG
	    2.0*alpha[1]*(    alpha[3])/alpha_sum/(1.0+alpha_sum), // CT
	        alpha[2]*(1.0+alpha[2])/alpha_sum/(1.0+alpha_sum), // GG
	    2.0*alpha[2]*(    alpha[3])/alpha_sum/(1.0+alpha_sum), // GT
	        alpha[3]*(1.0+alpha[3])/alpha_sum/(1.0+alpha_sum); // GG
    return ret;
}


bool dng::Pedigree::Initialize(params_t p) {
	using namespace Eigen;
	using namespace std;
	// TODO: Validate params_t p or return false???

	// Construct Genotype Prior
	double alpha[4] = {
		p.theta*p.nuc_freq[0], p.theta*p.nuc_freq[1],
		p.theta*p.nuc_freq[2], p.theta*p.nuc_freq[3]
	};
	// Use a parent-independent mutation model, which produces a
	// beta-binomial
	genotype_prior_[0] = population_prior(p.theta, p.nuc_freq, {p.ref_weight,0,0,0});
	genotype_prior_[1] = population_prior(p.theta, p.nuc_freq, {0,p.ref_weight,0,0});
	genotype_prior_[2] = population_prior(p.theta, p.nuc_freq, {0,0,p.ref_weight,0});
	genotype_prior_[3] = population_prior(p.theta, p.nuc_freq, {0,0,0,p.ref_weight});
	genotype_prior_[4] = population_prior(p.theta, p.nuc_freq, {0,0,0,0});
	
	meiosis_ = meiosis_matrix(p.mu,p.nuc_freq);
	meiosis_nomut_ = meiosis_matrix(p.mu,p.nuc_freq,true);

	mitosis_ = mitosis_matrix(p.mu_somatic,p.nuc_freq);
	pcr_ = mitosis_matrix(p.mu_pcr,p.nuc_freq);

	return true;
}

/*
RULES FOR LINKING READ GROUPS TO PEOPLE.

0) Read all read groups in bam files. Map RG to Library to Sample.

1) If tissue information is present, build a tissue tree connecting samples
   to zygotic genotypes.  Build mutation matrices for every unique branch
   length.

2) If no tissue information is present in pedigree, check to see if there is a
   sample with the same label as the individual.  If so, connect that sample to
   the individual with a somatic branch of length 0.  [TODO: Length 0 or 1???]

3) If a sample has multiple libraries, connect the libraries to sample with a 
   pcr-error branch.  Otherwise, connect a single library directly.

4) If a library has multiple read-groups, concat the read-groups.

*/

bool dng::Pedigree::Construct(const io::Pedigree& pedigree, const dng::ReadGroups& rgs) {
	using namespace boost;
	using namespace std;
	
  	// Reset Family Information
  	roots_.clear();
  	family_members_.clear();
  	family_members_.reserve(128);
  	peeling_op_.clear();
  	peeling_op_.reserve(128);
			
	num_members_ = pedigree.member_count();
	num_libraries_ = rgs.libraries().size();

	Graph pedigree_graph(num_members_+num_libraries_);
	graph_traits<Graph>::edge_iterator ei, ei_end;
	graph_traits<Graph>::vertex_iterator vi, vi_end;
	
	auto edge_types = get(edge_type, pedigree_graph);
	auto labels = get(vertex_label, pedigree_graph);
	auto groups = get(vertex_group, pedigree_graph);
	auto families = get(edge_family, pedigree_graph);

	// Go through rows and construct the graph
	vertex_t dummy_index = 0;
	for(auto &row : pedigree.table()) {
		// TODO: check for parent-child inbreeding
		vertex_t child = row.child;
		vertex_t dad = row.dad;
		vertex_t mom = row.mom;

		// check to see if mom and dad have been seen before
		auto id = edge(dad, mom, pedigree_graph);
		if(!id.second)
			add_edge(dad, mom, EdgeType::Spousal, pedigree_graph);
		// add the meiotic edges
		add_edge(mom, child, EdgeType::Meiotic, pedigree_graph);
		add_edge(dad, child, EdgeType::Meiotic, pedigree_graph);

		// Process newick file
		int res = newick::parse(row.sample_tree, child, pedigree_graph);
		if(res == -1) {
			throw std::runtime_error(
				"unable to parse somatic data for individual '" + pedigree.name(row.child) + "'." );
		} else if(res == 0) {
			vertex_t v = add_vertex(pedigree.name(row.child), pedigree_graph);
			add_edge(child,v,EdgeType::Mitotic,pedigree_graph);
		}
	}
	// Remove the dummy individual from the graph
	clear_vertex(dummy_index, pedigree_graph);

	// Connect Samples to Libraries
	for(tie(vi,vi_end) = vertices(pedigree_graph); vi != vi_end; ++vi) {
		if(labels[*vi].empty())
			continue;
		auto r = rgs.data().get<rg::sm>().equal_range(labels[*vi]);
		for(;r.first != r.second;++r.first) {
			vertex_t v = num_members_+rg::index(rgs.libraries(), r.first->library);
			if(!edge(*vi,v,pedigree_graph).second)
				add_edge(*vi,v,EdgeType::Library,pedigree_graph);
		}
		labels[*vi].clear();
	}
	
	// Calculate the connected components.  This defines independent sections
	// of the graph.
	std::size_t num_groups = connected_components(pedigree_graph, groups);
	
	// Calculate the biconnected components and articulation points.
	// This defines "nuclear" families and pivot indiviuduals.
	// Nodes which have no edges will not be part of any family.
	vector<vertex_t> articulation_vertices;
	std::size_t num_families = biconnected_components(pedigree_graph, families,
		back_inserter(articulation_vertices)).first;
		
	// Determine which edges belong to which nuclear families.
	typedef vector<vector<graph_traits<Graph>::edge_descriptor>> 
		family_labels_t;
	family_labels_t family_labels(num_families);
  	for(tie(ei, ei_end) = edges(pedigree_graph); ei != ei_end; ++ei)
  		family_labels[families[*ei]].push_back(*ei);
	
	// Determine the last family in each group.  All singleton groups will have
	// a value of -1 since they have no family assignment.
	typedef deque<std::size_t> root_families_t;
	root_families_t root_families(num_groups,-1);  	
  	for(std::size_t f = 0; f < family_labels.size(); ++f) {
  		// last one wins
  		auto first_edge = family_labels[f][0];
  		auto src_vertex = source(first_edge, pedigree_graph);
		root_families[groups[src_vertex]] = f;
  	}
  	
  	// Identify the pivot for each family.
  	// The pivot will be the last art. point that has an edge in
  	// the group.  The pivot of the last group doesn't matter.
  	vector<vertex_t> pivots(num_families,dummy_index);
  	for(auto a : articulation_vertices) {
  		graph_traits<Graph>::out_edge_iterator ei, ei_end;
		for(tie(ei, ei_end) = out_edges(a, pedigree_graph); ei != ei_end; ++ei) {
			// Just overwrite existing value so that the last one wins.
			pivots[families[*ei]] = a;
		}
  	}
  	// Root Pivots are special
  	for(auto f : root_families) {
  		if(f == -1) // skip all singleton groups
  			continue;
  		pivots[f] = dummy_index;
  	}

  	/*
	// Non-dummy singleton groups are root elements.
	// This has been removed because singleton groups have no data
	for(tie(vi, vi_end) = vertices(pedigree_graph); vi != vi_end; ++vi) {
		if(degree(*vi,pedigree_graph) > 0 || *vi == dummy_index)
			continue;
		roots_.push_back(*vi);
	}
	*/

/*	for(tie(ei, ei_end) = edges(pedigree_graph); ei != ei_end; ++ei) {
		cerr << "[" << (int)edge_types[*ei] << ","
			<< "" << families[*ei]  << "] "
			<< source(*ei,pedigree_graph) << " -> " << target(*ei,pedigree_graph)
			<< " " << labels[target(*ei,pedigree_graph)]
			<< "\n";
	}

 	for(std::size_t k = 0; k < family_labels.size(); ++k) {
 		cerr << "Family " << k+1 << ": " << pivots[k] << "\n";
 		for(auto &a : family_labels[k]) {
			cerr << "    [" << (int)edge_types[a] << "] "
				<< source(a,pedigree_graph) << " -> " << target(a,pedigree_graph)
				<< " " << labels[target(a,pedigree_graph)]
				<< "\n";
 		}
 	}
 	cerr << "Roots:";
 	for(auto a : roots_)
 		cerr << " " << a;
 	cerr << endl;
*/
 	num_nodes_ = num_vertices(pedigree_graph);
 	upper_.assign(num_nodes_, DNG_INDIVIDUAL_BUFFER_ASSIGN_TYPE);
 	lower_.assign(num_nodes_, DNG_INDIVIDUAL_BUFFER_ASSIGN_TYPE);
 	buffer_.resize(100,1);
 	full_transition_matrices_.resize(num_nodes_);
 	nomut_transition_matrices_.resize(num_nodes_);

 	vector<std::size_t> lower_written(num_nodes_,-1);

  	// Detect Family Structure and pivot positions
  	for(std::size_t k = 0; k < family_labels.size(); ++k) {
  		auto &family_edges = family_labels[k];
  		// Sort edges based on type and target
  		boost::sort(family_edges, [&](edge_t x, edge_t y) -> bool { return
  			(edge_types(x) < edge_types(y)) &&
  			(target(x, pedigree_graph) < target(y, pedigree_graph)); });
  		
  		// Find the range of the parent types
  		auto pos = boost::find_if(family_edges, [&](edge_t x) -> bool { 
  			return (edge_types(x) != EdgeType::Spousal); });
  		size_t num_parent_edges = distance(family_edges.begin(),pos);
  		
  		// Check to see what type of graph we have
  		if(num_parent_edges == 0) {
  			// If we do not have a parent-child single branch,
  			// we can't construct the pedigree.
  			// TODO: Write Error message
  			if(family_edges.size() != 1)
  				return false;
  			vertex_t parent = source(*pos, pedigree_graph);
  			vertex_t child = target(*pos, pedigree_graph);
  			family_members_.push_back({parent,child});

  			if(edge_types(*pos) == EdgeType::Library) {
				full_transition_matrices_[child] = pcr_;
				nomut_transition_matrices_[child] = pcr_;
  			} else {
				full_transition_matrices_[child] = mitosis_;
				nomut_transition_matrices_[child] = mitosis_;
			}

  			if(pivots[k] == dummy_index) {
  				roots_.push_back(family_members_.back()[0]);
  				if(lower_written[parent] != -1)
  					peeling_op_.emplace_back(&Op::PeelUp2);
  				else {
  					lower_written[parent] = k;
  					peeling_op_.emplace_back(&Op::PeelUp);  					
  				}
  			} else if(pivots[k] == parent) {
  				if(lower_written[parent] != -1)
  					peeling_op_.emplace_back(&Op::PeelUp2);
  				else {
  					lower_written[parent] = k;
  					peeling_op_.emplace_back(&Op::PeelUp);  					
  				}
  			} else {
  				peeling_op_.emplace_back(&Op::PeelDown);
  			}
  		} else if(num_parent_edges == 1) {
  			// We have a nuclear family with 1 or more children
  			family_members_.push_back({
  				source(family_edges.front(), pedigree_graph), // Dad
  				target(family_edges.front(), pedigree_graph)  // Mom
  			});
  			auto& family_members = family_members_.back();
  			while(pos != family_edges.end()) {
  				vertex_t child = target(*pos, pedigree_graph);
	 			full_transition_matrices_[child] = meiosis_;
				nomut_transition_matrices_[child] = meiosis_nomut_; 				
  				family_members.push_back(child); // Child
  				++pos; ++pos; // child edges come in pairs
  			}
  			auto pivot_pos = boost::range::find(family_members, pivots[k]);
  			size_t p = distance(family_members.begin(),pivot_pos);
  			// A family without a pivot is a root family
  			if(pivots[k] == dummy_index ) {
  				p = 0;
  				roots_.push_back(family_members[0]);
  			} else if(p > 2) {
  				swap(family_members[p], family_members[2]);
  				p = 2;
  			}
  			// TODO: Assert that p makes sense
			
			switch(p) {
			case 0:
				if(lower_written[family_members[0]] != -1) {
					peeling_op_.emplace_back(&Op::PeelToFather2);	
				} else {
					lower_written[family_members[0]] = k;
					peeling_op_.emplace_back(&Op::PeelToFather);
				}
				break;
			case 1:
				if(lower_written[family_members[1]] != -1) {
					peeling_op_.emplace_back(&Op::PeelToMother2);	
				} else {
					lower_written[family_members[1]] = k;
					peeling_op_.emplace_back(&Op::PeelToMother);
				}
				break;
			case 2:
			default:
				peeling_op_.emplace_back((family_members.size() == 3) ?
					&Op::PeelToChild : &Op::PeelToChild2);
				break;
			};
  		} else {
  			// Not a zero-loop pedigree
  			// TODO: write error message
  			return false;
  		}
   	}

/*	
	cerr << "Init Op\n";
 	for(int i=0;i<num_members_;++i) {
 		cerr << "\tw\tupper[" << i << "]\n"; 
 	}
 	for(int i=0;i<num_libraries_;++i) {
 		cerr << "\tw\tlower[" << num_members_+i << "]\n"; 
 	}
 	for(int i=0;i<peeling_op_.size();++i) {
 		cerr << "Peeling Op " << i+1;
 		if(peeling_op_[i] == &Op::PeelUp) {
 				cerr << " (PeelUp)\n";
 				cerr << "\tw\tlower[" << family_members_[i][0] << "]\n";
 				cerr << "\tr\tlower[" << family_members_[i][1] << "]\n";
 		} else if(peeling_op_[i] == &Op::PeelUp2) {
 				cerr << " (PeelUp2)\n";
 				cerr << "\trw\tlower[" << family_members_[i][0] << "]\n";
 				cerr << "\tr\tlower[" << family_members_[i][1] << "]\n";
 		} else if(peeling_op_[i] == &Op::PeelDown) {
 				cerr << " (PeelDown)\n";
 				cerr << "\tw\tupper[" << family_members_[i][1] << "]\n";
 				cerr << "\tr\tupper[" << family_members_[i][0] << "]\n";
 				cerr << "\tr\tlower[" << family_members_[i][0] << "]\n";
 		} else if(peeling_op_[i] == &Op::PeelToFather) {
 				cerr << " (PeelToFather)\n";
 				cerr << "\tw\tlower[" << family_members_[i][0] << "]\n";
 				cerr << "\tr\tupper[" << family_members_[i][1] << "]\n";
 				cerr << "\tr\tlower[" << family_members_[i][1] << "]\n";
 				for(int j=2;j<family_members_[i].size();++j)
	 				cerr << "\tr\tlower[" << family_members_[i][j] << "]\n";
 		} else if(peeling_op_[i] == &Op::PeelToFather2) {
 				cerr << " (PeelToFather2)\n";
 				cerr << "\trw\tlower[" << family_members_[i][0] << "]\n";
 				cerr << "\tr\tupper[" << family_members_[i][1] << "]\n";
 				cerr << "\tr\tlower[" << family_members_[i][1] << "]\n";
 				for(int j=2;j<family_members_[i].size();++j)
	 				cerr << "\tr\tlower[" << family_members_[i][j] << "]\n";
 		} else if(peeling_op_[i] == &Op::PeelToMother) {
 				cerr << " (PeelToMother)\n";
 				cerr << "\tw\tlower[" << family_members_[i][1] << "]\n";
 				cerr << "\tr\tupper[" << family_members_[i][0] << "]\n";
 				cerr << "\tr\tlower[" << family_members_[i][0] << "]\n";
 				for(int j=2;j<family_members_[i].size();++j)
	 				cerr << "\tr\tlower[" << family_members_[i][j] << "]\n";
 		} else if(peeling_op_[i] == &Op::PeelToMother2) {
 				cerr << " (PeelToMother2)\n";
 				cerr << "\trw\tlower[" << family_members_[i][1] << "]\n";
 				cerr << "\tr\tupper[" << family_members_[i][0] << "]\n";
 				cerr << "\tr\tlower[" << family_members_[i][0] << "]\n";
 				for(int j=2;j<family_members_[i].size();++j)
	 				cerr << "\tr\tlower[" << family_members_[i][j] << "]\n";
 		} else if(peeling_op_[i] == &Op::PeelToChild) {
 				cerr << " (PeelToChild)\n";
 				cerr << "\tw\tupper[" << family_members_[i][2] << "]\n";
 				cerr << "\tr\tupper[" << family_members_[i][0] << "]\n";
 				cerr << "\tr\tlower[" << family_members_[i][0] << "]\n";
 				cerr << "\tr\tupper[" << family_members_[i][1] << "]\n";
 				cerr << "\tr\tlower[" << family_members_[i][1] << "]\n";
 		} else if(peeling_op_[i] == &Op::PeelToChild2) {
 				cerr << " (PeelToChild2)\n";
 				cerr << "\tw\tupper[" << family_members_[i][2] << "]\n";
 				cerr << "\tr\tupper[" << family_members_[i][0] << "]\n";
 				cerr << "\tr\tlower[" << family_members_[i][0] << "]\n";
 				cerr << "\tr\tupper[" << family_members_[i][1] << "]\n";
 				cerr << "\tr\tlower[" << family_members_[i][1] << "]\n";
 				for(int j=3;j<family_members_[i].size();++j)
	 				cerr << "\tr\tlower[" << family_members_[i][j] << "]\n";
  		} else
 				cerr << " (Unknown)\n";
 	}
	cerr << "Root Op\n";
	for(int i=0;i<roots_.size();++i) {
		cerr << "\tr\tupper[" << roots_[i] << "]\n";
		cerr << "\tr\tlower[" << roots_[i] << "]\n";
	} 	
*/
	return true;
}


