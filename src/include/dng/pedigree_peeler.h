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

#pragma once
#ifndef DNG_PEDIGREE_PEELER_H
#define DNG_PEDIGREE_PEELER_H

#include <functional>

#include <dng/pedigree.h>
#include <dng/matrix.h>

/*
RULES FOR LINKING READ GROUPS TO PEOPLE.

0) Read all read groups in bam files. Map RG to Library to Sample.

1) If no tissue information is present in pedigree, check to see if there is a
   sample with the same label as the individual.  If so, connect that sample to
   the individual with a somatic branch of length 0.  [TODO: Length 0 or 1???]

2) If tissue information is present, build a tissue tree connecting samples
   to zygotic genotypes.  Build mutation matrices for every unique branch
   length.

3) If a sample has multiple libraries, connect the libraries to sample with a 
   pcr-error branch.  Otherwise, connect a sigle library directly.

4) If a library has multiple read-groups, concat the read-groups.

*/

namespace dng {

class PedigreePeeler {
public:
	typedef std::vector<std::vector<std::size_t>> family_members_t;

	bool Initialize(double theta, double mu);

	bool Construct(const Pedigree& pedigree);
	
	double CalculateLogLikelihood(const IndividualBuffer &penetrances) {
		// TODO: Assert that length of penetrances is equal to num_members_;
		// Copy Penetrance values into the lower buffer
		// TODO: Eliminate this copy????
		lower_ = penetrances;
		// Copy genotype Priors
		// TODO: Only update founders
		// TODO: use a different prior based on reference
		upper_.assign(num_members_, genotype_prior_);
		
		// Peel pedigree one family at a time
		for(std::size_t i = 0; i < peeling_op_.size(); ++i)
			(peeling_op_[i])(*this, i);
			
		// Sum over roots
		double ret = 0.0;
		for(auto r : roots_)
			ret += log((lower_[r]*upper_[r]).sum());
		return ret;
	}
	
protected:
	family_members_t family_members_;
	std::size_t num_members_;
	std::vector<std::size_t> roots_;
	
	IndividualBuffer upper_; // Holds P(Data & G=g)
	IndividualBuffer lower_; // Holds P(Data | G=g)
	
	Vector10d genotype_prior_; // Holds P(G | theta)

	MeiosisMatrix meiosis_;

	void PeelToFather(std::size_t id) {
		using namespace Eigen;
		auto family_members = family_members_[id];
		Vector100d buffer;
		buffer.setOnes();
		// Sum over children
		for(std::size_t i = 2; i < family_members.size(); i++) {
			buffer *= (meiosis_ * lower_[family_members[i]].matrix()).array();
		}
		// Include Mom
		Map<Matrix10d, Aligned> mat(buffer.data());
		lower_[family_members[0]] *= (mat *
			(upper_[family_members[1]]*lower_[family_members[1]]).matrix()).array();
	}

	void PeelToMother(std::size_t id) {
		using namespace Eigen;
		auto family_members = family_members_[id];
		Vector100d buffer;
		buffer.setOnes();
		// Sum over children
		for(std::size_t i = 2; i < family_members.size(); i++) {
			buffer *= (meiosis_ * lower_[family_members[i]].matrix()).array();
		}
		// Include Dad
		Map<RowMatrix10d, Aligned> mat(buffer.data());
		lower_[family_members[1]] *= (mat *
			(upper_[family_members[0]]*lower_[family_members[0]]).matrix()).array();
	}

	void PeelToChild(std::size_t id) {
		using namespace Eigen;
		auto family_members = family_members_[id];
		Vector100d buffer;
		buffer.setOnes();
		// Sum over children
		for(std::size_t i = 3; i < family_members.size(); i++) {
			buffer *= (meiosis_ * lower_[family_members[i]].matrix()).array();
		}
		// Parents
		buffer *= kroneckerProduct(
			(lower_[family_members[0]] * upper_[family_members[0]]).matrix(),
			(lower_[family_members[1]] * upper_[family_members[1]]).matrix()
		).array();
		
		upper_[family_members[2]].matrix() = meiosis_.transpose() * buffer.matrix();
	}

	typedef std::function<void(PedigreePeeler&, std::size_t)> PeelOp;

	std::vector<PeelOp> peeling_op_;	
};

}; // namespace dng


#endif // DNG_PEDIGREE_PEELER_H
