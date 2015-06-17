/*
 * Copyright (c) 2010, 2011 Genome Research Ltd.
 * Copyright (c) 2012, 2013 Donald Conrad and Washington University in St. Louis
 * Authors: Donald Conrad <dconrad@genetics.wustl.edu>,
 * Avinash Ramu <aramu@genetics.wustl.edu>
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

#ifndef LOOKUP_H_
#define LOOKUP_H_

#include <math.h>
#define WANT_STREAM       // include iostream and iomanipulators
//#include "newmatap.h"
//#include "newmatio.h"
#include <Eigen/Dense>


#define MIN_READ_DEPTH_SNP 10

// The old code was mixing floats and doubles when setting and processing matrices. Found it caused problems 
// when max/min was called on matrices (*Like.cc) - had multiple maxs/mins due to roundoff.
typedef double Real;

// Matrices for parameters
typedef Eigen::Matrix<Real, 100, 10> SNPMatrix;
typedef Eigen::Matrix<Real, 9, 3> IndelMatrix;
typedef Eigen::Matrix<Real, 10, 10> PairMatrix;
// TODO: A hack to save time. Should replace Dynamic size with appropiate sizes. Used in the *Like.cc methods.
typedef Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> Matrix;

// Used to convert array to Eigen matrices in makeLookup.cc
typedef Eigen::Map<Eigen::Matrix<Real, 100, 10, Eigen::RowMajor> > mapSNPMatrix;
typedef Eigen::Map<Eigen::Matrix<Real, 9, 3, Eigen::RowMajor> > mapIndelMatrix;
typedef Eigen::Map<Eigen::Matrix<Real, 10, 10, Eigen::RowMajor> > mapPairMatrix;

typedef std::vector<std::vector<std::string> > lookup_table_t;

// SNP Lookup Table
typedef struct {
  SNPMatrix aref; /*priors for "A" reference allele */
  SNPMatrix cref; /*priors for "C" reference allele */
  SNPMatrix gref; /*priors for "G" reference allele */
  SNPMatrix tref; /*priors for "T" reference allele */
  SNPMatrix snpcode; /*code for identifying a biallelic locus */
  SNPMatrix tp; /*code for identifying a biallelic locus */
  SNPMatrix code; /*code for identifying a biallelic locus */
  SNPMatrix mrate; /*mutation probability*/
  SNPMatrix denovo; /*code for identifying de novo events + mutation rate*/
  SNPMatrix norm; /*code for identifying de novo events + mutation rate*/
} lookup_snp_t;

// Indel Lookup Table
// TODO: mrate is not always being filled out in makeLookups.cc. Investigate
typedef struct {
  IndelMatrix priors; /*priors for Indel */
  IndelMatrix snpcode; /*code for identifying a biallelic locus */
  IndelMatrix tp; /*code for identifying a biallelic locus */
  IndelMatrix code; /*code for identifying a biallelic locus */
  IndelMatrix mrate; /*mutation probability*/
  IndelMatrix denovo; /*code for identifying de novo events + mutation rate*/
  IndelMatrix norm; /*code for identifying de novo events + mutation rate*/
  IndelMatrix hit; /* multiplication factor for indel mu */
} lookup_indel_t;

// Paired-sample Lookup Table
typedef struct {
  PairMatrix snpcode; /*code for identifying a biallelic locus */
  PairMatrix priors; /*priors */
  PairMatrix denovo; /*code for identifying de novo events + mutation rate*/
  PairMatrix norm; /*code for identifying de novo events + mutation rate*/
} lookup_pair_t;

#endif
