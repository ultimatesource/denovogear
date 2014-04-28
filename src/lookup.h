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

#include <math.h>
#define WANT_STREAM       // include iostream and iomanipulators
#include "newmatap.h"
#include "newmatio.h"

#define MIN_READ_DEPTH_SNP 10


using namespace std;

// SNP Lookup Table
typedef struct {
  Matrix aref; /*priors for "A" reference allele */
  Matrix cref; /*priors for "C" reference allele */
  Matrix gref; /*priors for "G" reference allele */
  Matrix tref; /*priors for "T" reference allele */
  Matrix snpcode; /*code for identifying a biallelic locus */
  Matrix tp; /*code for identifying a biallelic locus */
  Matrix code; /*code for identifying a biallelic locus */
  Matrix mrate; /*mutation probability*/
  Matrix denovo; /*code for identifying de novo events + mutation rate*/
  Matrix norm; /*code for identifying de novo events + mutation rate*/
} lookup_snp_t;

// Indel Lookup Table
typedef struct {
  Matrix priors; /*priors for Indel */
  Matrix snpcode; /*code for identifying a biallelic locus */
  Matrix tp; /*code for identifying a biallelic locus */
  Matrix code; /*code for identifying a biallelic locus */
  Matrix mrate; /*mutation probability*/
  Matrix denovo; /*code for identifying de novo events + mutation rate*/
  Matrix norm; /*code for identifying de novo events + mutation rate*/
  Matrix hit; /* multiplication factor for indel mu */
} lookup_indel_t;

// Paired-sample Lookup Table
typedef struct {
  Matrix snpcode; /*code for identifying a biallelic locus */
  Matrix priors; /*priors */
  Matrix denovo; /*code for identifying de novo events + mutation rate*/
  Matrix norm; /*code for identifying de novo events + mutation rate*/
} lookup_pair_t;
