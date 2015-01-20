/*
 * Copyright (c) 2015 Reed A. Cartwright
 * Authors:  kdai1 <kdai1@asu.edu>
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
#ifndef CXX_HTS_VCF_H
#define CXX_HTS_VCF_H

#include "dng/hts/hts.h"
#include "dng/matrix.h"

#include <vector>
#include <string>
#include <set>
#include <htslib/vcf.h>

#define VCF_VERSION "VCFv4.1"
#define SOURCE "Denovogear 1.0.1"

namespace hts { namespace vcf { 

/**
 * File - Class for writing VCF/BCF files to the output, either stdout or to a file
 * To write to stdout set file="-" in constructor. To write as a BCF file set mode="wb" in constructor.
 * To use:
 *  1. call addHdrMetadata() and addSample() to add information to the VCF header
 *  2. call writeHdr() to finish the header
 *  3. For each site that will be output call addRecord()
 *  4. The destructor needs to be called to write to the file.
 */
class File : hts::File
{
 public:
  File(const char *file, const char *mode, const char *pedfile);
  File(const File&) = delete;
  File& operator=(File&) = delete;
  virtual ~File();

  /** AddHeaderMetadata() - Adds a "##key=value" line to the VCF header */
  void AddHeaderMetadata(const char *key, const char *value);
  void AddHeaderMetadata(const char *key, double value);
  void AddHeaderMetadata(const char *key, int value);

  /** Add another sample/genotype field to the VCF file. */
  void AddSample(const char *sample);

  /** Creates the header field. Call only after adding all the sample fields */
  void WriteHeader();

  /**
   * AddRecord() - Add a new line to the body of the VCF/BCF file
   * @chrom: Value in CHROM field, the contig id for the sample
   * @pos: POS filed
   * @ref: REF field
   * @ll: log-likelihood for site, used in INFO field
   * @pmut: Probability of mutation, used in INFO field
   * @read_depths: contains information about site genotype for each sample
   * 
   * Builds a single line of data in the body of the VCF file. Should be called after writeHdr(). The ALT and genotype
   * fields build from read_depths.count[]. Each genotype field will list the number of bases in the sample - listed in 
   * the order the bases appear in the REF and ALT field.
   */
  void AddRecord(const char *chrom, int pos, const char ref, float ll, float pmut, std::vector<dng::depth5_t> &read_depths);

  
						       
 private:
  bcf_hdr_t *hdr;
  bcf1_t *rec;

  std::size_t nsamples = 0;
  std::set<std::string> contig_ids; // List of unique CHROM/contig values

};
  

}}

#endif /* CXX_HTS_VCF_H */
