/*
 * Copyright (c) 2015 Reed A. Cartwright <reed@cartwrig.ht>
 * Copyright (c) 2015 Kael Dai <kdai1@asu.edu>
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
#ifndef CXX_HTS_BCF_H
#define CXX_HTS_BCF_H

#include "hts.h"

#include <vector>
#include <string>
#include <set>
#include <htslib/vcf.h>

namespace hts { namespace bcf { 

/**
 * File - Class for writing VCF/BCF files to the output, either stdout or to a file
 * To write to stdout set file="-" in constructor. To write as a BCF file set mode="wb" in constructor.
 * To use:
 *  1. call AddHeaderMetadata() and AddSample() to add information to the VCF header
 *  2. call WriteHdr() to finish the header
 *  3. populate the first few columns of a VCF body line using SetID(), SetFilter(), SetAlleles(), UpdateInfo() and UpdateSamples() 
 *  4. call WriteRecord() to save/print the current VCF line and move on to the next
 */
class File : public hts::File {
public:
  
  /**
   * Initialize File/output stream for writing VCF/BCF
   * @file - file name to write to, use "-" for stdout
   * @mode: "w" = write VCF, "wb" = write BCF
   * @source: name of application writing the VCF file (optional)
   */
	File(const char *file, const char *mode, const char *source = nullptr) : hts::File(file, mode) {
		hdr = bcf_hdr_init("w"); // VCF header
		rec = bcf_init1(); // VCF record/body line

		if(source != nullptr) {
			std::string s = std::string("#source=\"") + source + "\"";
			bcf_hdr_append(hdr, s.c_str());
		}
    }
    File(const File&) = delete;
	File& operator=(const File&) = delete;

	File(File&& other) : hts::File(std::move(other)), hdr(other.hdr),
	  rec(other.rec)//, contig_ids(std::move(other.contig_ids))
    {
		other.hdr = nullptr;
		other.rec = nullptr;
    }


	File& operator=(File&& other) {
		if(this == &other)
      		return *this;
      
		hdr = other.hdr;
		other.hdr = nullptr;

		rec = other.rec;
		other.rec = nullptr;

		//contig_ids = std::move(other.contig_ids);

		return *this;
	}
  
	virtual ~File() {
		if(rec != nullptr)
			bcf_destroy1(rec);

		if(hdr != nullptr)
			bcf_hdr_destroy(hdr);
    }

	/** Use this version to add a new INFO/FILTER/FORMAT string to the header */
	int AddHeaderMetadata(const char *line) {
		if(line == nullptr)
			return -1; //bcf_hdr_append returns if line cannot be processed
		return bcf_hdr_append(hdr, line);
	}

	/** AddHeaderMetadata() - Adds a "##key=value" line to the VCF header */
	int AddHeaderMetadata(const char *key, const char *value) {
		if(key == nullptr || value == nullptr)
			return -1;
		std::string line = std::string("##") + key + "=" + value;
		return bcf_hdr_append(hdr, line.c_str());
	}

	int AddHeaderMetadata(const char *key, const std::string& value) {
		return AddHeaderMetadata(key, value.c_str());
	}


	template<typename T>
	int AddHeaderMetadata(const char *key, T value) {
		return AddHeaderMetadata(key, std::to_string(value));
	}

	/** Add another sample/genotype field to the VCF file. */
	int AddSample(const char *sample) {
		return bcf_hdr_add_sample(hdr, sample);
	}

	/** Add a "#contig=" metadata entry to the VCF header. */
	int AddContig(const char *contigid) {
	        if(contigid == nullptr)
		        return -1;
	        std::string conv = std::string("##contig=<ID=") + contigid + ",length=1>";
	        return bcf_hdr_append(hdr, conv.c_str());
	}

	/** Creates the header field. Call only after adding all the sample fields */
	int WriteHeader() {
		bcf_hdr_add_sample(hdr, nullptr); // libhts requires NULL sample before it will write out all the other samples
		return bcf_hdr_write(handle(), hdr);
	}


	//TODO: Split off rec into its own wrapper.

  /** 
   * SetID() - Sets the CHROM, POS, and ID fields in the VCF record.
   * chrom is required field. If id is NULL then ID value will default to "."
   */
	void SetID(const char *chrom, int pos, const char *id = nullptr) {
		// CHROM
		if(chrom == nullptr)
			return; // should display some error?

		/*
		std::string chrom_s(chrom);
		/*
		// TODO: Cache last contig_id, since it is likely to be the same??
		// TODO: Preprocess the "##congig strings based on bam contigs???"
		if(contig_ids.find(chrom_s) == contig_ids.end()) {
			// there needs to be a matching contig id in the header otherwise bcf_write1() will fault [tracked issue to vcf.cc vcf_format()]
			std::string conv = std::string("##contig=<ID=") + chrom_s + ",length=1>";
			bcf_hdr_append(hdr, conv.c_str());
			bcf_hdr_write(handle(), hdr);
			contig_ids.insert(chrom_s);
		}
		*/
		rec->rid = bcf_hdr_name2id(hdr, chrom);

		// POS
		rec->pos = pos;

		// ID
		if(id != nullptr) {
			bcf_update_id(hdr, rec, id);
		}
	}

 	void SetQuality(int quality) {
 		rec->qual = quality;
  	}


	int SetFilter(const char *filter) {
		if(filter == nullptr)
			return 0;
		int32_t fid = bcf_hdr_id2int(hdr, BCF_DT_ID, filter);
		return bcf_update_filter(hdr, rec, &fid, 1);
	}

  /**
   * SetAlleles() - Set the REF and ALT fields in the current record
   * @str: A comma separated list of all the alleles that show up in the sample. 
   *       the first element is the REF value.
   */
	int SetAlleles(const std::string &str) {
		if(str.empty())
			return 0;
		return bcf_update_alleles_str(hdr, rec, str.c_str());
	}


  /**
   * UpdateInfoField() - Add another key=Value pair to the INFO field
   * @key: must be defined in the Header or won't be added.
   * @value: if set to NULL then previous set key-value pairs will be removed.
   *
   * TODO: Add a vector<> version of flat, ints, and strings.
   */
	int UpdateInfo(const char *key, float value) {
		return bcf_update_info_float(hdr, rec, key, &value, 1);
	}

	int UpdateInfo(const char *key, int32_t value) {
		return bcf_update_info_int32(hdr, rec, key, &value, 1);
	}

	int UpdateInfo(const char *key, std::string &value) {
		return bcf_update_info_string(hdr, rec, key, value.c_str());
	}


  /**
   * UpdateSample() - Add a key to the FORMAT field and update each sample column
   *                  with the corresponding key values.
   * @name: the tag name that appears in FORMAT field. Must be defined in the header
   * @data: A list of values that will populate the sample/genotype fields. The vector
   *        should be a multiple of the number sample columns. 
   */
	int UpdateSamples(const char *name, const std::vector<float> &data) {
		assert(name != nullptr);
		return bcf_update_format_float(hdr, rec, name, &data[0], data.size());   
	}
	int UpdateSamples(const char *name, const std::vector<int32_t> &data) {
		assert(name != nullptr);
		return bcf_update_format_int32(hdr, rec, name, &data[0], data.size());
	}
	int UpdateSamples(const char *name, const std::vector<const char*> &data) {
		assert(name != nullptr);
		return bcf_update_format_string(hdr, rec, name, const_cast<const char**>(&data[0]), data.size());		
	}
	int UpdateSamples(const char *name, const std::vector<std::string> &data) {
		std::vector<const char*> v(data.size());
		for(decltype(data.size()) u = 0; u < data.size(); ++u )
			v[u] = data[u].c_str();
		return UpdateSamples(name, v);
	}


	/** Writes out the up-to-date info in the record and prepares for the next line */
	void WriteRecord() {
		// Add line to the body of the VCF
		bcf_write1(handle(), hdr, rec);
		// reset the record for the next line
		bcf_clear(rec);
	}
        
protected:
	bcf_hdr_t *hdr;
	bcf1_t *rec;
	
	//std::set<std::string> contig_ids; // List of unique CHROM/contig values
};
  

}}

#endif /* CXX_HTS_BCF_H */
