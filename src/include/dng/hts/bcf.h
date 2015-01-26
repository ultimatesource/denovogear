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
 *  1. call addHeaderMetadata() and addSample() to add information to the VCF header
 *  2. call WriteHdr() to finish the header
 *  3. populate the first few columns of a VCF body line using SetID(), SetFilter(), SetAlleles(), UpdateInfo() and UpdateSamples() 
 *  4. call WriteRecord() to save/print the current VCF line and move on to the next
 */
class File : hts::File
{
 public:
  
  /**
   * Initialize File/output stream for writing VCF/BCF
   * @file - file name to write to, use "-" for stdout
   * @mode: "w" = write VCF, "wb" = write BCF
   * @source: name of application writing the VCF file
   */
  File(const char *file, const char *mode, const char *source)
    : hts::File(file, mode)
    {
      hdr = bcf_hdr_init("w"); // VCF header
      rec = bcf_init1(); // VCF record/body line
      
      if(source != nullptr)
	{
	  // source field is not required but commonly used in VCF files. 
	  // htslib is not displaying line if source has spaces; fixed by wrapping in quotes
	  std::string s = std::string("#source=\"") + source + "\"";
	  bcf_hdr_append(hdr, s.c_str());
	}
    }

 File(File&& other) : hts::File(std::move(other)), 
    hdr(other.hdr), rec(other.rec), contig_ids(other.contig_ids)
    {
      other.hdr = nullptr;
      other.rec = nullptr;
      other.contig_ids.clear();
    }

  File(const File&) = delete;

  File& operator=(const File&) = delete;

  File& operator=(File&& other)
    {
      if(this == &other)
	return *this;
      
      hdr = other.hdr;
      other.hdr = nullptr;

      rec = other.rec;
      other.rec = nullptr;
      
      contig_ids = other.contig_ids;
      other.contig_ids.clear();

      return *this;
    }
      
  virtual ~File()
    {
      if(rec != nullptr)
	bcf_destroy1(rec);

      if(hdr != nullptr)
	bcf_hdr_destroy(hdr);

    }

  /** AddHeaderMetadata() - Adds a "##key=value" line to the VCF header */
  void AddHeaderMetadata(const char *key, const char *value)
  {
    if(key == nullptr || value == nullptr)
      // libhts will seg-fault if line is NULL
      return;

    std::string line = std::string("##") + key + "=" + value;
    bcf_hdr_append(hdr, line.c_str());
  }

  template<typename T>
  void AddHeaderMetadata(const char *key, T value)
  {
    AddHeaderMetadata(key, std::to_string(value).c_str());
  }

  /** Use this version to add a new INFO/FILTER/FORMAT string to the header */
  void AddHeaderMetadata(const char *line)
  {
    if(line != nullptr)
      bcf_hdr_append(hdr, line);
  }

  /** Add another sample/genotype field to the VCF file. */
  void AddSample(const char *sample)
  {
    bcf_hdr_add_sample(hdr, sample);
  }

  /** Creates the header field. Call only after adding all the sample fields */
  void WriteHeader()
  {
    bcf_hdr_add_sample(hdr, nullptr); // libhts requires NULL sample before it will write out all the other samples
    bcf_hdr_write(handle(), hdr);
  }


  /** 
   * SetID() - Sets the CHROM, POS, and ID fields in the VCF record.
   * chrom is required field. If id is NULL then ID value will default to "."
   */
  void SetID(const char *chrom, int pos, const char *id)
  {
    // CHROM
    if(chrom == NULL)
      return; // should display some error?

    std::string chrom_s(chrom);
    if(contig_ids.find(chrom_s) == contig_ids.end())
      {
	// there needs to be a matching contig id in the header otherwise bcf_write1() will fault [tracked issue to vcf.cc vcf_format()]
	std::string conv = std::string("##contig=<ID=") + chrom_s + ",length=1>";
	bcf_hdr_append(hdr, conv.c_str());
	contig_ids.insert(chrom_s);
	//stringstream conv;
	//conv << "##contig=<ID=" << chrom_s << ",length=1>"; // length is required by htslib but not sure what is the appropiate value?
	//bcf_hdr_append(hdr, conv.str().c_str());
	//bcf_hdr_write(fp,hdr); //TODO: unless we call write the ##contig header won't show up in the output
	//contig_ids.insert(chrom_s);
      }
    rec->rid = bcf_hdr_name2id(hdr, chrom);

    // POS
    rec->pos = pos;

    // ID
    if(id != nullptr)
      {
	bcf_update_id(hdr, rec, id);
      }
  }

  void SetQuality(int quality)
  {
    rec->qual = quality;
  }

  void setFilter(const char *filter)
  {
    if(filter != nullptr)
      {
	int32_t fid = bcf_hdr_id2int(hdr, BCF_DT_ID, "PASS");
	bcf_update_filter(hdr, rec, &fid, 1);
      }
  }

  /**
   * SetAlleles() - Set the REF and ALT fields in the current record
   * @alleles: A list of all the alleles that show up in the sample. 
   *          alleles[0] is used as the REF value.
   */
  void SetAlleles(std::vector<std::string> &alleles)
  {
    if(alleles.size() == 0)
      return;

    // convert to a comma separated list 
    std::string ss;
    ss += alleles[0];
    for(std::size_t a = 1; a < alleles.size(); a++)
      ss += std::string(",") + alleles[a];
    
    bcf_update_alleles_str(hdr, rec, ss.c_str());
  }


  /**
   * UpdateInfoField() - Add another key=Value pair to the INFO field
   * @key: must be defined in the Header or won't be added.
   * @value: if set to NULL then previous set key-value pairs will be removed.
   *
   * TODO: Add a vector<> version of flat, ints, and strings.
   */
  void UpdateInfo(const char *key, float value)
  {
    bcf_update_info_float(hdr, rec, key, &value, 1);
  }

  void UpdateInfo(const char *key, int32_t value)
  {
    bcf_update_info_int32(hdr, rec, key, &value, 1);
  }

  void UpdateInfo(const char *key, std::string &value)
  {
    bcf_update_info_string(hdr, rec, key, value.c_str());
  }


  /**
   * UpdateSample() - Add a key to the FORMAT field and update each sample column
   *                  with the corresponding key values.
   * @name: the tag name that appears in FORMAT field. Must be defined in the header
   * @data: A list of values that will populate the sample/genotype fields. The vector
   *        should be a multiple of the number sample columns. 
   */
  void UpdateSamples(const char *name, std::vector<float> &data)
  {
    bcf_update_format_float(hdr, rec, name, &data[0], data.size());   
  }

  void UpdateSamples(const char *name, std::vector<int32_t> &data)
  {
    bcf_update_format_int32(hdr, rec, name, &data[0], data.size());
  }

  
  void UpdateSamples(const char *name, std::vector<std::string> &data)
  {
    const char **values = new const char*[data.size()];
    for(std::size_t a = 0; a < data.size(); a++)
      values[a] = data[a].c_str();
    bcf_update_format_string(hdr, rec, name, values, data.size());

    delete [] values;
  }


  /** Writes out the up-to-date info in the record and prepares for the next line */
  void WriteRecord()
  {
    // Add line to the body of the VCF
    bcf_write1(handle(), hdr, rec);
    
    // reset the record for the next line
    bcf_clear(rec);
  }
  

  
						       
 protected:
  bcf_hdr_t *hdr;
  bcf1_t *rec;

  std::set<std::string> contig_ids; // List of unique CHROM/contig values

};
  

}}

#endif /* CXX_HTS_BCF_H */
