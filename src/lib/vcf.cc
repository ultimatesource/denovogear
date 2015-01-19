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

#include <iostream>
#include <string>
#include <boost/lexical_cast.hpp>
#include "dng/hts/vcf.h"
#include "dng/seq.h"

using namespace std;
using namespace dng;

namespace hts { namespace vcf {

    Vcfostream::Vcfostream(const char *file, const char *mode, const char *pedfile) : hts::File(file, mode)
    {
      hdr = bcf_hdr_init("w"); // VCF header
      rec = bcf_init1(); // VCF record/body line

      bcf_hdr_append(hdr, "##source=" SOURCE);
    }
    
    void Vcfostream::addHdrMetadata(const char *key, const char *value)
    {
      // convert key and value cstrings into one cstring "##key=value"
      size_t key_len = sizeof(key);
      size_t val_len = sizeof(value);
      char line[key_len + val_len + 3]; //
      strcpy(line, "##");
      strcat(line, key);
      strcat(line, "=");
      strcat(line, value);
      bcf_hdr_append(hdr, line);
    }

    void Vcfostream::addHdrMetadata(const char *key, double value)
    {
      string valstr = boost::lexical_cast<string>(value);
      addHdrMetadata(key, valstr.c_str());
    }

    void Vcfostream::addHdrMetadata(const char *key, int value)
    {
      string valstr = boost::lexical_cast<string>(value);
      addHdrMetadata(key, valstr.c_str());
    }

    void Vcfostream::addSample(const char *sample)
    {
      bcf_hdr_add_sample(hdr, sample);
      nsamples++;
    }

    void Vcfostream::writeHdr()
    {
      bcf_hdr_set_version(hdr, VCF_VERSION);
      bcf_hdr_append(hdr, "##INFO=<ID=LL,Number=1,Type=Float,Description=\"Log likelihood\">");
      bcf_hdr_append(hdr, "##INFO=<ID=PMUT,Number=1,Type=Float,Description=\"Probability of mutation\">");
      //bcf_hdr_append(hdr, "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">");
      //bcf_hdr_append(hdr, "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">");
      //bcf_hdr_append(hdr, "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">");
      //bcf_hdr_append(hdr, "##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">");
      //bcf_hdr_append(hdr, "##INFO=<ID=DB,Number=0,Type=Flag,Description=\"dbSNP membership, build 129\">");
      //bcf_hdr_append(hdr, "##INFO=<ID=H2,Number=0,Type=Flag,Description=\"HapMap2 membership\">");
      bcf_hdr_append(hdr, "##FILTER=<ID=PASS,Description=\"All filters passed\">");
      //bcf_hdr_append(hdr, "##FILTER=<ID=q10,Description=\"Quality below 10\">");
      //bcf_hdr_append(hdr, "##FILTER=<ID=s50,Description=\"Less than 50% of samples have data\">");
      //bcf_hdr_append(hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
      //bcf_hdr_append(hdr, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">");
      //bcf_hdr_append(hdr, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">");
      //bcf_hdr_append(hdr, "##FORMAT=<ID=HQ,Number=2,Type=Integer,Description=\"Haplotype Quality\">");

      // AD defined http://gatkforums.broadinstitute.org/discussion/1268/how-should-i-interpret-vcf-files-produced-by-the-gatk
      bcf_hdr_append(hdr, "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">");
  
      bcf_hdr_add_sample(hdr, NULL); // libhts requires NULL sample before it will write out all the other samples
      bcf_hdr_write(handle(), hdr);

    }
    
    void Vcfostream::addRecord(const char *chrom, int pos, const char ref, float ll, float pmut, std::vector<dng::depth5_t> &read_depths)
    {
      // CHROM
      // there needs to be a matching contig id in the header otherwise bcf_write1() will fault [tracked issue to vcf.cc vcf_format()]
      string chrom_s(chrom);
      if(contig_ids.find(chrom_s) == contig_ids.end())	{
	  stringstream conv;
	  conv << "##contig=<ID=" << chrom_s << ",length=1>"; // length is required by htslib but not sure what is the appropiate value?
	  bcf_hdr_append(hdr, conv.str().c_str());
	  //bcf_hdr_write(fp,hdr); //TODO: unless we call write the ##contig header won't show up in the output
	  contig_ids.insert(chrom_s);
	}

      rec->rid = bcf_hdr_name2id(hdr, chrom);
      
      // POS
      rec->pos = pos;

      // ID
      // None given, hts will default to '.'

      // QUAL
      // Currently newcaller does not display site quality nor collect the QUAL field from the bam file, skipping

      //FILTER
      // In newcaller site either PASSES min_basequal/min_prob or it doesn't. In the future consider adding options for sites that 
      // don't pass but should still be displayed.
      int32_t filter_array = bcf_hdr_id2int(hdr, BCF_DT_ID, "PASS");
      bcf_update_filter(hdr, rec, & filter_array, 1);


      // INFO
      // TODO: consider separating INFO into its own method to allow for abirtrary site 
      bcf_update_info_float(hdr, rec, "LL", &ll, 1);
      bcf_update_info_float(hdr, rec, "PMUT", &pmut, 1);
      
      // REF, ALT and genotype samples
      //first determine how many alleles (reference + alternatives) exists
      uint16_t order[4]; // Bases have to be listed in order they appear in both REF and ALT column
      size_t nalleles = 1;
      size_t ref_allele = seq::char_index(ref); 
      order[0] = ref_allele;
      order[1] = (ref_allele+1)%4;
      order[2] = (ref_allele+2)%4;
      order[3] = (ref_allele+3)%4;
      // for each nucleotide (other than the ref) see if the NT exists in at-least one of the samples, add to order[] map
      for(int index = 1; index < 4; index++)
	{
	  int nt = order[index];
	  for(int sample = 0; sample < read_depths.size(); sample++)
	    {
	      if(read_depths[sample].counts[nt] != 0)
		{
		  // Base exists, add it to the ordered list of bases and move on to looking for next base
		  order[nalleles++] = nt;
		  break;
		}
	    }
	}

      // Build a comma-seperated allele string used for REF/ALT field (first member of list is used as REF).
      const char *allelestr[nalleles];
      stringstream ss;
      ss << ref; // REF column
      string tmp;
      for(int a = 1; a < nalleles; a++)
	{
	  char nt = seq::indexed_char(order[a]); // ALT column
	  ss << "," << nt;
	}
      bcf_update_alleles_str(hdr, rec, ss.str().c_str());
      
      // Create the data for the N samples, array will be split into N fields and lists how many bases matched
      // REF and ALT alleles (in the given order)
      int32_t *gtcounts = new int32_t[bcf_hdr_nsamples(hdr)*nalleles];
      size_t index = 0;
      for(int sample = 0; sample < read_depths.size(); sample++)
	{
	  for(int nt = 0; nt < nalleles; nt++)
	    {
	      //cout << read_depths[sample].counts[order[nt]] << endl;
	      gtcounts[index++] = read_depths[sample].counts[order[nt]];
	    }
	} 
      bcf_update_format_int32(hdr, rec, "AD", gtcounts, bcf_hdr_nsamples(hdr)*nalleles);
      delete [] gtcounts;


      // Add line to the body
      bcf_write1(handle(), hdr, rec);

      // reset the record for the next line
      bcf_clear(rec);
    }
    
    Vcfostream::~Vcfostream()
    {
      if(rec != NULL)
	bcf_destroy1(rec);
      
      if(hdr != NULL)
	bcf_hdr_destroy(hdr);
      
      // hts_close(fp) is called in Base class destructor. This means stream has to be destroyed before it 
      // will be written to a file. A close() function may be preferable.

    }


}} // hts::vcf




/*
void makecounts(vector<dng::depth5_t> &read_depths, uint16_t *vals, size_t s)
{
  int rd_indx = 0;
  for(int a = 0; a < s; a++)
    {
      read_depths[rd_indx].counts[a%4] = vals[a];
      if(a != 0 && a%4 == 0)
	rd_indx++;
      //cout << "[" << rd_indx << ", " << a%4 << "] = " << vals[a] << endl;
    }

}

int main(int argc, char *argv[])
{
  hts::vcf::Vcfostream output("test5.vcf", "w", "ceu.ped");
  output.addHdrMetadata("qlen", 10);
  output.addHdrMetadata("region", "XX");
  output.addHdrMetadata("pop_diversity", 0.001);
  output.addSample("Solexa-18483.NA12878");
  output.addSample("Solexa-18484.NA12878");	
  output.addSample("Solexa-56411.NA12892");
  output.addSample("Solexa-56414.NA12891");	
  output.addSample("Solexa-59964.NA12891");
  output.addSample("Solexa-59965.NA12892");
  output.addSample("Solexa-59968.NA12891");
  output.addSample("Solexa-59969.NA12892");
  output.addSample("Solexa-61054.NA12892");
  output.addSample("Solexa-61060.NA12891");
  output.writeHdr();
  
  std::vector<depth5_t> read_depths(10,{0,0});
  uint16_t tmps[] = {0,0,11,27,0,0,19,21,0,0,0,18,0,0,0,18,0,0,0,13,0,0,0,10,0,0,0,11,0,0,0,28,0,0,0,32,0,0,0,32};
  //cout << endl <<  sizeof(tmps)/sizeof(uint16_t) << endl;
  makecounts(read_depths, tmps, sizeof(tmps)/sizeof(uint16_t));
  //0,0,11,27	0,0,19,21	0,0,0,18	0,0,0,18	0,0,0,13	0,0,0,10	0,0,0,11	0,0,0,28	0,0,0,32	0,0,0,32

  output.addRecord("10", 204931, 'G', -93.3296, 1, read_depths);

  //A = 0
  //C = 1
  //T = 3
  //G = 2

  //cout << "Hello" << endl;
}
*/
