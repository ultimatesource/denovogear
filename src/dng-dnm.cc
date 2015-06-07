#include "dng-dnm.h"
#include "denovogear.h"
#include "pedParser.h"

#include <iostream>
#include <string>
#include <vector>

#include <dng/hts/extra.h>
#include <dng/app.h>
#include <dng/task/dnm.h>
#include "htslib/synced_bcf_reader.h"

using namespace dng::task;


int DNM::operator()(std::string &model, DNM::argument_type &arg) {
  //int callDenovoFromBCF(std::string ped_file, std::string bcf_file, std::string op_vcf_f, 
  //		      std::string model, bool is_vcf, std::string region) {
  if(model != "auto" && model != "XS" && model != "XD")
    throw std::runtime_error("Invalid model option " + model + ". Use auto, XS, or XD.");
  
  if(arg.bcf.empty() || arg.ped.empty())
    throw std::runtime_error("ERROR ! Please specify both the PED file and BCF file ! Exiting!");


  snp_mrate = arg.snp_mrate;
  indel_mrate = arg.indel_mrate;
  poly_rate = arg.poly_rate;
  pair_mrate = arg.pair_mrate;
  mu_scale = arg.mu_scale;
  PP_cutoff = arg.pp_cutoff;
  RD_cutoff = arg.rd_cutoff;

  // Create SNP lookup
  lookup_snp_t lookupSNP;
  lookup_table_t tgtSNP;
  callMakeSNPLookup(tgtSNP, lookupSNP, model);

  // Create INDEL lookup
  lookup_indel_t lookupIndel;
  lookup_table_t tgtIndel;
  callMakeINDELLookup(tgtIndel, lookupIndel, model);

  // Create paired lookup
  lookup_pair_t lookupPair;
  lookup_table_t tgtPair;
  callMakePairedLookup(tgtPair, lookupPair);

  // TODO: Use Reed's PED parser
  // Iterate each position of BCF file
  Trio *trios;
  Pair *pairs;
  int trio_count = 0, pair_count = 0;
  parse_ped(arg.ped, &trios, &pairs, trio_count, pair_count);

  //create output vcf -- assumes theres only one trio/pair
  std::vector<std::string> samples;
  if(trio_count > 0) {
    samples.push_back(trios[0].cID);
    samples.push_back(trios[0].mID);
    samples.push_back(trios[0].fID);
  } 
  else {
    samples.push_back(pairs[0].tumorID);
    samples.push_back(pairs[0].normalID);
  }
   
  std::vector<hts::bcf::File> output_vcf;
  if(!arg.vcf.empty()) {
    output_vcf.emplace_back(arg.vcf.c_str(), "w");
    writeVCFHeader(output_vcf[0], arg.bcf, arg.ped, samples);
  }

  qcall_t mom_snp, dad_snp, child_snp;
  indel_t mom_indel, dad_indel, child_indel;
  pair_t tumor, normal;
  bcf_hdr_t *hout, *hin;
  int tid, begin, end;

  // Parse region
  hts::bcf::File vcf_input(arg.bcf.c_str(), "r");
  if(!arg.region.empty()) {
    tid = begin = end = -1;
    /*
    if(bcf_parse_region(str2id, region.c_str(), &tid, &begin, &end) >= 0) {
      bcf_idx_t *idx;
      idx = bcf_idx_load(bcf_file.c_str());
      if(idx) {
	uint64_t off;
	off = bcf_idx_query(idx, tid, begin);
	if(off == 0) {
	  fprintf(stderr, "[%s] no records in the query region.\n", __func__);
	  return 1; // FIXME: a lot of memory leaks...
	}
	bgzf_seek(bp->fp, off, SEEK_SET);
	bcf_idx_destroy(idx);
      }
    }
    */
  }
  int snp_total_count = 0, snp_pass_count = 0;
  int indel_total_count = 0, indel_pass_count = 0;
  int pair_total_count = 0, pair_pass_count = 0;
  
  const bcf_hdr_t *hdr = vcf_input.header();  
  bcf_srs_t *rec_reader = bcf_sr_init(); // used for iterating each rec in BCF/VCF

  // Open region if specified
  if(!arg.region.empty()) {
    int ret = bcf_sr_set_regions(rec_reader, arg.region.c_str(), 0);
    if(ret == -1)
      throw std::runtime_error("no records in the query region " + arg.region);
  }

  // Initialize the record reader to iterate through the BCF/VCF input
  int ret = bcf_sr_add_reader(rec_reader, arg.bcf.c_str());
  if(ret == 0) {
    int errnum = rec_reader->errnum;
    switch(errnum){
    case not_bgzf:
      throw std::runtime_error("Input file type does not allow for region searchs. Exiting!");
      break;
    case idx_load_failed:
      throw std::runtime_error("Unable to load query region, no index. Exiting!");
      break;
    case file_type_error:
      throw std::runtime_error("Could not load filetype. Exiting!");
      break;
    default:
	throw std::runtime_error("Could not load input sequence file into htslib. Exiting!");
    };
  }

  while(bcf_sr_next_line(rec_reader)) {
    bcf1_t *rec = bcf_sr_get_line(rec_reader, 0);
    int j = 0;
    int flag = 0;
    
    for(j = 0; j < trio_count; j++) {
      bcf_unpack(rec, BCF_UN_STR);
      int is_indel = bcf_2qcall(hdr, rec, trios[j], 
				&mom_snp, &dad_snp, &child_snp,
				&mom_indel, &dad_indel, &child_indel,
				flag);

      if(is_indel == 0) {
	snp_total_count++;
	trio_like_snp(child_snp, mom_snp, dad_snp, flag,
		      tgtSNP, lookupSNP, output_vcf,//vcf_out,
		      PP_cutoff, RD_cutoff, snp_pass_count);
      } else if(is_indel == 1) {
	
	indel_total_count++;
	trio_like_indel(&child_indel, &mom_indel, &dad_indel, flag,
			tgtIndel, lookupIndel, mu_scale, output_vcf,
			PP_cutoff, RD_cutoff, indel_pass_count, indel_mrate);
	
      } else if(is_indel < 0) {
	printf("\n BCF PARSING ERROR - Trios!  %d\n Exiting !\n", is_indel);
	exit(1);
      }
      
      // PROCESS  PAIRS
      if(model == "auto") { // paired sample model not developed for XS, XD yet
	for(j = 0; j < pair_count; j++) {
	  int is_indel = bcf2Paired(hout, rec, pairs[j], &tumor, &normal, flag);
	  if(is_indel == 0) {
	    pair_total_count++;
	    pair_like(tumor, normal, tgtPair, lookupPair, flag, output_vcf,
	    	      PP_cutoff, RD_cutoff, pair_pass_count);
	  } else if(is_indel < 0) {
	    printf("\n BCF PARSING ERROR - Paired Sample!  %d\n Exiting !\n", is_indel);
	    exit(1);
	  }
	}
      }     
    }
  }

  cerr << endl << "Total number of SNP sites interrogated: " << snp_total_count;
  cerr << endl << "Total number of SNP sites passing read-depth filters: " <<
    snp_pass_count;
  cerr << endl << "Total number of INDEL sites interrogated: " <<
    indel_total_count;
  cerr << endl << "Total number of INDEL sites passing read-depth filters: " <<
    indel_pass_count;
  cerr << endl << "Total number of Paired sample sites interrogated: " <<
    pair_total_count;
  cerr << endl <<
    "Total number of Paired sample sites passing read-depth filters: " <<
    pair_pass_count;

  std::cerr << std::endl << "Done !" << std::endl;
  exit(0);
}



/****************************************************************************************
int DNM::operator()(std::string &model, DNM::argument_type &arg) {
  std::cout << "HERE with " << model << std::endl;
}
*/

typedef dng::CommandLineApp<dng::task::DNM> CallApp;

class DNGApp : CallApp {
public:
  std::string model;
  
  DNGApp(int argc, char *argv[]) : CallApp(argc, argv)
  {
    if(argc > 1)
      model = argv[1];
  }

  int operator()() 
  {
        using namespace std;
        if(arg.version) {
            return CmdVersion();
        }
        if(arg.help || arg.input.empty()) {
            return CmdHelp();
        }
        return task_(model, arg);
  }

protected:
  virtual int CmdHelp() const {
    string usage_name(arg.run_name);
    if(usage_name.substr(0, 4) == "dng-") {
      usage_name[3] = ' ';
    }
    cerr << "Usage:" << std::endl
	 << "Autosomes:" << std::endl
	 << "\tdng dnm auto --bcf bcf_f --ped ped_f [OR] dng dnm auto --vcf vcf_f --ped ped_f\n"
	 << "X chromosome in male offspring:\n"
	 << "\tdng dnm XS --bcf bcf_f --ped ped_f [OR] dng dnm XS --vcf vcf_f --ped ped_f\n"
	 << "X chromosome in female offspring:\n"
	 << "\tdng dnm XD --bcf bcf_f --ped ped_f [OR] dng dnm XD --vcf vcf_f --ped ped_f\n"
	 << endl;
    cerr << ext_desc_ << endl;
    return EXIT_SUCCESS;    
  }

};

int main(int argc, char *argv[]) {
  try {
    return DNGApp(argc, argv)();
    //return CallApp(argc, argv)();
  } catch(std::exception &e) {
    std::cerr << e.what() << std::endl;
  }

  return EXIT_FAILURE;
}
