#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include "bgzf.h"
#include <math.h>

#define WANT_STREAM       // include iostream and iomanipulators
#include "newmatap.h"
#include "newmatio.h"

using namespace std;

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

#ifdef use_namespace
using namespace RBD_LIBRARIES;
#endif

#include "glf.h"

#ifndef VERSION
#define VERSION "dummy"
#endif

#define MRATE 5e-7
#define MIN_READ_DEPTH 0
#define MIN_MAPQ 40 


int autLength[24] = { 0, 247249719, 242951149, 199501827, 191263066, 180837866, 
		      170896992, 158821424, 146274818, 140273252, 135374737, 
		      134452384, 132289542, 114127979, 106360583, 100338915, 
		      88822253, 78654741, 76117153, 63806651, 62435958, 
		      46944323, 49591432, 154913754} ;
int X_length = 154913754 ;


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
} lookup_t;


static void trio_like(glf3_t *g1_c,glf3_t *g1_m,glf3_t *g1_d, lookup_t & lookup, vector<vector<string > >  & tgt,int coor,int flag, char* ref_name)
{




   Real a[10];   
   Real maxlike_null,maxlike_denovo,pp_null,pp_denovo,denom,numer;   

  Matrix M(1,10);
  Matrix C(10,1);
  Matrix D(10,1);
  Matrix P(10,10);
  Matrix F(100,10);
  Matrix L(100,10);
  Matrix T(100,10);
  Matrix DN(100,10);
  Matrix PP(100,10);
 int i,j,k,l,z;
 double resc_m,resc_d,resc_c,mymin,mind;

  //Rescale Likelihoods
  //find minimum value among samples
 
 // if (g1_m->min_lk > g1_d->min_lk){ mymin= g1_c->min_lk> g1_m->min_lk ? g1_c->min_lk : g1_m->min_lk;}
 //else {mymin= g1_c->min_lk> g1_d->min_lk ? g1_c->min_lk : g1_d->min_lk;}

 // mind= pow(10,-mymin/10.);
 //resc_m=pow(10,- (double) g1_m->min_lk/10.)/mind;
 //resc_d=pow(10,- (double) g1_d->min_lk/10.)/mind;
 //resc_c=pow(10,- (double) g1_c->min_lk/10.)/mind;

 // cout<<"MIN LK "<<g1_m->min_lk<<" "<<resc_m<<" "<<resc_d<<" "<<resc_c<<" "<<mymin<<endl;

  //Load vectors

 
        for (j = 0; j != 10; ++j) a[j]=pow(10,-g1_m->lk[j]/10.);
 M<<a;
     for (j = 0; j != 10; ++j) a[j]=pow(10,-g1_d->lk[j]/10.);
D<<a;
for (j = 0; j != 10; ++j) a[j]=pow(10,-g1_c->lk[j]/10.);
 C<<a;
        

	P=KP(M,D);

	F=KP(P,C);

	

	T=SP(F,lookup.tp); //combine with transmission probs 


        switch("XACMGRSVTWYHKDBN"[g1_m->ref_base]){
         
	case 'A':
	  L=SP(T,lookup.aref); break;

        case 'C':
	  L=SP(T,lookup.cref); break;

        case 'G':
	  L=SP(T,lookup.gref); break;

        case 'T':
	  L=SP(T,lookup.tref); break;

	default: L=T; break;
	}

	DN=SP(L,lookup.mrate);

        PP=SP(DN,lookup.norm);   //zeroes out configurations with mendelian error
        maxlike_null = PP.maximum2(i,j);       

	//Find max likelihood of de novo trio configuration
	PP=SP(DN,lookup.denovo);   //zeroes out configurations with mendelian inheritance
        maxlike_denovo=PP.maximum2(k,l); 

  

	  //make proper posterior probs

	  denom=DN.sum();
          numer=PP.sum();
          pp_denovo=maxlike_denovo/denom;
          pp_null=1-pp_denovo;

//cout << setw(10) << setprecision(5) << scientific << PP(k,l) << endl;
// cout<<" *****************"<<endl;
//cout << setw(10) << setprecision(5) << scientific << L(k,l) << endl;	  
//cout<<" *****************"<<endl;
//cout << setw(10) << setprecision(5) << scientific << lookup.mrate(k,l) << endl;	  
//cout<<" *****************"<<endl;
//cout << setw(10) << setprecision(5) << scientific << lookup.denovo(k,l) << endl;	  


//cout << setw(10) << setprecision(5) << scientific << PP(i,j) << endl;
// cout<<" *****************"<<endl;
//cout << setw(10) << setprecision(5) << scientific << L(i,j) << endl;	  
//cout<<" *****************"<<endl;
//cout << setw(10) << setprecision(5) << scientific << lookup.mrate(i,j) << endl;	  
//cout<<" *****************"<<endl;
//cout << setw(10) << setprecision(5) << scientific << lookup.denovo(i,j) << endl;	  

			     	     	  
  //                    for (z = 0; z != 10; ++z) printf("\t%d", g1_c->lk[z]);
//		      cout<<" blah "<<g1_c->min_lk <<"\n";
//		      for (z = 0; z != 10; ++z) printf("\t%d", g1_m->lk[z]);
//		          cout<<" blah "<<g1_m->min_lk <<"\n";
//		       for (z = 0; z != 10; ++z) printf("\t%d", g1_d->lk[z]);
//		          cout<<" blah "<<g1_d->min_lk <<"\n";


	  if (pp_denovo>0.01){
	  cout<<ref_name<<" "<<coor<<" "<<"XACMGRSVTWYHKDBN"[g1_m->ref_base]<<" "<<maxlike_null<<" "<<pp_null<<" "<<tgt[i-1][j-1]<<" "<<lookup.snpcode(i,j)<<" "<<lookup.code(i,j)<<" "<<maxlike_denovo<<" "<<pp_denovo<<" "<<tgt[k-1][l-1]<<" "<<lookup.code(k,l)<<" "<<flag;
	  //printf(" %d\t%d\t%c\t",chr,coor,"XACMGRSVTWYHKDBN"[g1_m.ref_base]);
	  // printf("%lf\t%s\t",maxlike_null,tgt[i-1][j-1]);
	  //printf("%d\t%d\t%g\t",lookup.snpcode(i,j),lookup.code(i,j),maxlike_denovo);
 //printf("%s\t%d\t",tgt[k-1][l-1],lookup.code(k,l));
      printf(" depth %d\t%d\t%d\t",g1_c->depth,g1_d->depth,g1_m->depth);
      printf(" Qs %d\t%d\t%d\n",g1_c->rms_mapQ,g1_d->rms_mapQ,g1_m->rms_mapQ);
	  }
     

}

static void trio_scan(const char *fn_samp1, const char *fn_samp2, const char *fn_samp3, vector<vector<string > >  & tgt, lookup_t & 
lookup)
{
	FILE *fs;
        
	char str[1024];
	glf3_t *g1_d, *g1_c, *g1_m;
	int coor=1,mcoor=1,dcoor=1,i,saw=0,len,flag=0;
        glf3_header_t *h;
	glfFile g_c,g_d,g_m;
	//	char *ref_name;
	char * ref_name;
          

  g1_c = glf3_init1();
  g1_d = glf3_init1();
  g1_m = glf3_init1();
  
if (!(g_c = bgzf_open(fn_samp1,"r")))
  {fprintf (stderr, "failed to open .glz file %s\n", fn_samp1); exit(1);}

       h = glf3_header_read(g_c);
  


       if (! (	g_m = bgzf_open(fn_samp2,"r")))
  {fprintf (stderr, "failed to open .glz file %s\n", fn_samp2); exit(1);}

        h = glf3_header_read(g_m);

	if(! ( g_d = bgzf_open(fn_samp3,"r")))
  {fprintf (stderr, "failed to open .glz file %s\n", fn_samp3); exit(1);}

        h = glf3_header_read(g_d);


	while ((ref_name = glf3_ref_read(g_c, &len)) != 0) {





		
	  cerr<<"read value of l "<<len<<" refname "<<ref_name<<endl;
	 		free(ref_name);

                ref_name = glf3_ref_read(g_d, &len);


	  cerr<<"read value of l "<<len<<" refname "<<ref_name<<endl;
	  free(ref_name);
		
		//		ref_name = (char*)calloc(l, 1);


		ref_name = glf3_ref_read(g_m, &len);
	      	cerr<<"read value of len "<<len<<" refname "<<ref_name<<endl;	


		//read through all l bases in the chromosome
	      	//for (i=0;i<len;i++){		
		//		cout<<" TEST ** "<<g1_c->rtype<<" END "<<endl;  
		  while(glf3_read1(g_c, g1_c) && g1_c->rtype !=  GLF3_RTYPE_END){



  coor = coor + g1_c->offset;

 while (g1_c->rtype ==  GLF3_RTYPE_INDEL){
   //cout<<"hit "<<endl;
glf3_read1(g_c, g1_c);
   //  cout<<"I Position "<<g1_c->offset << g1_d->offset <<g1_m->offset<<" Type "<<g1_c->rtype<<" "<<g1_d->rtype<<" "<<g1_m->rtype<<" "<< GLF3_RTYPE_INDEL<< endl;

 coor = coor + g1_c->offset;
	 
 }




		

		  while( mcoor <coor || g1_m->rtype ==  GLF3_RTYPE_INDEL){	  	 
          	  glf3_read1(g_m, g1_m);
                  mcoor = mcoor + g1_m->offset;

		  }
		  
	

		  while( dcoor <coor || g1_d->rtype ==  GLF3_RTYPE_INDEL){	  	 
		  glf3_read1(g_d, g1_d);
                  dcoor = dcoor + g1_d->offset;}
     
	
		  if (dcoor != coor || mcoor !=coor){continue;} //

		  //		  		  cout<<"Position "<<coor<<" Type "<<g1_c->rtype<<" "<<g1_d->rtype<<" "<<g1_m->rtype<<" "<< GLF3_RTYPE_INDEL<< endl;

		 
        if (g1_c->rms_mapQ < MIN_MAPQ || g1_d->rms_mapQ < MIN_MAPQ || g1_m->rms_mapQ < MIN_MAPQ){flag=1;}
	if (g1_c->depth < MIN_READ_DEPTH || g1_d->depth < MIN_READ_DEPTH || g1_m->depth < MIN_READ_DEPTH){flag=1;}


	//	if (coor % 100000 == 0){cerr<<"Made it to "<<coor<<endl;}

	
		trio_like(g1_c,g1_m,g1_d,lookup,tgt,coor,flag,ref_name);


		flag=0;
                saw++;

		  }
		    free(ref_name);
		cerr<<" parsed "<<saw<<" good sites out of "<<i<<endl;

}



	//clean everything up

        glf3_header_destroy(h);
//        glf3_header_destroy(h2);
//      glf3_header_destroy(h3);

	glf3_destroy1(g1_c);
	glf3_destroy1(g1_d);
	glf3_destroy1(g1_m);

	bgzf_close(g_c);
	bgzf_close(g_d);
	bgzf_close(g_m);

}






/* commands */


static int read_lookup( vector<vector<string > >  & tgt, lookup_t & lookup)
{
  string line;
  int tmp[1000];
  float tmp2[1000];
  float tmp3[1000];
  float tmp4[1000];
  float tmp5[1000];
  float tmp6[1000];
  float tmp7[1000];
  float tmp8[1000];
  float tmp9[1000];
  float tmp10[1000];

  //  vector<string> tmp4; 
  float blah;
  int i=0,j=0,k=0;
  string blah2;

  ifstream lookup_file("lookup.txt");
  if (!lookup_file){cerr <<"cannot open lookup table"<<endl; exit(1);}
  while (getline(lookup_file,line)){

    istringstream iss(line);
    iss >> blah;
    tmp[k]= (int) blah;
    iss >> blah;
    tmp2[k]=blah;
    iss >> blah;
    tmp3[k]=blah;
    iss >> blah2;
    tgt[i][j++]=blah2;

    if (j==10){i++;j=0;}
      
    iss >> blah;
    tmp4[k]=blah;
    iss >> blah;
    tmp5[k]=blah;
    iss >> blah;
    tmp6[k]=blah;
    iss >> blah;
    tmp7[k]=blah;
    iss >> blah;
    tmp8[k]=blah;
    iss >> blah;
    tmp9[k]=blah;
    iss >> blah;
    tmp10[k]=blah;
    k++;

     }
      lookup_file.close();
      lookup.code << tmp2;
      lookup.tp << tmp3;
      lookup.snpcode << tmp;
      lookup.mrate << tmp4;
      lookup.denovo << tmp5;
      lookup.norm << tmp6;
      lookup.aref << tmp7;
      lookup.cref << tmp8;
      lookup.gref << tmp9;
      lookup.tref << tmp10;

      cerr<<" Read "<<k <<" elements from lookup"<<endl;
      cerr <<" First mrate"<< lookup.mrate(1,1) <<" Last "<< lookup.mrate(100,10)<<endl;
      cerr <<" First "<< lookup.code(1,1) <<" Last "<< lookup.code(100,10)<<endl;
      cerr <<" First "<< tgt[0][0] <<" Last "<< tgt[99][9] <<endl;
      cerr <<" First "<< lookup.tref(1,1) <<" Last "<< lookup.tref(100,10) <<endl;
      return 0; 
}


static int usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program: glf_utils (GLF utilities)\n");
	fprintf(stderr, "Version: %s\n", VERSION);
	fprintf(stderr, "Contact: Heng Li <lh3@sanger.ac.uk>\n\n");
	fprintf(stderr, "Usage:   glf_utils <command> [options]\n\n");
	fprintf(stderr, "Command: info       get reference information from .glz\n");
	fprintf(stderr, "         glf2glt    convert .glz to .glt, the text format\n");
	fprintf(stderr, "\n");
	return 1;
}

int main(int argc, char *argv[])
{

  

  Matrix code(100,10);

  Matrix tp(100,10);
 
  vector<string> tmp;
  vector<vector<string > > tgt;

  lookup_t lookup;
 
lookup.aref.resize(100,10);
lookup.cref.resize(100,10);
lookup.gref.resize(100,10);
lookup.tref.resize(100,10);  
lookup.snpcode.resize(100,10);  
lookup.code.resize(100,10);  
lookup.tp.resize(100,10);  
lookup.mrate.resize(100,10);
lookup.denovo.resize(100,10);
lookup.norm.resize(100,10);

  for (int j=0;j<10;j++){tmp.push_back("NA");}
  for (int i=0;i<100;i++){tgt.push_back(tmp);}



  //  printf("%s\n",argv[1]);
  read_lookup(tgt,lookup);
  trio_scan(argv[1],argv[2],argv[3],tgt,lookup);
 
	return 0;
}
