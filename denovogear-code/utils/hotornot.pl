#!/usr/bin/perl
#Hot Or Not PL
#script to take genome features froms one file and 
#check for overlaps with features in a second file
#output regurgitates original file with new column 
#indicating annotation overlap
#"Overlap" can be qualitative or quantitative
#Update History: Sept 2, 2008 added reciprocal overlap
#Sept 5, 2008 : small bugfix concerning calculation of quantitative overlap
#Sept 18, 2008: cleaned up warnings
#July 14 2009: add functions read_snps and snpsort to read and sort SNPs specifically

use warnings;


$annotation=$ARGV[0];
$cnvfile=$ARGV[1];
$outfile=$ARGV[2];
$QQ=$ARGV[3]; #4 counts #3 SNPs, #2 reciprocal overlap #1 Quantitative overlap, 0 Qualitative overlap 

#should revamp this to provide both counts (qual) and quantitative overlap??

$RECIPROCAL_OVERLAP=$ARGV[4];
$PRINT_PAIR=$ARGV[5];

if ($#ARGV<3){ usage(); exit;}
if ($ARGV[3]==2 && $#ARGV==3){ die "Please specify Reciprocal overlap for mode 2\n";}

if ($QQ==2){         $my_analysis =   \&compare_cnvs;}
elsif ($QQ==3){      $my_analysis =   \&compare_snps;}
else {               $my_analysis =   \&compare_annotation;}

%first=();

$first{1}=0;
$mychr=1;
$win_size=0;
$match=0;
$curr_chr=0;

print "USING QQ $QQ\n";


if ($QQ!=3){
    @get_back=read_cnvs($cnvfile); }
else {    @get_back=read_snps($cnvfile); }

if ($#get_back==1){
    @header=@{$get_back[1]};    
    push @header,$annotation;
}

@cnv=@{$get_back[0]};



@get_back=read_data($annotation);
@ann=@{$get_back[1]};
%first=%{$get_back[0]};

open(OUT1,">$outfile") or die "Can't open outfile $outfile[0]\n";

#print header if need be
if (defined @header){
foreach $el (@header){
    if ($el=~/\S+\/(\S+)/){$el=$1;}
print OUT1 $el."\t";}
print OUT1 "\n";} 

print "screening $#cnv cnvs\n";
for ($i=0;$i<=$#cnv;$i++){

$nmatch=0;
$match=0;    
$density=0;
@lmatch=();
@rmatch=();
#case that chromosome does not exist in annotation but exists in cnv file
if (!defined $first{$cnv[$i][0]}){

for ($z=0;$z<=$#{$cnv[0]};$z++){	 

print OUT1 "$cnv[$i][$z] ";

}

print OUT1 " $match\n";

next;

    }

$temp=$first{$cnv[$i][0]};

if ($curr_chr!=$cnv[$i][0]){print "Chr $cnv[$i][0]\n";$curr_chr=$cnv[$i][0];}


while ($temp < $#ann && $ann[$temp][2]<$cnv[$i][1] && $ann[$temp][0]==$cnv[$i][0]){$temp++;}

#hack to deal with SNPS until better implementation
if ($QQ==3){$event_end=$cnv[$i][1];}
else{$event_end=$cnv[$i][2];}


while($temp<=$#ann && $ann[$temp][1]<$event_end && $ann[$temp][0]==$cnv[$i][0] ){

#   $match=compare_annotation(\@{$ann[$temp]},\@{$cnv[$i]},$win_size);
   $match=$my_analysis->(\@{$ann[$temp]},\@{$cnv[$i]},$win_size);

     if ($match!=0){
	 $nmatch++;

    if ($QQ==1){ #calculate density for quantitative overlap analysis
	     $density+= calc_overlap($ann[$temp][1],$ann[$temp][2],$cnv[$i][1],$cnv[$i][2]);}

push(@lmatch,$ann[$temp][1]);     
push(@rmatch,$ann[$temp][2]);
     
     }

   if ($temp==$#ann){last;}
   $temp++;

}


#print CNVE info
if ($nmatch!=0){

for ($z=0;$z<=$#{$cnv[0]};$z++){	 
print OUT1 "$cnv[$i][$z] ";

}


#Print overlap summary
  if ($QQ==1){
    
    $density/=($cnv[$i][2]-$cnv[$i][1])     ;
    $rdensity = sprintf("%.3f",$density);
    print OUT1 "$rdensity\n";

  } 

elsif ($QQ==4){ print OUT1 "$nmatch \n";}
#elsif ($QQ==2){ print OUT1 "$nmatch\n";}
elsif ($QQ==2){ print OUT1 "$nmatch\n";}

else {print OUT1 "1 \n";}
 
}

elsif ($nmatch==0){

for ($z=0;$z<=$#{$cnv[0]};$z++){	 
print OUT1 "$cnv[$i][$z] ";

}
if ($QQ!=2){
    print OUT1 "$nmatch\n";}
else {print OUT1 "$nmatch\n";}
#else {print OUT1 "0 0 $nmatch\n";}
}

}

close(OUT1);

##################################
## COMPARE SNPS
#################################

sub compare_snps{
    my $matched=0;
    my @a=@{$_[0]};
    my @b=@{$_[1]};
    my $win=$_[2];

    if ($a[1]<= $b[1] && $b[1]<=$a[2]){$matched=1;}

    return $matched;
}


sub compare_annotation{
    my $matched=0;
    my @a=@{$_[0]};
    my @b=@{$_[1]};
    my $win=$_[2];


#contained within 
if (($b[1]-$win)<= $a[1] && $a[2] <= ($b[2]+$win)) {
 #     print "type III ";
    $matched=1;$overlap=$a[2]-$a[1];
#print "$a[1] $a[2] $b[1] $b[2]\n";
}

#contained within 
elsif ($a[1] <= ($b[1]-$win)  &&  ($b[2]+$win) <= $a[2]) {
 #      print "type IV ";
    $matched=1; $overlap=$b[2]-$b[1];
#print "$a[1] $a[2] $b[1] $b[2]\n";
}

#overlap left edge 
elsif ($b[1]< ($a[2]-$win) && ($b[2]-$win) >$a[2]) {
 #      print "type I ";
    $matched=1; $overlap=$a[2]-$b[1];
#print "$a[1] $a[2] $b[1] $b[2]\n";
}

#overlap right edge
elsif ($a[1]> ($b[1]+$win) && ($b[2]+$win) >$a[1]) {
 #    print "type II ";
    $matched=1;$overlap=$b[2]-$a[1];
#print "$a[1] $a[2] $b[1] $b[2]\n";
}


 return $matched;
}



sub calc_overlap{
###Arg 0 = left edge of annot, Arg 1 = right edge of annotation, Arg 2= left edge of bin, Arg 3= right edge of bin
my  @data=@_;
my $overlap=0;

$overlap=$data[3]-$data[2]; # start with total bin size

if ($data[1]<$data[3]){$overlap-=($data[3]-($data[1]+1));} # if right
							   # edge of
							   # annotation
							   # inside
							   # bin

if ($data[0]>$data[2]){$overlap-=($data[0]-$data[2]);} # if left edge
						       # of annotation
						       # inside bin

#print "co $overlap $data[0] $data[1] $data[2] $data[3] $data[4]******\n";

return $overlap;  

}



###########################
#     COMPARE CNVS     ###
###########################


sub compare_cnvs{
#return 0 if overlap is not > thresh
#otherwise, return overlap

    my $matched=0;
    my @a=@{$_[0]};
    my @b=@{$_[1]};
    my $win=$_[2];
    my $overlap=0;
    my $cnv1_length=0;
    my $cnv2_length=0;
#print "$a[1] $a[2] $b[1] $b[2] win $win\n";

#contained within 
if (($b[1]-$win)<= $a[1] && $a[2] <= ($b[2]+$win)) {
 #     print "type III ";
    $matched=1;$overlap=$a[2]-$a[1];
#print "$a[1] $a[2] $b[1] $b[2]\n";
}

#contained within 
elsif ($a[1] <= ($b[1]-$win)  &&  ($b[2]+$win) <= $a[2]) {
 #      print "type IV ";
    $matched=1; $overlap=$b[2]-$b[1];
#print "$a[1] $a[2] $b[1] $b[2]\n";
}

#overlap left edge 
elsif ($b[1]< ($a[2]-$win) && ($b[2]-$win) >$a[2]) {
 #      print "type I ";
    $matched=1; $overlap=$a[2]-$b[1];
#print "$a[1] $a[2] $b[1] $b[2]\n";
}

#overlap right edge
elsif ($a[1]> ($b[1]+$win) && ($b[2]+$win) >$a[1]) {
 #    print "type II ";
    $matched=1;$overlap=$b[2]-$a[1];
#print "$a[1] $a[2] $b[1] $b[2]\n";
}



 if ($matched==1){

     $cnv1_length=$a[2]-$a[1];
     $cnv2_length=$b[2]-$b[1];

 if ($cnv1_length==0 || $cnv2_length==0){print "division by 0 in function compare CNVS\n $a[1] $a[2]\n$b[1] $b[2]\n";}

 $cnv1_overlap=$overlap/$cnv1_length;
 $cnv2_overlap=$overlap/$cnv2_length;


if ($overlap<0){die "length of overlap less than 0: $overlap \n $cnv1_length $cnv2_length\n $b[1] $b[2] $a[1] $a[2]\n";}

# print "Overlaps $overlap $cnv1_overlap $cnv2_overlap $cnv1_length $cnv2_length\n";
#If greatest overlap for both event is less than threshold, kill it
if ($cnv1_overlap<$RECIPROCAL_OVERLAP || $cnv2_overlap<$RECIPROCAL_OVERLAP){$matched=0;}

 }

    return $matched;
}



sub read_genes(){
    my $genefile=$_[0];
    my @genes=();
    my $row=0;
open(IN,"< $genefile") or die "Can't open gene file $genefile \n";
while(<IN>) {

 chomp($line=$_);
    if ($line=~/#/){next;}  
    if ($line=~/random/){next;}
    @data=split('\t',$line);
# if ($data[1]!~/chr1/){last;}
 if ($data[1]=~/chr[XY]/){next;}
 if ($data[1]=~/hla/){next;}



#if (grep(/$data[10]/,@protein_id)){next;}
 #push(@protein_id,$data[10]);
# if ($row>0 && $data[10] eq $genes[($row-1)][10]){next;}

if ($row>0 && $data[3] eq $genes[($row-1)][1]){next;} #skip genes with identical trxn start site

   $genes[$row]=[split('\t',$line)];

	if ($genes[$row][1]=~/chr(\d+)/){$genes[$row][0]=$1;   }
 $genes[$row][1]=$genes[$row][3];
 $genes[$row][2]=$genes[$row][4];
    
  if ($genes[$row][0] != $mychr){$mychr=$genes[$row][0]; $first{$mychr}=$row;}   
      $row++;
    
}

print "Read $row entries from gene file\n";
close(IN);

return @genes;}


sub read_data{
    my $file= $_[0];
    my %first=();
    my $chr=1;
    my @arr;
    $first{1}=0;
    my $scol=1;
    my $ecol=2;
    my $row=0; 
   if ($file=~/gff/){$scol=3;$ecol=4;}

open(IN,"< $file") or die "cannot open input file $file\n";
while(<IN>) {
chomp($line=$_); 
if ($line=~/#/){next;}
   if ($line=~/track/){next;}
  if ($line=~/chrY/){next ;}    
   if ($line=~/chrM/){next ;}    
   if ($line=~/hap/){next ;}
   if ($line=~/random/){next ;}        
#print "Read $line\n";
   $line =~ s/\t/ /g; 
   $arr[$row]=[split(' ',$line)];

if ($arr[$row][0]=~/X/){$arr[$row][0]=23;}
if ($arr[$row][0]=~/chr(\d+)/){$arr[$row][0]=$1;} 

  $arr[$row][1]=$arr[$row][$scol]; #put start/end into 2,3 column
  $arr[$row][2]=$arr[$row][$ecol];

#if ($row>0 && $arr[$row][1]==$arr[($row-1)][1] && $arr[$row][2]==$arr[($row-1)][2]){next;} #remove entries with identical start/stop position;

   $row++;

}

close(IN);

#SORT 
# my @as= sort multisort @arr;    
    if ($#arr<1000000){
     @as= sort { for my $ix (0 .. 2){
    my $cmp = $a->[$ix] <=> $b->[$ix];
    return $cmp if $cmp;
    } 
		return 0;
 } @arr;

    } else {@as=@arr;}

    for ($row=0;$row<=$#as;$row++){
    if ($as[$row][0]!=$chr){$first{$as[$row][0]}=$row; $chr=$as[$row][0];}}

print "read $row lines for $file first obs $as[0][0] $as[0][1] $as[0][2]\n";
print  "read $row lines for $file last obs $as[$row-1][0] $as[$row-1][1] $as[$row-1][2]\n";
return (\%first,\@as);
}


sub read_cnvs{
my    $cnvfile=$_[0];
my @cnv=();
my $row=-1;
my $mychr=1;
open(IN,"< $cnvfile") or die "Can't open CNV file $cnvfile\n";
while(<IN>) {
 chomp($line=$_);
   $line =~ s/\t/ /g; 
 if ($line=~/!/){next;}
   if ($line=~/#/){@header=split(' ',$line); next;}
   if ($line=~/chrY/){next ;}    
   if ($line=~/chrM/){next ;}    
   if ($line=~/random/){next ;}        
#   if ($line=~/hap/){next ;}        
$row++;
    $cnv[$row]=[split(' ',$line)];
  if ($cnv[$row][0]=~/chrX/){$cnv[$row][0]=23;}
  if ($cnv[$row][0] eq "X"){$cnv[$row][0]=23;}
  if ($cnv[$row][0]=~/chr(\d+)/){$cnv[$row][0]=$1;   }
  if ($cnv[$row][0] != $mychr){$mychr=$cnv[$row][0]; $first{$mychr}=$row;}  

   }

print "read $row lines for $cnvfile first obs $cnv[0][0] $cnv[0][1]\n";
print  "read $row lines for $cnvfile last obs $cnv[$row-1][0] $cnv[$row-1][1]\n";

 my @arr= sort multisort @cnv;    
   @cnv=@arr;
close(IN);
if (defined @header){
    return(\@cnv,\@header);}
else{ return \@cnv;}
}

sub read_snps{
my    $cnvfile=$_[0];
my @cnv=();
my $row=-1;
my $mychr=1;
open(IN,"< $cnvfile") or die "Can't open CNV file $cnvfile\n";
while(<IN>) {
 chomp($line=$_);
   $line =~ s/\t/ /g; 
 if ($line=~/!/){next;}
   if ($line=~/#/){@header=split(' ',$line); next;}
   if ($line=~/chrY/){next ;}    
   if ($line=~/chrM/){next ;}    
   if ($line=~/random/){next ;}        
#   if ($line=~/hap/){next ;}        
$row++;
    $cnv[$row]=[split(' ',$line)];
  if ($cnv[$row][0]=~/chrX/){$cnv[$row][0]=23;}
  if ($cnv[$row][0] eq "X"){$cnv[$row][0]=23;}
  if ($cnv[$row][0]=~/chr(\d+)/){$cnv[$row][0]=$1;   }
  if ($cnv[$row][0] != $mychr){$mychr=$cnv[$row][0]; $first{$mychr}=$row;}  

   }

print "read $row lines for $snpfile first obs $cnv[0][0] $cnv[0][1]\n";
print  "read $row lines for $snpfile last obs $cnv[$row-1][0] $cnv[$row-1][1]\n";

 my @arr= sort snpsort @cnv;    
   @cnv=@arr;
close(IN);
if (defined @header){
    return(\@cnv,\@header);}
else{ return \@cnv;}
}


sub multisort {
    
# first split $a and $b into fields
($field1a, $field2a, $field3a)=@{$a};
($field1b, $field2b, $field3b)=@{$b};

# now do the actual sort:
$field1a <=> $field1b
or
$field2a <=> $field2b
or
$field3a <=> $field3b;
}
sub snpsort {
    
# first split $a and $b into fields
($field1a, $field2a, $field3a)=@{$a};
($field1b, $field2b, $field3b)=@{$b};

# now do the actual sort:
$field1a <=> $field1b
or
    $field2a <=> $field2b;


}


sub usage{
print "Usage:
annotation=ARGV[0];
cnvfile=ARGV[1];
outfile=ARGV[2];
QQ=ARGV[3]; #4 use CNV but print counts of hit to each CNV #3 use SNP instead of CNV, 2 Reciprocal overlap, 1 Quantitative overlap, 0 Qualitative overlap
Reciprocal Overlap Threshold = ARGV[4]; (Optional) if QQ==2
\n";}
