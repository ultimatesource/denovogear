#!/usr/bin/perl -w

use strict;
use v5.10;

my @X = ('A','C','G','T', 'N');

# http://perldoc.perl.org/perlfaq4.html#How-do-I-permute-N-elements-of-a-list%3f
# http://stackoverflow.com/a/28712605

sub next_permutation {
	my @idx = @_;
    my $p = $#idx;
    --$p while $idx[$p-1] > $idx[$p];
    my $q = $p or return;
    push @idx, reverse splice @idx, $p;
    ++$q while $idx[$p-1] > $idx[$q];
    @idx[$p-1,$q]=@idx[$q,$p-1];
	return @idx;
}

my @gt_list = (0, 1, 2, 3, 1, 4, 5, 6, 2, 5, 7, 8, 3, 6, 8, 9);
sub gt10 {
	my @s = @_;
	my @r = ();
	for(my $i = 0; $i < @s; ++$i) {
		for(my $j = 0; $j <= $i; ++$j) {
			push(@r, $gt_list[4*$s[$j]+$s[$i]]);
		}
	}
	return @r;
}

my @types = ();
my @types4 = ();
my @genotypes10 = ();
for(my $k = 1;$k<=4;++$k) {
	my @n = (0,1,2,3);
	do {
		my @front = @n[0..($k-1)];
		@n = (@front, reverse(@n[$k..$#n]));
		push(@types,[@front]);
		push(@types4,[@n]);
		push(@genotypes10,[gt10(@n)]);
	} while(@n = next_permutation(@n))
}

my @ntypes = map { [4,@{$_}] } @types;
push(@types,@ntypes);

my @strings = map { join("", @X[@{$_}]) } @types;

my $command = shift || 'tsv';

if($command eq 'tsv' ) {
	for(my $u=0;$u<@types;++$u) {
		my @s = @{$types[$u]};
		shift(@s) if($s[0] == 4);
		my $list = join(",", @s);
		my $list4 = join(",", @{$types4[$u % 64]});
		my $num = @s;
		my $lc = lc($strings[$u]);
		say("$u\t$strings[$u]\t$lc\t$num\t$list4\t$list");
	}
} elsif($command eq 'cxx') {
	for(my $u=0;$u<@types;++$u) {
		my @s = @{$types[$u]};
		my $ref = $s[0]; 
		shift(@s) if($s[0] == 4);
		my $list = join(",", @s);
		my $list4 = join(",", @{$types4[$u % 64]});
		my $gt10 = join(",", @{$genotypes10[$u % 64]});
		my $num = @s;
		my $uc = $strings[$u];
		my $lc = lc($uc);
		my $comma = ($u < $#types) ? ',' : '';
		my $id = sprintf("% -4s","$u,");
		$uc = sprintf("% -8s","\"$uc\",");
		$lc = sprintf("% -8s","\"$lc\",");
		say("    {$id $num, $uc $lc $ref, {$list4}}$comma");
	}
	say("");
	for(my $u=0;$u<@types;++$u) {
		my $comma = ($u < $#types) ? ',' : '';

		my @s = @{$types[$u]};
		shift(@s) if($s[0] == 4);
		my $num = @s;
		$num = $num*($num+1)/2;
		$num = sprintf("% -3s","$num,");
		my $id = sprintf("% -4s","$u,");

		my $gt10 = join(",", @{$genotypes10[$u % 64]});
		say("    {$id $num {$gt10}}$comma");
	}	
} else {
	say("Unknown command output format '$command'");
}