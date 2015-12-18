#!/usr/bin/perl -w

use strict;
use v5.10;

my @X = ('a','c','g','t', 'n');

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

my $number = 0;
for(my $z = 0;$z<2;++$z) {
	for(my $k = 1;$k<=4;++$k) {
		my @n = (0,1,2,3);
		unshift(@n,4) if($z);
		do {
			my @front = @n[0..($k+$z-1)];
			my $str = "$number\t" . join("",@X[@front]);
			say($str);
			++$number;
			@n = (@front, reverse(@n[($k+$z)..$#n]));
		} while(@n = next_permutation(@n))
	}
}
