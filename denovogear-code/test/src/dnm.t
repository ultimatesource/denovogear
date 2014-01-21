#! /usr/bin/perl

use strict;
use warnings;
use Test::More;
use File::Compare;

my $op_dir = "./test/op";
ok(-e $op_dir, "test output directory exists");

my $exe = "./build/src/denovogear";
ok(-e $exe, "denovogear executable exists");

my $expected_auto_op = "./test/op/dnm_auto_op";
ok(-e $expected_auto_op, "expected auto op exists");

my $auto_op = "dnm_auto_test_current_op";
my $auto_cmd = "./build/src/denovogear dnm auto --bcf ./examples/sample_CEU.bcf --ped ./examples/sample_CEU.ped --rd_cutoff 0 --pp_cutoff 0  > $auto_op 2>&1";
my $rv_auto = system($auto_cmd);
is($rv_auto, 0, "auto return value is zero");
is(compare($auto_op, $expected_auto_op), 0, "auto expected output matches"); 
system("rm -rf $auto_op");

my $expected_xs_op = "./test/op/dnm_XS_op";
ok(-e $expected_xs_op, "XS expected output exists");

my $XS_op = "dnm_XS_test_current_op";
my $XS_cmd = "./build/src/denovogear dnm XS --bcf ./examples/sample_CEU.bcf --ped ./examples/sample_CEU.ped --rd_cutoff 0 --pp_cutoff 0  > $XS_op 2>&1";
my $rv_xs = system($XS_cmd);
is($rv_xs, 0, "XS return value is zero");
is(compare($XS_op, $expected_xs_op), 0, "XS expected output matches"); 
system("rm -rf $XS_op");

my $expected_xd_op = "./test/op/dnm_XD_op";
ok(-e $expected_xd_op, "XD expected output exists");

my $XD_op = "dnm_XD_test_current_op";
my $XD_cmd = "./build/src/denovogear dnm XD --bcf ./examples/sample_CEU.bcf --ped ./examples/sample_CEU.ped --rd_cutoff 0 --pp_cutoff 0  > $XD_op 2>&1";
my $rv_xd = system($XD_cmd);
is($rv_xd, 0, "XD return value is zero");
is(compare($XD_op, $expected_xd_op), 0, "XD expected output matches"); 
system("rm -rf $XD_op");

done_testing();
