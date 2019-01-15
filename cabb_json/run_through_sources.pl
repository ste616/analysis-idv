#!/usr/bin/perl

use strict;

# This script is a wrapper script to call the source investigation script for
# a whole bunch of sources.

my $rootdir = $ARGV[0];
my $cheatdir = "cheatplots";
print $rootdir."\n";

# Get the list of sources we can find below that directory.

my @allsources = `find -L $rootdir -name '*.json' -exec basename {} .json \\; | sort | uniq`;
for (my $i = 0; $i <= $#allsources; $i++) {
    chomp($allsources[$i]);
    print "Source ".($i + 1).": ".$allsources[$i]."\n";

    # And run source_investigation.
    my $cmd = "/n/ste616/usr/pve/bin/python source_investigation.py $allsources[$i] $rootdir $cheatdir";
    print $cmd."\n";
    if (system $cmd) {
	die "\n\nSomething has gone wrong!\n";
    }

    # Link the output directory to the web pages.
    my $link_target = "/n/ste616/usr/src/analysis-idv/cabb_json/".$cheatdir."/".$allsources[$i];
    my $link_endpoint = "/n/ste616/www/c2914_ihv_".$cheatdir."/".$allsources[$i];
    if (!-l $link_endpoint) {
	my $link_cmd = "ln -s ".$link_target." ".$link_endpoint;
	print $link_cmd."\n";
	system $link_cmd;
    }
}
