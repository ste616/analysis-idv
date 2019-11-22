#!/usr/bin/perl

# After normal initial reduction, this script will figure out the
# time ranges to do time-dependent recalibration, and do the calibration.

use Astro::Time;
use DateTime;
use Data::Dumper;
use JSON;
use POSIX;
use Math::Round qw(round);
use strict;

# Globals.
my %station_positions;
my @configuration_strings;
# ATCA location.
my ($atca_lat, $atca_long, $atca_alt) = ( -30.3128846, 149.55013, (236.87 / 1000.0) );
my ($lat, $lon, $alt) = ( deg2rad($atca_lat), 
			  deg2rad($atca_long), ($atca_alt) );



# The argument is the name of the Miriad dataset that came straight out
# of atlod.
my $top_dataset = $ARGV[0];

if (!-d $top_dataset) {
    die "The supplied argument $top_dataset is not found.\n";
}

# The second optional argument is the name of the set to use as the calibration
# set, and defaults to 1934-638.
my $calset = "1934-638";
if ($#ARGV == 1) {
    $calset = $ARGV[1];
}

# The amount of data in each time segment, in minutes.
# The amount we would prefer.
my $segtime = 2.0;
# But if the amount is just over, we include it too.
my $extratime = 1.0;

# We're going to get the uvindex from this file and manipulate it.
my $uvindex_cmd = "uvindex vis=".$top_dataset;
my @uvindex_output = `$uvindex_cmd`;

# Parse all the information.
my @freqconfigs;
my %sources;
my $f = -1;
my $osrc = "";
for (my $i = 0; $i <= $#uvindex_output; $i++) {
    # Get rid of the newline.
    chomp($uvindex_output[$i]);
    # Split the line.
    my @els = split(/\s+/, $uvindex_output[$i]);

    # Work out the frequency configurations.
    if ($uvindex_output[$i] =~ /^Frequency Configuration (\d+)$/) {
	$f = $1 - 1;
	$freqconfigs[$f] = [];
    } elsif ($f >= 0) {
	if ($uvindex_output[$i] eq "") {
	    $f = -1;
	} else {
	    if ($els[1] ne "Channels") {
		push @{$freqconfigs[$f]}, [ $els[1], $els[2], $els[3] ];
	    }
	}
    }

    # Parse our way through the time lists.
    if (($els[0] =~ /JAN/) || ($els[0] =~ /FEB/) || ($els[0] =~ /MAR/) ||
	($els[0] =~ /APR/) || ($els[0] =~ /MAY/) || ($els[0] =~ /JUN/) ||
	($els[0] =~ /JUL/) || ($els[0] =~ /AUG/) || ($els[0] =~ /SEP/) ||
	($els[0] =~ /OCT/) || ($els[0] =~ /NOV/) || ($els[0] =~ /DEC/)) {

	if ($osrc ne "") {
	    my $enddt = &miriad2datetime($els[0]);
	    $sources{$osrc}->{'times'}->[$#{$sources{$osrc}->{'times'}}]->[1] = $els[0];
	    $sources{$osrc}->{'datetimes'}->[$#{$sources{$osrc}->{'datetimes'}}]->[1] = $enddt;
	    $osrc = "";
	}

	my $src = $els[1];
	if (($src ne "1934-638") && ($src !~ /focus/) && ($src ne "Total")) {
	    # This is a valid source.
	    if (!defined $sources{$src}) {
		# Add the necessary structure.
		$sources{$src} = { 'ra' => "", 'dec' => "", 'times' => [], 'datetimes' => []  };
	    }
	    my $startdt = &miriad2datetime($els[0]);
	    push @{$sources{$src}->{'times'}}, [ $els[0], "", $els[6] ];
	    push @{$sources{$src}->{'datetimes'}}, [ $startdt, 0, $els[6] ];
	    $osrc = $src;
	}
    }

    # Get source coordinates.
    if (($#els == 5) && (defined $sources{$els[0]})) {
	$sources{$els[0]}->{'ra'} = $els[2];
	$sources{$els[0]}->{'dec'} = $els[3];
    }

}

# Get the cycle time.
my $ct_file = "cycletime.log";
my $varplt_cmd = "varplt vis=".$top_dataset." device=/null log=".$ct_file." xaxis=time ".
    "yaxis=inttime";
if (-e $ct_file) {
    system "rm ".$ct_file;
}
system $varplt_cmd;
open(C, $ct_file) || die "Unable to determine cycle time.\n";
my %ctimes;
while(<C>) {
    chomp;
    my @els = split(/\s+/);
    if (($els[0] ne "#") && ($#els == 3)) {
	if (!defined $ctimes{$els[3]}) {
	    $ctimes{$els[3]} = 1;
	} else {
	    $ctimes{$els[3]} += 1;
	}
    }
}
close(C);
my @ckeys = keys %ctimes;
my $cycletime = 0;
if ($#ckeys == 0) {
    # Good there's only one cycle time.
    print "cycle time is ".$ckeys[0]."\n";
    $cycletime = $ckeys[0];
} else {
    # Fuck.
    print "multiple cycle times detected.\n";
}

# Now properly round the start and end datetimes.
my $rnd_cycle = round($cycletime * 5) / 10;
print "rounded half cycle time is ".$rnd_cycle."\n";
my $halfcycle = DateTime::Duration->new(
    seconds => floor($rnd_cycle),
    nanoseconds => ($rnd_cycle - floor($rnd_cycle)) * 1e9
    );
my $datadir = "data";
if (!-d $datadir) {
    system "mkdir $datadir";
}

foreach my $s (keys %sources) {
    print "examining source ".$s."\n";
    my $jfile = $datadir."/".$s.".json";
    my %jobj = (
	'sourceName' => $s,
	'rightAscension' => $sources{$s}->{'ra'},
	'declination' => $sources{$s}->{'dec'},
	'arrayConfiguration' => "",
	'timeSeries' => []
	);
    for (my $i = 0; $i <= $#{$sources{$s}->{'datetimes'}}; $i++) {
	$sources{$s}->{'datetimes'}->[$i]->[0]->subtract($halfcycle);
	$sources{$s}->{'datetimes'}->[$i]->[1]->add($halfcycle);
	print "  time range  ".$sources{$s}->{'times'}->[$i]->[0]." to ".
	    $sources{$s}->{'times'}->[$i]->[1]."\n";
	print "  is actually ".&datetime2miriad($sources{$s}->{'datetimes'}->[$i]->[0])." to ".
	    &datetime2miriad($sources{$s}->{'datetimes'}->[$i]->[1])."\n";
	my $scanlength = $sources{$s}->{'datetimes'}->[$i]->[1]->delta_ms(
	    $sources{$s}->{'datetimes'}->[$i]->[0]);
	my $scanlength_minutes = $scanlength->minutes() + $scanlength->seconds() / 60.0;
	print "    this is a delta of ".$scanlength_minutes." minutes \n";
	my $iter_starttime = $sources{$s}->{'datetimes'}->[$i]->[0]->clone();
	my $iter_endtime;
	my $finished = 0;
	if ($scanlength_minutes < ($segtime + $extratime)) {
	    $iter_endtime = $sources{$s}->{'datetimes'}->[$i]->[1]->clone();
	    $finished = 1;
	} else {
	    $iter_endtime = $iter_starttime->clone();
	    $iter_endtime->add( minutes => floor($segtime), 
				seconds => floor(($segtime - floor($segtime)) * 60.0) );
	}
	my $rc = &calibrate($s, $freqconfigs[$sources{$s}->{'datetimes'}->[$i]->[2] - 1],
			    $iter_starttime, $iter_endtime);
	if ($rc->{'code'} == 1) {
	    # We were able to do the calibration, so we can now make a measurement.
	    my $mset = &measure($s, $sources{$s},
				$freqconfigs[$sources{$s}->{'datetimes'}->[$i]->[2] - 1],
				$iter_starttime, $iter_endtime, $rc);
	    #print Dumper $mset;
	    if ($#{$mset->{'fluxDensityData'}->[0]->{'data'}} > -1) {
		# Copy the required parameters.
		my $tso = { 'fluxDensityData' => $mset->{'fluxDensityData'},
			    'flaggedFraction' => $mset->{'flaggedFraction'},
			    'defect' => $mset->{'defect'},
			    'mjdRange' => $mset->{'mjdRange'},
			    'closurePhase' => $mset->{'closurePhase'},
			    'hourAngleRange' => $mset->{'hourAngleRange'},
			    'centreFreqs' => $rc->{'frequencies'}
		};
		push @{$jobj{'timeSeries'}}, $tso;
		if ($i == 0) {
		    # Copy the array configuration.
		    $jobj{'arrayConfiguration'} = $mset->{'arrayConfiguration'};
		}
	    }
	}
	while ($finished == 0) {
	    $iter_starttime = $iter_endtime;
	    $scanlength = $sources{$s}->{'datetimes'}->[$i]->[1]->delta_ms(
		$iter_starttime);
	    $scanlength_minutes = $scanlength->minutes() + $scanlength->seconds() / 60.0;
	    if ($scanlength_minutes < ($segtime + $extratime)) {
		$iter_endtime = $sources{$s}->{'datetimes'}->[$i]->[1]->clone();
		$finished = 1;
	    } else {
		$iter_endtime = $iter_starttime->clone();
		$iter_endtime->add( minutes => floor($segtime),
				    seconds => floor(($segtime - floor($segtime)) * 60.0) );
	    }
	    $rc = &calibrate($s, $freqconfigs[$sources{$s}->{'datetimes'}->[$i]->[2] - 1],
			     $iter_starttime, $iter_endtime);
	    if ($rc->{'code'} == 1) {
		# We were able to do the calibration, so we can now make a measurement.
		my $mset = &measure($s, $sources{$s},
				    $freqconfigs[$sources{$s}->{'datetimes'}->[$i]->[2] - 1],
				    $iter_starttime, $iter_endtime, $rc);
		#print Dumper $mset;
		if ($#{$mset->{'fluxDensityData'}->[0]->{'data'}} > -1) {
		    # Copy the required parameters.
		    my $tso = { 'fluxDensityData' => $mset->{'fluxDensityData'},
				'flaggedFraction' => $mset->{'flaggedFraction'},
				'defect' => $mset->{'defect'},
				'mjdRange' => $mset->{'mjdRange'},
				'closurePhase' => $mset->{'closurePhase'},
				'hourAngleRange' => $mset->{'hourAngleRange'},
				'centreFreqs' => $rc->{'frequencies'}
		    };
		    push @{$jobj{'timeSeries'}}, $tso;
		}
	    }
	}
    }
    # Write out the JSON.
    open(J, ">".$jfile);
    print J to_json \%jobj;
    close(J);

}

#print Dumper %sources;

sub freqconfig2band {
    my $fconfig = shift;

    # What band?
    my $band = "";
    if ($fconfig->[0]->[1] < 3.8) {
	$band = "16cm";
    } elsif ($fconfig->[0]->[1] < 12) {
	$band = "4cm";
    } elsif ($fconfig->[0]->[1] < 30) {
	$band = "15mm";
    }

    return $band;
}

sub changebanddir {
    my $band = shift;

    # Change to that directory.
    chomp(my $cdir = `find . -type d -name '$band'`);
    if (-d $cdir) {
	print "changing directory to $cdir\n";
	chdir $cdir;
	return 0;
    } else {
	$cdir = "v3_redo/".$band;
	if (-d $cdir) {
	    print "changing directory to $cdir\n";
	    chdir $cdir;
	    return 0;
	}	    
	print "Unable to find appropriate band directory $cdir.\n";
	return -1;
    }

}

sub calibrate {
    my $srcname = shift;
    my $fconfig = shift;
    my $stime = shift;
    my $etime = shift;

    my $rv = { 'code' => 0, 'frequencies' => [] };
    
    print "calibrating source $srcname\n";
    my $band = &freqconfig2band($fconfig);
    
    if (&changebanddir($band) == -1) {
	$rv->{'code'} = -1;
	return $rv;
    }
    
    # Find the frequencies present for this source.
    my @datasets = glob "$srcname.*/visdata";
    my @cfreqs;
    for (my $i = 0; $i <= $#datasets; $i++) {
	print "examining dataset $datasets[$i]\n";
	if ($datasets[$i] =~ /^.*\.(.*)\/visdata$/) {
	    # Is this one of the frequencies in this configuration?
	    my $f = $1;
	    my $fcfg = -1;
	    for (my $j = 0; $j <= $#{$fconfig}; $j++) {
		print Dumper($fconfig->[$j]);
		my $c = ($fconfig->[$j]->[1] + ($fconfig->[$j]->[0] - 1) *
			 $fconfig->[$j]->[2] / 2) * 1000;
		print "c = $c  f = $f\n";
		my $diff = round(abs($c - $f));
		#print "c - f = ".(abs($c - $f))."\n";
		#print "cmp = ".(abs($fconfig->[$j]->[2] * 1000))."\n";
		if ($diff <= abs($fconfig->[$j]->[2] * 1000)) {
		    print "match found\n";
		    $fcfg = $j;
		    last;
		}
	    }
	    if ($fcfg != -1) {
		print "adding frequency $f which is cfg $fcfg\n";
		push @cfreqs, [ $f, $fcfg ];
		push @{$rv->{'frequencies'}}, $f * 1;
	    }
	}
    }

    
    # Do the calibration per frequency.
    for (my $i = 0; $i <= $#cfreqs; $i++) {
	my $dset = $srcname.".".$cfreqs[$i]->[0];
	my $cset = $calset.".".$cfreqs[$i]->[0];
	# Check that we actually have data in this time range.
	my $logout = "checkn.log";
	if (-e $logout) {
	    system "rm ".$logout;
	}
	my $smtime = &datetime2miriad($stime);
	my $emtime = &datetime2miriad($etime);
	print "dealing with time range $smtime to $emtime\n";
	my $uvplt_cmd = "uvplt vis=".$dset." axis=time,amp ".
	    "device=/null options=nopol,nocal,nopass ".
	    "stokes=xx,yy \"select=time(".$smtime.",".
	    $emtime.")\" > ".$logout;
	system $uvplt_cmd;
	my $npoints = 0;
	if (-e $logout) {
	    open(L, $logout);
	    while(<L>) {
		chomp(my $line = $_);
		if ($line =~ /^Read (.*) visibilities from all files/) {
		    $npoints = $1;
		    last;
		}
	    }
	    close(L);
	}
	if ($npoints > 0) {
	    print "found $npoints points for this time range\n";
	    # Do the calibration.
	    my $gpcal_cmd = "gpcal vis=".$dset." interval=0.1 ".
		"options=xyvary,nopol nfbin=2 refant=3";
	    $gpcal_cmd .= " \"select=time(".
		$smtime.",".$emtime.")\"";
	    print $gpcal_cmd."\n";
	    system $gpcal_cmd;
	    my $gpboot_cmd = "gpboot vis=".$dset." cal=".$cset;
	    $gpboot_cmd .= " \"select=time(".
		$smtime.",".$emtime.")\"";
	    print $gpboot_cmd."\n";
	    system $gpboot_cmd;
	    $rv->{'code'} = 1;
	} else {
	    print "unable to find any data for this time range\n";
	    $rv->{'code'} = 0;
	}
    }

    # Change our directory back to where we were.
    chdir "../..";

    return $rv;
}

sub measure {
    my $srcname = shift;
    my $coords = shift;
    my $fconfig = shift;
    my $stime = shift;
    my $etime = shift;
    my $calresult = shift;
    
    print "measuring source $srcname\n";
    my $band = &freqconfig2band($fconfig);

    if (&changebanddir($band) == -1) {
	return;
    }

    my %ro = (
	'source' => $srcname,
#	'epochName' => $ename,
	'rightAscension' => $coords->{'ra'},
	'declination' => $coords->{'dec'},
	'mjdRange' => { 
	    'low' => $stime->mjd(), 
	    'high' => $etime->mjd() },
	'arrayConfiguration' => "",
	'hourAngleRange' => { 'low' => 0, 'high' => 0 },
	'closurePhase' => [],
	'fluxDensityFits' => [],
	'fluxDensityData' => [],
	'defect' => [],
	'flaggedFraction' => []
	);

    my @stokes = ( 'i' );
    my @orders = ( 3 );
    my @options = ( ",log" );

    # Make some necessary directories.
    my $plotdir = "plots";
    if (!-d $plotdir) {
	system "mkdir $plotdir";
    }

    for (my $i = 0; $i <= $#stokes; $i++) {
	my $p = $plotdir."/".$srcname;
	my $o = $datadir."/".$srcname;
	if ($stokes[$i] ne "i") {
	    $o .= ".".$stokes[$i];
	}
	$o .= ".plotgen";

	if ($i == 0) {
	    # Get the closure phase and flagging statistic.
	    for (my $k = 0; $k <= $#{$calresult->{'frequencies'}}; $k++) {
		my $p2 = $srcname.".".$calresult->{'frequencies'}->[$k];
		my %clop = &measure_closure_phase($p2, $stime, $etime);
		push @{$ro{'closurePhase'}}, {
		    'IF' => $calresult->{'frequencies'}->[$k] * 1,
		    'average_value' => $clop{'closure_phase'}->{'average_value'} * 1.0,
		    'measured_rms' => $clop{'closure_phase'}->{'measured_rms'} * 1.0,
		    'theoretical_rms' => $clop{'closure_phase'}->{'theoretical_rms'} * 1.0
		};
		push @{$ro{'flaggedFraction'}}, &measure_flagging_statistic($p2, $stime, $etime);
	    }

	    # Calculate the hour angle ranges.
	    my $minlst = mjd2lst($ro{'mjdRange'}->{'low'}, $atca_long);
	    my $maxlst = mjd2lst($ro{'mjdRange'}->{'high'}, $atca_long);
	    my $raturns = str2turn($ro{'rightAscension'}, "H");
	    my $loha = (($minlst - $raturns) * 24.0);
	    if ($loha < -12) {
		$loha += 24.0;
	    } elsif ($loha > 12) {
		$loha -= 24.0;
	    }
	    my $hiha = (($maxlst - $raturns) * 24.0);
	    if ($hiha < -12) {
		$hiha += 24.0;
	    } elsif ($loha > 12) {
		$hiha -= 24.0;
	    }
	    if ($hiha < $loha) {
		my $tha = $loha;
		$loha = $hiha;
		$hiha = $tha;
	    }

	    my $halow = sprintf "%.4f", $loha;
	    my $hahigh = sprintf "%.4f", $hiha;
	    $ro{'hourAngleRange'}->{'low'} = $halow * 1.0;
	    $ro{'hourAngleRange'}->{'high'} = $hahigh * 1.0;

	    # Get the array configuration.
	    my $ap = $srcname.".".$calresult->{'frequencies'}->[0];
	    $ro{'arrayConfiguration'} = &determine_array($ap);
	}
	
	# Measure the flux densities.
	my $outplot = sprintf("%s/%s_%.4f.ps/cps", $plotdir, $srcname,
			      (($ro{'mjdRange'}->{'low'} +
				$ro{'mjdRange'}->{'high'}) / 2.0));
	if (-e $outplot) {
	    system "rm $outplot";
	}
	my $pcmd = "uvfmeas order=".$orders[$i]." stokes=".$stokes[$i].
	    " device=".$outplot." \"select=time(".
	    &datetime2miriad($stime).",".&datetime2miriad($etime).")\" ".
	    "vis=";
	for (my $j = 0; $j <= $#{$calresult->{'frequencies'}}; $j++) {
	    if ($j > 0) {
		$pcmd .= ",";
	    }
	    $pcmd .= $srcname.".".$calresult->{'frequencies'}->[$j];
	}
	$pcmd .= "";
	my $pvcmd = $pcmd." options=plotvec,mfflux,machine,malpha".$options[$i].
	    " log=fluxdensities.txt > tmp.log";
	system "rm tmp.log fluxdensities.txt";
	system $pvcmd;
	    
	push @{$ro{'fluxDensityFits'}}, &readlog("tmp.log");
	my $vecdens = &readdata("fluxdensities.txt");
	push @{$ro{'fluxDensityData'}}, {
	    'stokes' => 'I', 'mode' => 'vector',
	    'data' => $vecdens
	};
	    
	my $pscmd = $pcmd." options=machine".$options[$i].
	    " log=fluxdensities2.txt > tmp2.log";
	system "rm tmp2.log fluxdensities2.txt";
	system $pscmd;
	my $scadens = &readdata("fluxdensities2.txt");
	my $dsum = 0.0;
	my $dnum = 0;
	my $defect = 0.0;
	for (my $j = 0; $j <= $#{$vecdens}; $j++) {
	    if ($vecdens->[$j]->[1] > 0) {
		$dsum += ($scadens->[$j]->[1] / $vecdens->[$j]->[1]);
		$dnum++;
	    }
	}
	if ($dnum > 0) {
	    $defect = sprintf "%.1f", (($dsum / $dnum) - 1.0) * 100.0;
	}
	push @{$ro{'defect'}}, {
	    'stokes' => 'I', 'defect' => ($defect * 1.0)
	};
    }

    chdir "../..";
    
    return \%ro;

}

sub load_configurations {
    if ($#configuration_strings > -1) {
	# Already loaded and cached.
	return;
    }

    open(ARRAYS, "/home/jstevens/usr/share/configuration_stations.file");
    while(<ARRAYS>) {
	chomp;
	push @configuration_strings, $_;
    }
    close(ARRAYS);

    return;
}

sub readdata {
    my $fname = shift;

    my @d;
    open(F, $fname);
    while(<F>) {
	chomp;
	my @e = split(/\s+/);
	my $frq = sprintf "%.3f", $e[1];
	my $fld = sprintf "%.3f", $e[2];
	push @d, [ $frq * 1.0, $fld * 1.0 ];
    }
    close(F);

    my @sd = sort { $a->[0] <=> $b->[0] } @d;
    
    return \@sd;
}

sub readlog {
    my $logname = shift;

    my %rv = (
	'fitCoefficients' => [],
	'alphaCoefficients' => [],
	'alphaReference' => { 'fluxDensity' => 0, 'frequency' => 0 },
	'fitScatter' => 0,
	'mode' => "",
	'stokes' => ""
	);
    
    open(L, $logname);
    while(<L>) {
	chomp(my $line = $_);
	my @e = split(/\s+/, $line);
	if ($e[0] eq "Coeff:") {
	    for (my $i = 1; $i < $#e; $i++) {
		push @{$rv{'fitCoefficients'}}, $e[$i] * 1.0;
	    }
	    push @{$rv{'fitCoefficients'}}, $e[$#e];
	} elsif ($e[0] eq "MFCAL") {
	    $rv{'alphaReference'}->{'fluxDensity'} = &removecomma($e[2]) * 1.0;
	    $rv{'alphaReference'}->{'frequency'} = &removecomma($e[3]) * 1.0;
	} elsif ($e[0] eq "Alpha:") {
	    for (my $i = 1; $i <= $#e; $i++) {
		push @{$rv{'alphaCoefficients'}}, $e[$i] * 1.0;
	    }
	} elsif ($e[0] eq "Scatter") {
	    my $fsc = sprintf "%.3f", $e[3];
	    $rv{'fitScatter'} = $fsc * 1.0;
	} elsif ($e[3] eq "Coefficients:") {
	    $rv{'mode'} = lc($e[0]);
	} elsif ($e[0] eq "Stokes") {
	    $rv{'stokes'} = $e[1];
	}
    }
    close(L);

    return \%rv;
}

sub removecomma {
    my $s = shift;

    $s =~ s/\,//g;

    return $s;
}

sub determine_array {
    my $set = shift;

    # Load the required data.
    &load_configurations();

    # Get the positions of the antennas.
    my $cmd = "uvlist vis=".$set." options=full,array";
    my @cout = &execute_miriad($cmd);

    my %antpos = (
	'x' => [], 'y' => [], 'z' => [] );
    for (my $i=0; $i<=$#cout; $i++) {
	my @els = split(/\s+/, $cout[$i]);
	if ($els[1] > 0 && $els[1] < 7) {
	    $antpos{'x'}->[$els[1] - 1] = $els[2];
	    $antpos{'y'}->[$els[1] - 1] = $els[3];
	    $antpos{'z'}->[$els[1] - 1] = $els[4];
	}
    }

    # Adjust to make antenna 6 the reference.
    for (my $i=0; $i<6; $i++) {
	$antpos{'x'}->[$i] -= $antpos{'x'}->[5];
	$antpos{'y'}->[$i] -= $antpos{'y'}->[5];
	$antpos{'z'}->[$i] -= $antpos{'z'}->[5];
	$antpos{'x'}->[$i] *= -1;
	$antpos{'y'}->[$i] *= -1;
	$antpos{'z'}->[$i] *= -1;
    }

    # The station interval is 15.3m.
    my $station_interval = 15.3;
    my @array_stations;
    for (my $i=0; $i<6; $i++) {
	my $ew_offset = floor(($antpos{'y'}->[$i] / $station_interval) + 0.5) + 392;
	my $ns_offset = floor(($antpos{'x'}->[$i] / $station_interval) + 0.5) + 0;
	if ($ns_offset == 0) {
	    push @array_stations, "W".$ew_offset;
	} else {
	    push @array_stations, "N".$ns_offset;
	}
    }

    # Find the best match to the array.
    my $max_matches = 0;
    my $match_array = '';
    for (my $i=0; $i<=$#configuration_strings; $i++) {
	my $curr_match_count = 0;
	for (my $j=0; $j<=$#array_stations; $j++){
	    if ($configuration_strings[$i] =~ /$array_stations[$j]/){
		$curr_match_count++;
	    }
	}
	if ($curr_match_count > $max_matches){
	    $max_matches = $curr_match_count;
	    $match_array = $configuration_strings[$i];
	}
    }

    return $match_array;
}

sub measure_flagging_statistic {
    my $set = shift;
    my $stime = shift;
    my $etime = shift;
    
    my $cmd = "uvfstats vis=".$set." mode=channel options=absolute,unflagged ".
	"\"select=time(".&datetime2miriad($stime).",".
	&datetime2miriad($etime).")\"";
    my $nchans = 0;
    my $nflagged = 0;
    my @cout = &execute_miriad($cmd);

    my $strt = 0;
    for (my $i = 0; $i <= $#cout; $i++) {
	my $line = $cout[$i];
	$line =~ s/^\s+//;
	my @els = split(/\s+/, $line);
	if ($strt == 1) {
	    $nchans += 1;
	    if ($els[1] < 15) {
		# This is a completely flagged channel basically.
		$nflagged += 1;
	    }
	} else {
	    if ($els[0] eq "-------") {
		$strt = 1;
	    }
	}
    }

    my $ffrac = $nflagged / $nchans;
    my $sffrac = sprintf("%.2f", $ffrac);
    $ffrac = $sffrac * 1.0;
    return $ffrac;
}

sub measure_closure_phase {
    my $set = shift;
    my $stime = shift;
    my $etime = shift;
    
    my $closurelog = "closure_log.txt";
    my $cmd = "closure vis=".$set." stokes=i device=/null options=log ".
	"\"select=time(".&datetime2miriad($stime).",".
	&datetime2miriad($etime).")\"";
    if (-e $closurelog) {
	system "rm -f ".$closurelog;
    }
    my @cout = &execute_miriad($cmd);

    my %rv = (
	'closure_phase' => { 'theoretical_rms' => 0, 
			     'measured_rms' => 0,
			     'average_value' => -999 }
	);
    for (my $i=0; $i<=$#cout; $i++) {
	my @els = split(/\s+/, $cout[$i]);
	if ($els[0] eq "Actual") {
	    $rv{'closure_phase'}->{'measured_rms'} = $els[$#els];
	} elsif ($els[0] eq "Theoretical") {
	    $rv{'closure_phase'}->{'theoretical_rms'} = $els[$#els];
	}
    }

    if (-e $closurelog) {
	my @pvals;
	open(F, "closure_log.txt");
	while(<F>) {
	    chomp;
	    my @els = split(/\s+/);
	    if ($#els == 2) {
		push @pvals, $els[2];
	    }
	}
	close(F);
	$rv{'closure_phase'}->{'average_value'} = 
	    sprintf "%.3f", &average(@pvals);
    }

    return %rv;
}

sub average {
    my @a = @_;
    
    my $s = 0;
    for (my $i=0; $i<=$#a; $i++) {
	$s += $a[$i];
    }
    $s /= ($#a + 1);

    return $s;
}

sub execute_miriad {
    my ($miriad_command)=@_;

    my @miriad_output;
    print "executing $miriad_command\n";
    open(MIRIAD,"-|")||exec $miriad_command." 2>&1";
    while(<MIRIAD>){
	chomp;
	my $line=$_;
	push @miriad_output,$line;
    }
    close(MIRIAD);

    return @miriad_output;
}

sub datetime2miriad {
    my $dtime = shift;

    my $outtime = sprintf("%02d%3s%02d:%02d:%02d:%02d.%1d", $dtime->year() - 2000,
			  uc($dtime->month_abbr()), $dtime->day(), $dtime->hour(),
			  $dtime->minute(), $dtime->second(), $dtime->nanosecond * 1e-8);
    
    return $outtime;
}

sub miriad2datetime {
    my $mtime = shift;
    
    #    print $mtime."\n";
    if ($mtime =~ /^(..)(...)(..)\:(..)\:(..)\:(..)\.(.)$/) {
	my $dtime = DateTime->new(
	    year => 2000 + $1,
	    month => &monthname2number($2),
	    day => $3,
	    hour => $4,
	    minute => $5,
	    second => $6,
	    nanosecond => ($7 * 1e8),
	    time_zone => 'UTC'
	    );
	return $dtime;
    }
    return 0;
}

sub monthname2number {
    my $mname = shift;

    if ($mname =~ /^jan/i) {
	return 1;
    } elsif ($mname =~ /^feb/i) {
	return 2;
    } elsif ($mname =~ /^mar/i) {
	return 3;
    } elsif ($mname =~ /^apr/i) {
	return 4;
    } elsif ($mname =~ /^may/i) {
	return 5;
    } elsif ($mname =~ /^jun/i) {
	return 6;
    } elsif ($mname =~ /^jul/i) {
	return 7;
    } elsif ($mname =~ /^aug/i) {
	return 8;
    } elsif ($mname =~ /^sep/i) {
	return 9;
    } elsif ($mname =~ /^oct/i) {
	return 10;
    } elsif ($mname =~ /^nov/i) {
	return 11;
    } elsif ($mname =~ /^dec/i) {
	return 12;
    }
    return -1;
}
