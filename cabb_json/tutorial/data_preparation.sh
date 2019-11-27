#!/bin/bash
# Data preparation script
# Variability reduction tutorial
# Jamie Stevens 2019

echo "Data preparation script"
# Check for the necessary files.
RPFITS="2017-02-22_1147.C2914 2017-02-22_1542.C2914 2017-02-22_1937.C2914"
for f in ${RPFITS}
do
	if [ ! -f $f ]; then
		echo "Required RPFITS file $f is not present"
		exit 1
	fi
done
echo "All RPFITS files present, continuing."

# Load the data into Miriad.
MIRFILE=c2914_2017-02-22.uv
RPFIN="$(echo "$RPFITS" | sed 's/[[:space:]]/,/g')"
if [ -d $MIRFILE ]; then
	echo "Removing existing dataset"
	rm -rf $MIRFILE
fi
echo "Loading data into Miriad."
atlod in=$RPFIN out=$MIRFILE options=birdie,rfiflag,xycorr,opcorr,noauto nopcorr=32

# Split the data.
FREQS="4928 6720 8512 10304"
for f in ${FREQS}
do
	if [ -d uvsplit.${f} ]; then
		rm -rf uvsplit.${f}
	fi
done
echo "Splitting the data into frequencies."
uvsplit vis=$MIRFILE "select=-shadow(23)" options=nosource

# Make a directory for the individual sources.
SRCDIR="4cm"
if [ -d $SRCDIR ]; then
	echo "Removing existing source directory ${SRCDIR}"
	rm -rf $SRCDIR
fi
echo "Making source directory ${SRCDIR}"
mkdir $SRCDIR

# Now we operate per frequency.
for f in ${FREQS}
do
	echo "Working on frequency ${f} MHz"
	echo "Automatically flagging RFI"
	pgflag vis=uvsplit.${f} stokes=i,q,u,v flagpar=8,5,5,3,6,3 "command=<b" options=nodisp
	cd $SRCDIR

	echo "Splitting into sources"
	uvsplit vis=../uvsplit.${f}

	echo "Creating bandpass solution"
	mfcal vis=1934-638.${f} interval=0.1

	echo "Copying bandpass to stable calibrator"
	gpcopy vis=1934-638.${f} out=1349-145.${f}

	echo "Calculating gains and polarisation leakages"
	gpcal vis=1349-145.${f} interval=0.1 nfbin=2 options=xyvary,qusolve

	echo "Copy leakages to flux density calibrator"
	gpcopy vis=1349-145.${f} out=1934-638.${f}

	echo "Calculate gains for flux density calibrator"
	gpcal vis=1934-638.${f} interval=0.1 nfbin=2 options=xyvary,qusolve,nopol

	echo "Scaling calibrator gains to flux density standard"
	gpboot vis=1349-145.${f} cal=1934-638.${f}
	cd ..
done

# Make a plot of the stable calibrator to see if it looks properly calibrated.
echo "Plotting flux density of stable calibrator as check"
uvfmeas vis=4cm/1349-145.* stokes=i order=2 options=plotvec,log device=check_1349-145.png/png

