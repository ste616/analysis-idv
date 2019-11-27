TIME-VARIABLE CALIBRATION AND ANALYSIS FOR ATCA DATA
Tutorial by Jamie Stevens, 2019 November
Jamie.Stevens@csiro.au

This tutorial covers how to use the time-variable calibration tools available in the GitHub
repository:

https://github.com/ste616/analysis-idv

It does this by doing a reduction and analysis of an epoch of C2914 data. This data was looking
at the intra-hour variable (IHV) source J1325-111 discovered by the C2914 team.

1. Getting the data.

Use ATOA to access the now public-domain C2914 data from February 22 2017.

https://atoa.atnf.csiro.au

You will want to get the files 2017-02-22_1147.C2914, 2017-02-22_1542.C2914, and
2017-02-22_1937.C2914. Put all of these files in an otherwise empty directory. In total, this
is 9.3 GB of data.

2. Prepare the data.

You should first do a pretty standard reduction, assuming that 1934-638 is the flux density
and bandpass calibrator, and that the source 1349-145 is a gain calibrator. The idea here is
that we should establish that our calibration strategy works sufficiently well for the gain
calibrator. We do that by determining the flux density model of this calibrator with the task
uvfmeas, and getting a smooth, believable, high-quality result.

There is nothing special about the reduction you should perform here, and everything is what
the Miriad or ATCA users guides would recommend. Included with this tutorial is the shell
script "data_preparation.sh" which will do this standard reduction and output a plot of the
flux density model for 1349-145. Please examine this script for guidance on your own reductions.

3. Prepare the code.

You should obtain the analysis-idv package from GitHub in the normal way, and install it as
follows:

cd analysis-idv/cabb_json
python setup.py install

Should this not work, and break on the installation of other dependencies, use pip to install
the following packages:

astropy
docopt
dominate
ephem
matplotlib
mirpy
scipy

Ensure that the two Python scripts "timevariation_cal.py" and "source_investigation.py" are
in your executables path.

4. Do the time-variable calibration.

In the C2914 data directory, after you have run the calibration script, or done your own
calibration, you should now run the timevariation_cal.py script.

It should generally only need a single argument, that being the output from the Miriad
atlod stage. In this particular case, if you used the "data_preparation.sh" script,
you should give the command:

timevariation_cal.py c2914_2017-02-22.uv

The script works by assuming the following:
1. There is one source which has accurate bandpass/flux density and polarisation calibration
   tables.
2. That all other sources of interest can be self-calibrated over time, and flux density
   can be obtained using the source from assumption 1.

You can see how this is exactly the same as the normal process used for 1349-145 in section
2 above. The script:
- copies over the calibration tables from the reference source (which by default is 1934-638)
  to the source of interest,
- does a gain calibration to determine the time-dependent gains, assuming the leakages from
  the reference source are accurate,
- rescales the gains to the flux density standard of the reference source.

For 1349-145, which is assumed to have a non-varying flux density over the timescale of the
observation, this process is easily done manually, as evidenced by the preparation script.
All the time-variable calibration script does is handle the repetitive splitting of the
data into time chunks, and calibrating them in the same way.

After the calibration, the flux densities are measured for each channel from the time
chunk that was just calibrated, along with other parameters like the flux density model
(which gives the spectral index and curvature), the closure phase and the defect. All of this
data is then saved into a JSON file, labelled with time as MJD. The script does this for
all sources that it finds in the top-level file (except for the reference source); if there
is a source which had been assumed to be non-varying, the output of the script can be used to
check that assumption.

For the IHV source J1325-111, the observations were made in two frequency configurations,
and these were alternated every 2 minutes. The first configuration had CABB central frequencies
of 4928 and 6720 MHz, while the second configuration had 8512 and 10304 MHz. It would thus
seem natural to produce a flux density every 2 minutes. By default, this is what the time-variable
calibration script does, but this can be easily changed with the switch "--segment", which
accepts as an argument the number of minutes included for each flux density measurement. For
example, if you would rather (for sensitivity purposes perhaps) that a measurement include 10
minutes of data, you could run the script like so:

timevariation_cal.py --segment 10 c2914_2017-02-22.uv

There is an associated switch "--slop", which allows the script to include slightly more data
than the segment, if it makes sense to do so. By default, slop is 1 minute. For example, if you
had observed 31 minutes of data on your source, and set the segment to be 10 minutes, the
script would split the data into three chunks, two with 10 minutes of data, and the last with 11 minutes
(since slop of 1 minute allows the script to do this). If you had have observed the source for 32
minutes instead, and kept slop at 1 minute, the script would be forced to split the data into
four chunks, three with 10 minutes of data, and the fourth with 2 minutes (and thus be much less
sensitive). To keep the three chunks in this case, you could:

timevariation_cal.py --segment 10 --slop 2 c2914_2017-02-22.uv

At the end of the execution of the script, you will have two extra directories: "data" and "plots".
The "plots" directory will be full of PS files showing the spectra of each source at each
time chunk, as they came out of uvfmeas during the measurement stage. There will be one plot for
the vector-averaged data, and one for the scalar-averaged data, for each chunk, for each source.
There is no need to keep these plots, but they can be handy for diagnosing problems with your dataset,
although it can be unweildy to go through them all.

The "data" directory will contain a single JSON file for each source that the script found
and measured. The JSON is in the top-level an object, and the members of that object are described
below.
{ "sourceName": the name of the source, which should be the same as the file name (string)
  "rightAscension": the J2000 right ascension of the phase centre of the observation (string, sexagesimal hours)
  "declination": the J2000 declination of the phase centre of the observation (string, sexagesimal degrees)
  "arrayConfiguration": the name of the array configuration during the observation (string)
  "timeSeries": [ an array of all the measurements made, one per time chunk, each element being an object with
  		  the members listed below
    "closurePhase": [ an array of all the measurements made with the closure task on this chunk, one per
    		      IF, each element being an object with the members listed below
      "theoretical_rms": the RMS thermal noise expected for this chunk (float, Jy)
      "measured_rms": the RMS thermal noise actually measured for this chunk (float, Jy)
      "average_value": the average closure phase over all baselines and for the whole chunk (float, deg)
      "IF": the central frequency of this band (int, MHz)
    ]
    "centreFreqs": [ ] an array of all the centre frequencies in this configuration (array of int, MHz)
    "flaggedFraction": [ ] an array of the fraction of data flagged in this chunk in each band (array of float;
    		           has same order as "centreFreqs")
    "defect": [ an array of all the measurements made of the defect, one per Stokes parameter, each element
    	        being an object with the members listed below
      "stokes": the Stokes parameter (string)
      "defect": the measured defect (float)
    ]
    "hourAngleRange": { an object describing the hour-angle limits of this chunk, with the members
      "low": the earliest hour angle of this chunk (float, hours)
      "high": the latest hour angle of this chunk (float, hours)
    }
    "mjdRange": { an object describing the MJD limits of this chunk, with the members
      "low": the earliest MJD of this chunk (float)
      "high": the latest MJD of this chunk (float)
    }
    "fluxDensityData": [ an array of all the measurements made of the flux density of this source in this
    		         chunk, one per Stokes parameter per averaging method, each element being an object
			 with the members listed below
      "stokes": the Stokes parameter (string)
      "mode": the averaging mode used (string, one of "vector" or "scalar")
      "data": [ an array of all the measurements, one per unflagged channel, each element being a two
      	        element array
	[ frequency (float, GHz), flux density (float, Jy) ]
      ]
    ]
  ]
}

5. Use the JSON to make things humans can understand.

The "source_investigation.py" code uses the analysis_idv library to read in the JSON files,
produce plots, and calculate autocorrelation timescales.

To run this code, you need to give it a source name (it works only on one source at a time).
For example, to analyse the J1325-111 source in this dataset, you might only need to run:

source_investigation.py j1325-111

The script needs to know where to get the JSON data for this source. By default, it searches
anywhere below the current directory for anything called j1325-111.json (or whatever source
you asked it to look for). If multiple files are found, it will treat them all as separate
epochs. If you want to find the data in a directory different to the working directory, use
the switch "--root-dir", with its argument being the top-level directory containing all the
data; the code only accepts one directory as its root directory.

For example, if all your data is under the directory /home/jbloggs/timevariabledata/, then
give the command:

source_investigation.py --root-dir /home/jbloggs/timevariabledata j1325-111

The code will output several plots, data files, and make some web pages, and by default all
of these will be output into the working directory. For neatness, it may be better to specify
a target directory for each, with the switches "--plot-dir" for each of the plots and data
files, and "--web-dir" for the web pages (and for the plots that appear on the web pages). For
example, to output the non-web plots into the "graphics" directory and the web pages and plots
into the "web" directory:

source_investigation.py --plot-dir graphics --web-dir web j1325-111

If the directories don't exist, the code will create them for you, if it can.

There are a few other switches available, which are best understood by describing the plots
and other files that this code will create. These will be given as <path>/<filename>, where
<path> will be where they end up at the end of the process, and "graphics" will be the path
for the non-web outputs, and "web" will be the path for the web outputs; <src> will be whatever
the source name is, <mjd> is the rounded MJD of the epoch, <year> is the year of the epoch,
<doy> is the day of the year of the epoch, and <freq> is the frequency being analysed.

web/<src>/daily_<src>_dynamic_mjd<mjd>_doy<year>-<doy>.png: this plot is a dynamic spectrum
of a single epoch. It has time on the x-axis, frequency on the y-axis, and a colour scaling
for the flux density. It is a waterfall plot, with the chunk spacing as the resolution in
the time direction, while the resolution in the frequency direction can be controlled by the
user with the "--spectral-averaging" switch. By default, this is 1 MHz, but it can be set larger
to increase the signal-to-noise ratio of each measurement, and to reduce the amount of time
that the code needs to produce this plot. For example, to average to 16 MHz resolution:

source_investigation.py --spectral-averaging 16 --plot-dir graphics --web-dir web j1325-111

graphics/full_<src>_timeSeries.png: this plot is the flux density vs time for all the files
that were included in the execution. There will be several lines on this plot, one for each
frequency that is analysed, and each point on that line is one of the time chunk measurements.
The user can specify the frequencies to analyse in one of two ways. The first is to specify
them directly, with multiple uses of the "--frequency" switch. For example, to ask the code
to analyse the frequencies 5604 MHz and 9133 MHz:

source_investigation.py --frequency 5604 --frequency 9133 --plot-dir graphics --web-dir web j1325-111

As many frequencies as desired can be given. The second way is to specify a comb of analysis
frequencies using the three switches "--low-frequency", "--high-frequency" and
"--frequency-interval". For example, to make the code look at the frequencies 5000, 5500, 6000,
6500 and 7000 MHz:

source_investigation.py --low-frequency 5000 --high-frequency 7001 --frequency-interval 500 \
                        --web-dir web --plot-dir graphics j1325-111

The high frequency is specified slightly higher than 7000, because the the frequency generator
doesn't include the end point (think of the behaviour of Python's range function).

The actual frequencies that are analysed may not necessarily be those you ask for though. The code
will find the frequency that is closest to your desired frequency that occurs frequently enough
to make sense. It does this so that you don't need to worry about accidentally selecting a channel
which has been heavily flagged due to RFI. The frequency chosen will also depend on the setting
for the "--spectral-averaging" switch, as described for the previous plot. So it may be that should
you give 16 MHz for the spectral averaging, and ask for frequencies between 5000 and 10000 MHz,
at a 500 MHz interval, you may end up with the frequencies 5008, 5504, 6000, 6496, 7008, 7504,
8000, 8496, 9008 and 9504 MHz (these are central frequencies for the averaged channel that was
selected in each case). There will be a legend on this plot which tells you which frequencies
were selected.

graphics/animation_frames/daily_<src>_mjd<mjd>_doy<year>-<doy>_animation_frame*.png: these
plots show how the whole spectrum evolves over the course of the epoch, as a series of
plots which can be assembled as animation frames, should you wish. In each frame, the current
time chunk is shown as a red line, with the MJD of the chunk shown in the top left of the frame.
All old chunks are shown in grey, and the grey colour lightens as the age of the chunk gets
larger.

web/<src>/daily_<src>_timeSeries_mjd<mjd>_doy<year>-<doy>.png: this plot is very similar to
the "full_<src>_timeSeries.png" plot as described above, except this one is only for a single
epoch, and doesn't have a frequency legend on it. The lines it has are exactly the same frequencies
as the former plot.

web/<src>/daily_<src>_acf_mjd<mjd>_doy<year>-<doy>_f<freq>.png: these plots show the autocorrelation
lag spectra for each frequency that was analysed. The correlation strength should always be 1
at lag 0 (of course, it is an autocorrelation). The strength decreases with increasing lag, and a
Gaussian is fit to the correlation strength (if possible) and plotted here. The full width at 1/e
is plotted as a horizontal dashed line, and range of the data range used for fitting the Gaussian is shown
as two vertical solid lines. There is a pink-shaded region, and all data points that lie here
are considered uncorrelated.

web/<src>/<src>_doy<year>-<doy>_timescaleFrequency.png: these plots show how the timescale
(in hours) varies with respect to frequency, for each epoch.

web/<src>/<src>_<freq>_timescales.png: these plots show how the timescale at each frequency
changes over the epochs, where each epoch is plotted as its day-of-year, and multiple years
will be shown in different colours. This is good if you are expecting an annual cycle for the
variability timescales.

web/<src>/index.html: this web page presents a link to each epoch individually, and shows how
the timescales vary across the epochs.

web/<src>/<src>_doy<year>-<day>_epochinfo.html: these pages show the dynamic spectrum and the
time series plot for each epoch individually, along with the autocorrelation lag plots, and the
plot showing how the timescale varies with frequency.

