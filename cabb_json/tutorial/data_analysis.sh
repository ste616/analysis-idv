#!/bin/bash
# Data analysis script
# Variability reduction tutorial
# Jamie Stevens 2019

echo "Data analysis script"

# Run the time-variable calibrator.
echo "Calibrating and measuring for time-variability"
timevariation_cal.py --segment 2 --slop 1 c2914_2017-02-22.uv

echo "Analysing and plotting"
# The source that should vary.
OPTIONS="--root-dir . --plot-dir graphics --web-dir web --spectral-averaging 16 --low-frequency 5000 --high-frequency 10001 --frequency-interval 500"
source_investigation.py $OPTIONS j1325-111
# And the source that shouldn't.
source_investigation.py $OPTIONS 1349-145

echo "Calibration, measurement and analysis complete."

