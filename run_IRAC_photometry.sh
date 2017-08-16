#!/bin/bash
##############################################################
# WISPIPE - IRAC PHOTOMETRY
# Ivano Baronchelli, 2017

# Call in bash using: > source ./run_IRAC_photometry.py FIELD >& log#.log
# Example:
# if working on Par68:
# source ./run_IRAC_photometry.py 68 I >& LOG/log68.log
# or
# source ./run_IRAC_photometry.py 68 A >& LOG/log68.log
# Using 'I' or 'A', an interactive or automatic run can be performed
##############################################################

ur_setup

cwd=$(pwd) # current directory
cd $WISPDATA/aXe/Par$1/ # go to working directory
python $IRAC_PHOT/IRAC_photometry.py $1 $2
cd $cwd # return to current directory
