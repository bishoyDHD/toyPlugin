#!/bin/sh
#nice ./bin/visco ../recipe/CsI/test.xml $path_decoded/hdw_3vf48.root $path_cooked/calib47Runs.root
echo "Please enter run number to be sanity checked"
read num;
nice ./bin/visco ../recipe/CsI/test.xml $path_decoded/hdw_${num}.root $path_cooked/testRealTimeRun${num}.root
#---> For target calib files
#nice ./bin/visco ../recipe/CsI/test.xml $path_decoded/hdw_${num}.root:$path_converted/run${num}tgtConv.root $path_cooked/sanityRun${num}.root
