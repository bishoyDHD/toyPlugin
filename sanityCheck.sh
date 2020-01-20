#!/bin/sh
#nice ./bin/visco ../recipe/CsI/test.xml $path_decoded/hdw_3vf48.root $path_cooked/calib47Runs.root
echo "Please enter run number to be analyzed (sanity-wise)"
read num;
nice ./bin/visco ../recipe/CsI/sanityCheck.xml $path_cooked/sanityRun${num}.root $path_marinate/cookedSanity${num}.root
