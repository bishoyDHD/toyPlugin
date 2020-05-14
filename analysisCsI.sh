#!/bin/sh
echo "**********************************************************"
echo "**************| BEGIN ANALYSIS ROUTINE |******************"
echo "**********************************************************"
echo "Select type of analysis (1) or (2)";
echo "  1.) Calibration analysis";
echo "  2.) Clustering analysis";
read antype
if [ $antype == 1 ]; then
    nice ./bin/visco ../recipe/CsI/calibCsI.xml $path_cooked/calib65timeRuns.root $path_marinate/mar65TimeRuns.root
    #nice ./bin/visco ../recipe/CsI/calibCsI.xml $path_cooked/calib47Runs.root $path_marinate/marinCalib47Run.root
else
    if [ $antype == 2 ]; then
        echo "Enter run number to be analyzed"
        read runNo
        nice ./bin/visco ../recipe/CsI/marinate.xml $path_cooked/cookedRun${runNo}.root $path_marinate/marinRun${runNo}.root
        #nice ./bin/cooker ../recipe/CsI/marinate.xml $path_cooked/cookedRun${runNo}.root $path_marinate/marinRun${runNo}.root
    fi
fi
