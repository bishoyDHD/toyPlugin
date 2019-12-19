#!/bin/sh
echo "**********************************************************"
echo "**************| BEGIN ANALYSIS ROUTINE |******************"
echo "**********************************************************"
echo "Enter run number to be analyzed"
read runNo
nice ./bin/cooker ../trekCkr/recipe/CsI/marinate.xml ${path_cooked}/cookedRun${runNo}.root ${path_marinate}/marinRun${runNo}.root
