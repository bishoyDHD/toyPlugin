#!/bin/sh
echo " *** Entering Cooker: Did you add a new file(s) (yes or no)"
read ans;
ans1="yes";
ans2="no";
if [ $ans1 == $ans ]; then
    cmake ../
    make -j8 && make install
else # case no: No new file(s) added 
   if [ $ans2 == $ans ]; then
      echo " Did you make any modifications (yes or no)?"
      read newAns;
      if [ $ans1 == $newAns ]; then # if existing file(s) added
         make -j8 && make install
      else
         echo " Moving onto running the plugin...";
         echo " without changes of couse";
      fi
   fi
fi
echo " Please enter run number to be analysed: ";
read num;
echo " Would you like to use visco?"
read visCo;
if [ $ans1 == $visCo ]; then
    nice ./bin/visco ../recipe/CsI/clustering.xml $path_decoded/run${num}.root:$path_converted/run${num}conv.root $path_cooked/cookedRun${num}.root
else
    nice ./bin/cooker ../recipe/CsI/clustering.xml $path_decoded/run${num}.root:$path_converted/run${num}conv.root $path_cooked/cookedRun${num}.root
fi
