#!/bin/bash
# Some helpers
bold=$(tput bold)
normal=$(tput sgr0)

# Test number (folder name)
name_test="003"
name_prof_file="profile.dat"

# Path to the input files
path_to_input_files="../bin/bin/data"
name_of_input_files="calculation_input.dat"

# Re-compile
echo ''
echo "${bold}Building...${normal}"
cd ../build
./configure.sh
cd ../profiling

# Create a dummy of the output file
echo "# Elapsed CPU time of the main loop, based on the provided test cases (see tests/data/(0,1,2) folders)" > $name_test/$name_prof_file
echo "# -O5 option" >> $name_test/$name_prof_file
echo "#test   time, [s]" >> $name_test/$name_prof_file

# Process reference data
ref_time1=$(sed -n 4p 000/profile.dat | awk '{print $2;}')
ref_time2=$(sed -n 5p 000/profile.dat | awk '{print $2;}')
ref_time3=$(sed -n 6p 000/profile.dat | awk '{print $2;}')

echo ''
echo "${bold}Running...${normal}"

# Run the code, pase the output and write it to the "profile" file
echo "  Test #1"
../bin/bin/mesoapprox.x $path_to_input_files/0/$name_of_input_files > bpec.dat
time_test1=$(cat bpec.dat | grep "Elapsed CPU-time" | awk '{print $3;}')
echo '1' $time_test1 >> $name_test/$name_prof_file
improvement=$(echo "scale=2; ($ref_time1-$time_test1)/$ref_time1*100" | bc)
echo "         ...done in "$time_test1 "seconds. Improvement: "$improvement '%'


echo "  Test #2"
../bin/bin/mesoapprox.x $path_to_input_files/1/$name_of_input_files > k2nnc2.dat
time_test2=$(cat k2nnc2.dat | grep "Elapsed CPU-time" | awk '{print $3;}')
echo '2' $time_test2 >> $name_test/$name_prof_file
improvement=$(echo "scale=2; ($ref_time2-$time_test2)/$ref_time2*100" | bc)
echo "         ...done in "$time_test2 "seconds. Improvement: "$improvement '%'


echo '  Test #3'
../bin/bin/mesoapprox.x $path_to_input_files/2/$name_of_input_files > k3nnc2.dat
time_test3=$(cat k3nnc2.dat | grep "Elapsed CPU-time" | awk '{print $3;}')
echo '3' $time_test3 >> $name_test/$name_prof_file
improvement=$(echo "scale=2; ($ref_time3-$time_test3)/$ref_time3*100" | bc)
echo "         ...done in "$time_test3 "seconds. Improvement: "$improvement '%'


# Move generated output to the folder
mv bpec.dat $name_test
mv BPEC_Fortran_Theta_vs_Mu.txt $name_test
mv k2nnc2.dat $name_test
mv K2NNC2_Fortran_Theta_vs_Mu.txt $name_test
mv k3nnc2.dat $name_test
mv K3NNC2_Fortran_Theta_vs_Mu.txt $name_test

# Perform the sanity check
echo ''
echo "${bold}Testing...${normal}"
sanity_tests_file=sanity_check.dat
python compare.py > $sanity_tests_file
if grep -q "Failed!" "$sanity_tests_file"; then
  echo "Warning! Some test has failed! See '"$sanity_tests_file"' file)"
else
  echo "  All tests have been passed"
  rm $sanity_tests_file
fi

echo ''
