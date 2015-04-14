#!/bin/bash

# Set coordinates type coor (1 for Mercator, 2 for lat-lon)
coor=1
sed -i "s/^coor *= *[0-9]\+/coor=$coor/" control_file.txt

# Run baro_v1 for all following vals of nph
for nph in 720 1080 1440 2160 2880 #3600
do
sed -i "s/^nph *= *[0-9]\+/nph=$nph/" control_file.txt
unbuffer ./baro_v1 2>&1 | tee -a logs/itd_runs2/log_${coor}_${nph}.txt
done

# Set coordinates type coor (1 for Mercator, 2 for lat-lon)
coor=2
sed -i "s/^coor *= *[0-9]\+/coor=$coor/" control_file.txt

# Run baro_v1 for all following vals of nph
for nph in 720 1080 1440 2160 2880 3600 4320 #5040
do
sed -i "s/^nph *= *[0-9]\+/nph=$nph/" control_file.txt
unbuffer ./baro_v1 2>&1 | tee -a logs/itd_runs2/log_${coor}_${nph}.txt
done





sed -i "s/^itd_coeff *= *[0-9]\+/itd_coeff=1.4/" control_file.txt







# Set coordinates type coor (1 for Mercator, 2 for lat-lon)
coor=1
sed -i "s/^coor *= *[0-9]\+/coor=$coor/" control_file.txt

# Run baro_v1 for all following vals of nph
for nph in 720 1080 1440 2160 2880 #3600
do
sed -i "s/^nph *= *[0-9]\+/nph=$nph/" control_file.txt
unbuffer ./baro_v1 2>&1 | tee -a logs/itd_runs2/log_${coor}_${nph}.txt
done

# Set coordinates type coor (1 for Mercator, 2 for lat-lon)
coor=2
sed -i "s/^coor *= *[0-9]\+/coor=$coor/" control_file.txt

# Run baro_v1 for all following vals of nph
for nph in 720 1080 1440 2160 2880 3600 4320 #5040
do
sed -i "s/^nph *= *[0-9]\+/nph=$nph/" control_file.txt
unbuffer ./baro_v1 2>&1 | tee -a logs/itd_runs2/log_${coor}_${nph}.txt
done
