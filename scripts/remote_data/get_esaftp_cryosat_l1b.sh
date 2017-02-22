#!/bin/bash

# server settings (this should not change)
radar_modes=(SAR SIN)
remote_ftp=ftp://science-pds.cryosat.esa.int/SIR_NRT/LATEST/
wget_opt="-nv -c -r -nd"

# Check Arguments
if [ "$#" -lt "3" ]
  then
    echo "invalid call (get_cryosat2_nrt_data.sh \$local_dir \$year \$month [\$day])"
    exit
fi
local_base_dir=$1
year=$2
month=$3
day=$4


local_base_dir=$1
local_dir=$local_base_dir/$year/$month

# Get the login detail for the ESA ftp
echo "Specify login for $remote_ftp"
echo -n user:
read -s user
echo
echo -n passw:
read -s passw
echo

# radar modes are in same directory on remote ftp but need to be in 
# separate dirs for pysiral
for radar_mode in ${radar_modes[@]}; do
    
    # get the local directory
    local_dir=$local_base_dir/SIR_"$radar_mode"_L1/$year/$month
    echo "Downloading (SIR_${radar_mode}_L1) -> $local_dir"

    # transfer files using wget
    remote_selection="CS_NRT__SIR_${radar_mode}_1B_${year}${month}${day}*.TGZ"
    wget $wget_opt -A $remote_selection -P $local_dir --user $user --password $passw $remote_ftp

    # unzip all archives
    echo "unzip TGZ files"
    tgs_file_list=(${local_dir}/*.TGZ)
    for a in "${tgs_file_list[@]}"; do tar -xzf $a -C $local_dir; done

    # clean up
    echo "clean up"
    rm ${local_dir}/*.TGZ
done
