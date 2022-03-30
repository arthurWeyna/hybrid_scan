#!/bin/bash

##This script should help in creating necessary conda environments, and in setting up the Snakefile


main=`dirname "$(readlink -f "$0")"`


datadir=$main/data
echo "Where do you want to store data files produced or used by this pipeline?"
echo "Default is "$datadir
read -p "Give another path or press enter to use default: " input
if [ ! -z "$input" ]; 
then
datadir=`realpath $input`
fi
mkdir -p $datadir
awk -v FS="=" -v OFS="=" -v datad=$datadir -v toold=$main/tools -v ucef=$main/uce_list -v buscof=$main/busco_list -v q="\"" '/^DATADIR=/{$2=q""datad""q}; /^TOOLDIR=/{$2=q""toold""q}; /^UCETARGETS=/{$2=q""ucef""q}; /^BUSCOTARGETS=/{$2=q""buscof""q}; {print}' $main/Snakefile > Snakefile2

mv Snakefile2 Snakefile

echo "Snakefile is set. Constructing conda environments..."

if ! command -v conda &> /dev/null
then
	echo "Conda could not be found. Get the installer at https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html"
	exit
fi

if ! command -v mamba &> /dev/null
then
	echo "Mamba could not be found. Install it with 'conda install mamba -n base -c conda-forge'."
	exit
fi

$main/tools/config/mamba_setup.sh $main/tools/config/*.yaml

