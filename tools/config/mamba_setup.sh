#!/bin/bash

#This script construct conda environments from a list of yaml files given as parameters 
mamba env list > env_list

for yaml in "$@"
do
	echo $yaml
	env=`sed -e 's/\.yaml$//g' -e 's/^.*\///g' <<< $yaml`
	echo $env
	if grep -q "$env"_env env_list; then
		mamba remove --name "$env"_env --all
	fi
	mamba env create --name "$env"_env --file "$yaml" 
done

rm env_list
