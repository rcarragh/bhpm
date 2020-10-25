#!/bin/bash

cd ./bhpm.npm/changed_parameters/
./cleanup.sh
cd -

cd ./bhpm.npm/default_parameters/
./cleanup.sh
cd -

cd ./bhpm.npm/initial_values/
./cleanup.sh
cd -

cd ./bhpm.npm/single_chain/changed_parameters/
./cleanup.sh
cd -

cd ./bhpm.npm/single_chain/default_parameters/
./cleanup.sh
cd -

cd ./bhpm.npm/single_chain/initial_values/
./cleanup.sh
cd -

cd ./bhpm.pm/changed_pm_weights/
./cleanup.sh
cd -

cd ./bhpm.pm/changed_sim_params/
./cleanup.sh
cd -

cd ./bhpm.pm/default_parameters/
./cleanup.sh
cd -

cd ./bhpm.pm/initial_values/
./cleanup.sh
cd -

cd ./bhpm.pm/single_chain/changed_pm_weights/
./cleanup.sh
cd -

cd ./bhpm.pm/single_chain/changed_sim_params/
./cleanup.sh
cd -

cd ./bhpm.pm/single_chain/default_parameters/
./cleanup.sh
cd -

cd ./bhpm.pm/single_chain/initial_values/
./cleanup.sh
cd -
