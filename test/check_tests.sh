#!/bin/bash


H=`pwd`

echo "#####################################"
echo "Checking bhpm.npm...."

cd $H/bhpm.npm/changed_parameters
./check_test.sh
if [ $? -ne 0 ]
then
	echo "bhpm.npm(1)...FAILED <----"
else
echo "bhpm.npm(1)...PASSED"
fi


cd $H/bhpm.npm/default_parameters
./check_test.sh
if [ $? -ne 0 ]
then
	echo "bhpm.npm(2)...FAILED <----"
else
echo "bhpm.npm(2)...PASSED"
fi


cd $H/bhpm.npm/initial_values/
./check_test.sh
if [ $? -ne 0 ]
then
	echo "bhpm.npm(3)...FAILED <----"
else
echo "bhpm.npm(3)...PASSED"
fi


cd $H/bhpm.npm/single_chain/changed_parameters/
./check_test.sh
if [ $? -ne 0 ]
then
	echo "bhpm.npm(4)...FAILED <----"
else
echo "bhpm.npm(4)...PASSED"
fi


cd $H/bhpm.npm/single_chain/default_parameters/
./check_test.sh
if [ $? -ne 0 ]
then
	echo "bhpm.npm(5)...FAILED <----"
else
echo "bhpm.npm(5)...PASSED"
fi


cd $H/bhpm.npm/single_chain/initial_values/
./check_test.sh
if [ $? -ne 0 ]
then
	echo "bhpm.npm(6)...FAILED <----"
else
echo "bhpm.npm(6)...PASSED"
fi

echo "#####################################"
echo "Checking bhpm.pm...."

cd $H/bhpm.pm/changed_pm_weights/
./check_test.sh
if [ $? -ne 0 ]
then
	echo "bhpm.pm(1)...FAILED <----"
else
echo "bhpm.pm(1)...PASSED"
fi



cd $H/bhpm.pm/changed_sim_params/
./check_test.sh
if [ $? -ne 0 ]
then
	echo "bhpm.pm(2)...FAILED <----"
else
echo "bhpm.pm(2)...PASSED"
fi



cd $H/bhpm.pm/default_parameters/
./check_test.sh
if [ $? -ne 0 ]
then
	echo "bhpm.pm(3)...FAILED <----"
else
echo "bhpm.pm(3)...PASSED"
fi


cd $H/bhpm.pm/initial_values/
./check_test.sh
if [ $? -ne 0 ]
then
	echo "bhpm.pm(4)...FAILED <----"
else
echo "bhpm.pm(4)...PASSED"
fi

cd $H/bhpm.pm/single_chain/changed_pm_weights/
./check_test.sh
if [ $? -ne 0 ]
then
	echo "bhpm.pm(5)...FAILED <----"
else
echo "bhpm.pm(5)...PASSED"
fi

cd $H/bhpm.pm/single_chain/changed_sim_params/
./check_test.sh
if [ $? -ne 0 ]
then
	echo "bhpm.pm(6)...FAILED <----"
else
echo "bhpm.pm(6)...PASSED"
fi

cd $H/bhpm.pm/single_chain/default_parameters/
./check_test.sh
if [ $? -ne 0 ]
then
	echo "bhpm.pm(7)...FAILED <----"
else
echo "bhpm.pm(7)...PASSED"
fi

cd $H/bhpm.pm/single_chain/initial_values/
./check_test.sh
if [ $? -ne 0 ]
then
	echo "bhpm.pm(8)...FAILED <----"
else
echo "bhpm.pm(8)...PASSED"
fi



exit 0
