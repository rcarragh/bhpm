echo "#####################################"
echo "Checking bhpm.npm...."

cd bhpm.npm\changed_parameters
CALL check_test.bat
cd..\..

cd bhpm.npm\default_parameters
CALL check_test.bat
cd..\..

cd bhpm.npm\initial_values\
CALL check_test.bat
cd..\..

cd bhpm.npm\single_chain\changed_parameters
CALL check_test.bat
cd..\..\..

cd bhpm.npm\single_chain\default_parameters
CALL check_test.bat
cd..\..\..

cd bhpm.npm\single_chain\initial_values\
CALL check_test.bat
cd..\..\..

echo "#####################################"
echo "Checking bhpm.pm...."

cd bhpm.pm\changed_pm_weights\
CALL check_test.bat
cd..\..

cd bhpm.pm\changed_sim_params\
CALL check_test.bat
cd..\..

cd bhpm.pm\default_parameters\
CALL check_test.bat
cd..\..

cd bhpm.pm\initial_values\
CALL check_test.bat
cd..\..

cd bhpm.pm\single_chain\changed_pm_weights\
CALL check_test.bat
cd..\..\..

cd bhpm.pm\single_chain\changed_sim_params\
CALL check_test.bat
cd..\..\..

cd bhpm.pm\single_chain\default_parameters\
CALL check_test.bat
cd..\..\..

cd bhpm.pm\single_chain\initial_values\
CALL check_test.bat
cd..\..\..
