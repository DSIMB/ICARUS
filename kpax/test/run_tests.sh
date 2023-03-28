#!/bin/sh
#
#  run_tests.sh
#
echo "--------------------------------------------------------------------"
echo "Running test #1 - align two structures... "
echo "--------------------------------------------------------------------"

command="kpax d1bhga1.ent d1cs6a3.ent"

echo $command
eval $command
result=$?

if test $result -ne 0; then
   echo "The command \'$command\' returned an error code."
   exit $result
fi

echo "--------------------------------------------------------------------"
echo "Running test #2 - make multiple alignment of three structures... "
echo "--------------------------------------------------------------------"

command="kpax -multi d1bhga1.ent d1cs6a3.ent 3ullA00"

echo $command
eval $command
result=$?

if test $result -ne 0; then
   echo "The command \'$command\' returned an error code."
   exit $result
fi

echo "--------------------------------------------------------------------"
echo "Running test #3 - flexible multiple alignment of three structures... "
echo "--------------------------------------------------------------------"

command="kpax -multi -flex d1bhga1.ent d1cs6a3.ent 3ullA00"

echo $command
eval $command
result=$?

if test $result -ne 0; then
   echo "The command \'$command\' returned an error code."
   exit $result
fi


echo "--------------------------------------------------------------------"
echo "Running test #4 - build a small database... "
echo "--------------------------------------------------------------------"

KPAX_DATABASE=/tmp
export KPAX_DATABASE

command="kpax -build=test d1bhga1.ent d1cs6a3.ent 1qb5D00  3ullA00"

echo $command
eval $command
result=$?

if test $result -ne 0; then
   echo "The command \'$command\' returned an error code."
   exit $result
fi

echo "--------------------------------------------------------------------"
echo "Running test #5 - search a small database... "
echo "--------------------------------------------------------------------"

command="kpax -db=test d1bhga1.ent"

echo $command
eval $command
result=$?

if test $result -ne 0; then
   echo "The command \'$command\' returned an error code."
   exit $result
fi

echo "--------------------------------------------------------------------"

