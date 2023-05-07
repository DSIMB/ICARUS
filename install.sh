#!/bin/bash

# Run SWORD once to compile its dependency: DSSP
./bin/SWORD &>/dev/null
if [ -f ./bin/SWORD_bin/Dssp/dsspcmbi ]
then
    echo "Installed SWORD"
else
    echo "Error: unable to compile necessary dependancies for SWORD"
    exit 1
fi

echo "All set, good to go !"
