#!/bin/bash

# Run SWORD once to compile its dependency: DSSP
# ./bin/SWORD &>/dev/null
# if [ -f ./bin/SWORD_bin/Dssp/dsspcmbi ]
# then
#     echo "Installed SWORD in ./bin"
# else
#     echo "Error: unable to compile necessary dependancies for SWORD"
#     exit 1
# fi

# # Compile TMalign on Mac or Linux
# unameOut="$(uname -s)"
# case "${unameOut}" in
#     Linux*)     g++ -static -O3 -ffast-math -lm -o ./bin/TMalign ./bin/sources/TMalign.cpp;;
#     Darwin*)    g++ -O3 -ffast-math -lm -o ./bin/TMalign ./bin/sources/TMalign.cpp;;
#     *)          echo "Unknown system: neither a Mac nor Linux."; exit 1;;
# esac

# if [ -f ./bin/TMalign ]
# then
#     echo "Compiled & installed TMalign in ./bin"
# else
#     echo "An error occured, TMalign could not be compiled and installed in ./bin"
#     exit 1
# fi

# Install KPAX
pushd bin >/dev/null 2>&1
wget http://kpax.loria.fr/download/kpax-5.1.3-x64-mint18.3.tgz >/dev/null 2>&1
tar -xf kpax-5.1.3-x64-mint18.3.tgz
rm -f kpax-5.1.3-x64-mint18.3.tgz

if [ -f kpax/bin/kpax5.1.3.x64 ]
then
    ln -s kpax/bin/kpax5.1.3.x64 kpax/bin/kpax5.0.2.x64
    echo "Installed KPAX"
else
    echo "An error occured, KPAX could not be installed in bin/"
    exit 1
fi
popd >/dev/null 2>&1
echo "All set, good to go !"