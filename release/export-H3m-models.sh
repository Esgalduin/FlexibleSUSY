#!/bin/sh

# creates MSSM models with H3m corrections

models="MSSMNoFVH3m,MSSMNoFVatMGUTH3m,NUHMSSMNoFVH3m"
version=$(cat config/version)
output="FlexibleH3m-${version}"

for m in $(echo "$models" | tr ',' ' ') ; do
    ./createmodel --name=$m -f
done

./configure --with-models="${models}" --with-install-dir="$output" --disable-compile
make clean
make -j4
make install-src

tar czf "${output}".tar.gz "${output}"
