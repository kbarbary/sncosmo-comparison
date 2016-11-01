# sncosmo-comparison

Compare fit results of sncosmo and snfit

Setup
=====


Download code and JLA data:

```
./download.sh
```

Build snfit:

```
DIR=$(pwd)
cd src/snfit-2.4.2
export LIBS="-lgfortran"
./configure --prefix=$DIR
make
make install
```

Environment variables (add to `.bashrc` or similar):

```
export LD_LIBRARY_PATH=$DIR/lib:$LD_LIBRARY_PATH
export SALTPATH=$DIR/snfit_data
```

Test:

```
mkdir test
cd test
../bin/snfit ../jla_light_curves/lc-03D4ag.list
```