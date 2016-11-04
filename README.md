# sncosmo-comparison

Compare fit results of sncosmo and snfit for the [SDSS-SNLS Joint Lightcurve Analaysis](http://supernovae.in2p3.fr/sdss_snls_jla/ReadMe.html)

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

Run snfit:

```
mkdir results_snfit
make all  # or, e.g., make 03D4ag
```