#!/bin/bash
URL=http://supernovae.in2p3.fr/salt/lib/exe/fetch.php?media=
VERSION=2.4.2

# download
mkdir -p download
wget ${URL}salt2-4_data.tgz -O download/salt2-4_data.tgz
wget ${URL}snfit-${VERSION}.tar.gz -O download/snfit-${VERSION}.tar.gz
wget http://supernovae.in2p3.fr/sdss_snls_jla/jla_light_curves.tgz -O download/jla_light_curves.tgz

# unpacking data
tar xvzf download/salt2-4_data.tgz
tar xvzf download/jla_light_curves.tgz

# unpack source code
mkdir src
tar xvzf download/snfit-${VERSION}.tar.gz
mv snfit-${VERSION} src/
