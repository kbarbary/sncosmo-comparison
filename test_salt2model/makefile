SRC=../src/snfit-2.4.2/src
LIB=../lib

test-salt2model : test_salt2model.o
	/bin/bash $(SRC)/../libtool --tag=CXX   --mode=link g++ -DUSELAPACK -DgFortran -g -Wall -O3 -R $(shell realpath $(LIB)) -L$(LIB) -o test-salt2model test_salt2model.o -lsnfit -llapack -lblas -llapack -lgfortran

test_salt2model.o : test_salt2model.cc
	g++ -DHAVE_CONFIG_H -I$(SRC) -I$(SRC)/.. -g -Wall -O3 -DUSELAPACK -DgFortran -g -Wall -O3 -MT test_salt2model.o -MD -MP -MF $(SRC)/.deps/snfit.Tpo -c -o test_salt2model.o test_salt2model.cc
