# `?=` means assigning the value to the variable if it is not already defined
# Set `armadillo` path
ARMADILLO_INCLUDE ?= /Users/yangjunjie/work/armadillo-code/include
# Set `Eigen` path
EIGEN_INCLUDE     ?= /Users/yangjunjie/Downloads/eigen-3.4.0

# Set compiler
CXX      ?= g++
# Set c++ compiler flags
CXXFLAGS ?= -Wall -Werror -Wno-sign-compare -Wno-comment -std=c++11 -O3 -I $(EIGEN_INCLUDE)

# Run the tests
test: main test_utils
	./bin/test_utils.x
	./bin/main.x ./input/h2o/STO-3G/ 5 5 > ./output/h2o-sto-3g.out
	./bin/main.x ./input/h2o/DZ/     5 5 > ./output/h2o-dz.out
	./bin/main.x ./input/ch4/STO-3G/ 5 5 > ./output/ch4-sto-3g.out
	./bin/main.x ./input/h2o/DZP/    5 5 > ./output/h2o-dzp.out
	grep "SCF energy =" ./output/*.out > ./output/tmp.out
	diff ./output/tmp.out ./output/scf.log
	grep "MP2 energy =" ./output/*.out > ./output/tmp.out
	diff ./output/tmp.out ./output/mp2.log
	rm ./output/tmp.out

# Compile the main executable
main: ./src/main.cc ./bin/utils.o
	$(CXX) $(CXXFLAGS) -o ./bin/main.x $^

test_utils: ./src/test/test_utils.cc ./bin/utils.o
	$(CXX) $(CXXFLAGS) -o ./bin/test_utils.x $^

./bin/utils.o: ./src/utils.cc
	$(CXX) $(CXXFLAGS) -o ./bin/utils.o -c $^

test_ccsd_imds: ./src/test/test_ccsd_imds.cc ./bin/ccsd_imds.o ./bin/utils.o
	$(CXX) $(CXXFLAGS) -o ./bin/test_ccsd_imds.x $^

./bin/ccsd_imds.o: ./src/ccsd_imds.cc
	$(CXX) $(CXXFLAGS) -o ./bin/ccsd_imds.o -c $^

# Remove automatically generated files
clean :
	rm -rf ./bin/* ./output/*.log
