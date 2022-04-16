ARMADILLO_INCLUDE ?= /Users/yangjunjie/work/armadillo-code/include
EIGEN_INCLUDE     ?= /Users/yangjunjie/Downloads/eigen-3.4.0

CXX      ?= g++
CXXFLAGS ?= -Wall -Werror -Wno-sign-compare -Wno-comment -std=c++11 -O3 -I $(ARMADILLO_INCLUDE) -I $(EIGEN_INCLUDE)
# -DARMA_DONT_USE_WRAPPER -framework Accelerate

# Run a test
test: main test_utils
	./bin/test_utils.x
	./bin/main.x ./input/h2o/STO-3G/ 5 5 > ./output/h2o-sto-3g.log
	./bin/main.x ./input/h2o/DZ/     5 5 > ./output/h2o-dz.log
	./bin/main.x ./input/ch4/STO-3G/ 5 5 > ./output/ch4-sto-3g.log
	# ./bin/main.x ./input/h2o/DZP/    5 5 > ./output/h2o-dzp.log

# Compile the main executable
main: ./src/main.cc ./bin/utils.o
	$(CXX) $(CXXFLAGS) -o ./bin/main.x $^

# Compile the main executable
test_utils: ./src/test/test_utils.cc ./bin/utils.o
	$(CXX) $(CXXFLAGS) -o ./bin/test_utils.x $^

./bin/utils.o : ./src/utils.cc
	$(CXX) $(CXXFLAGS) -o ./bin/utils.o -c ./src/utils.cc

# Remove automatically generated files
clean :
	rm -rf ./bin/* ./output/*.log
