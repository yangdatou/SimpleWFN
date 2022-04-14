ARMADILLO_INCLUDE ?= /Users/yangjunjie/work/armadillo-code/include
EIGEN_INCLUDE     ?= /Users/yangjunjie/Downloads/eigen-3.4.0
CXX      ?= g++
CXXFLAGS ?= -Wall -Werror -Wno-sign-compare -Wno-comment -std=c++11 -O2 -I $(ARMADILLO_INCLUDE) -I $(EIGEN_INCLUDE) -DARMA_DONT_USE_WRAPPER -framework Accelerate

# Run a regression test
test: main
	./bin/main.x ./input/h2o/STO-3G/ > ./output/out.log

# Compile the main executable
main: ./src/main.cc ./src/utils.cc
	$(CXX) $(CXXFLAGS) $^ -o ./bin/main.x

# Remove automatically generated files
clean :
	rm -rvf ./bin/* ./output/*.out