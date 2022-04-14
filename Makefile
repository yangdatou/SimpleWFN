ARMADILLO_INCLUDE ?= /Users/yangjunjie/work/armadillo-code/include
EIGEN_INCLUDE     ?= /Users/yangjunjie/Downloads/eigen-3.4.0

CXX      ?= g++
CXXFLAGS ?= -Wall -Werror -Wno-sign-compare -Wno-comment -std=c++11 -O2 -I $(ARMADILLO_INCLUDE) -I $(EIGEN_INCLUDE) -DARMA_DONT_USE_WRAPPER -framework Accelerate

# Run a test
test: test_utils
	./bin/test_utils.x

# Compile the main executable
main: ./src/main.cc ./src/utils.cc
	$(CXX) $(CXXFLAGS) $^ -o ./bin/main.x

# Compile the main executable
test_utils: ./src/test/test_utils.cc ./src/utils.cc
	$(CXX) $(CXXFLAGS) $^ -o ./bin/test_utils.x

# Remove automatically generated files
clean :
	rm -rf ./bin/* ./output/*.out
