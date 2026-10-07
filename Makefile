CXX ?= g++
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra -Iinclude

all: test_fft bench_fft

build:
	mkdir -p build

test_fft: tests/test_fft.cpp src/fft.cpp include/fft.h | build
	$(CXX) $(CXXFLAGS) tests/test_fft.cpp src/fft.cpp -o build/test_fft

bench_fft: bench/bench.cpp src/fft.cpp include/fft.h | build
	$(CXX) $(CXXFLAGS) bench/bench.cpp src/fft.cpp -o build/bench_fft

test: test_fft
	./build/test_fft

bench: bench_fft
	./build/bench_fft

clean:
	rm -rf build

.PHONY: all test bench clean
