CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -pedantic

BUILD := build
LIB   := $(BUILD)/fft.o

.PHONY: all test bench clean

all: $(BUILD)/test_fft $(BUILD)/bench_fft

$(BUILD):
	mkdir -p $(BUILD)

$(LIB): src/fft.cpp src/fft.h | $(BUILD)
	$(CXX) $(CXXFLAGS) -c src/fft.cpp -o $@

$(BUILD)/test_fft: tests/test_fft.cpp $(LIB)
	$(CXX) $(CXXFLAGS) tests/test_fft.cpp $(LIB) -o $@

$(BUILD)/bench_fft: bench/bench_fft.cpp $(LIB)
	$(CXX) $(CXXFLAGS) bench/bench_fft.cpp $(LIB) -o $@

test: $(BUILD)/test_fft
	./$(BUILD)/test_fft

bench: $(BUILD)/bench_fft
	./$(BUILD)/bench_fft

clean:
	rm -rf $(BUILD)
