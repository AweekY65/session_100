#include "fft.h"

#include <chrono>
#include <iomanip>
#include <iostream>
#include <random>

using numfft::Complex;
using numfft::Signal;

namespace {

template <typename F>
double time_ms(F&& fn, int reps) {
    auto t0 = std::chrono::steady_clock::now();
    for (int i = 0; i < reps; ++i) fn();
    auto t1 = std::chrono::steady_clock::now();
    return std::chrono::duration<double, std::milli>(t1 - t0).count() / reps;
}

} // namespace

int main() {
    std::mt19937_64 rng(99);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);

    std::cout << "N\tFFT_ms\tDFT_ms\tFFT_ns/NlogN\tDFT_ns/N^2\n";
    std::cout << std::scientific << std::setprecision(3);

    for (std::size_t n = 256; n <= (1u << 18); n <<= 1) {
        Signal x(n);
        for (auto& v : x) v = {dist(rng), dist(rng)};

        int fft_reps = static_cast<int>(std::max<std::size_t>(3, (1u << 22) / n));
        double fft_ms = time_ms([&] { volatile Signal r = numfft::fft_forward(x); }, fft_reps);

        double fft_norm = fft_ms * 1e6 / (static_cast<double>(n) * std::log2(static_cast<double>(n)));

        std::cout << n << "\t" << fft_ms << "\t";

        if (n <= 16384) { // DFT 在大 N 下太慢，只测到 2^14
            int dft_reps = static_cast<int>(std::max<std::size_t>(1, (1u << 18) / n));
            double dft_ms = time_ms([&] { volatile Signal r = numfft::dft_forward(x); }, dft_reps);
            double dft_norm = dft_ms * 1e6 / (static_cast<double>(n) * n);
            std::cout << dft_ms << "\t" << fft_norm << "\t" << dft_norm << "\n";
        } else {
            std::cout << "-\t" << fft_norm << "\t-\n";
        }
    }
    return 0;
}
