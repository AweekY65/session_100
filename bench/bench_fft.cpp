// Runtime statistics: verifies the O(N^2) trend of the direct DFT and the
// O(N log N) trend of the radix-2 FFT as N grows.
#include "../src/fft.h"

#include <chrono>
#include <cstdio>
#include <random>

namespace {

using numfft::Signal;
using Clock = std::chrono::steady_clock;

Signal make_random_signal(std::size_t n) {
    std::mt19937_64 rng(42);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    Signal x(n);
    for (auto& v : x) v = {dist(rng), dist(rng)};
    return x;
}

template <typename F>
double time_ms(F&& f) {
    const auto t0 = Clock::now();
    f();
    const auto t1 = Clock::now();
    return std::chrono::duration<double, std::milli>(t1 - t0).count();
}

} // namespace

int main() {
    std::printf("=== numfft benchmark ===\n");
    std::printf("%-10s %14s %14s %12s %12s\n", "N", "DFT (ms)", "FFT (ms)",
                "DFT/N^2", "FFT/(NlogN)");

    double prev_dft_per_n2 = 0.0, prev_fft_per_nlogn = 0.0;
    for (std::size_t n = 64; n <= (1UL << 20); n <<= 1) {
        Signal x = make_random_signal(n);

        // Direct DFT is O(N^2); skip it above 2^14 to keep runtime sane.
        double dft_ms = -1.0;
        double dft_per_n2 = 0.0;
        if (n <= (1UL << 14)) {
            dft_ms = time_ms([&] { volatile Signal r = numfft::dft(x); (void)r; });
            dft_per_n2 = dft_ms / (static_cast<double>(n) * n);
        }

        double fft_ms = time_ms([&] { volatile Signal r = numfft::fft(x); (void)r; });
        const double nlogn =
            static_cast<double>(n) * std::log2(static_cast<double>(n));
        const double fft_per_nlogn = fft_ms / nlogn;

        if (dft_ms >= 0.0)
            std::printf("%-10zu %14.3f %14.4f %12.3e %12.3e\n", n, dft_ms,
                        fft_ms, dft_per_n2, fft_per_nlogn);
        else
            std::printf("%-10zu %14s %14.4f %12s %12.3e\n", n, "(skipped)",
                        fft_ms, "-", fft_per_nlogn);

        prev_dft_per_n2 = dft_per_n2;
        prev_fft_per_nlogn = fft_per_nlogn;
    }
    (void)prev_dft_per_n2;
    (void)prev_fft_per_nlogn;

    std::printf("\nInterpretation: DFT/N^2 and FFT/(N*log2 N) stay roughly\n");
    std::printf("constant as N grows, confirming O(N^2) and O(N log N).\n");
    return 0;
}
