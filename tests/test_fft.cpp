// Automated tests for the local DFT/FFT library.
// Every test prints its maximum absolute error; exit code is non-zero on failure.
#include "../src/fft.h"

#include <cmath>
#include <cstdio>
#include <random>

namespace {

using numfft::Complex;
using numfft::Signal;

constexpr double kPi = 3.14159265358979323846264338327950288;
constexpr double kTol = 1e-9;

int g_failures = 0;

double max_abs_diff(const Signal& a, const Signal& b) {
    double m = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i)
        m = std::max(m, std::abs(a[i] - b[i]));
    return m;
}

void report(const char* name, double max_err, double tol) {
    const bool ok = max_err <= tol;
    std::printf("[%-4s] %-38s max_err = %.3e (tol %.1e)\n",
                ok ? "PASS" : "FAIL", name, max_err, tol);
    if (!ok) ++g_failures;
}

void report_bool(const char* name, bool ok) {
    std::printf("[%-4s] %-38s\n", ok ? "PASS" : "FAIL", name);
    if (!ok) ++g_failures;
}

// 1. Impulse: x = [1, 0, 0, ...] -> X[k] = 1 for all k.
void test_impulse() {
    const std::size_t n = 64;
    Signal x(n, {0.0, 0.0});
    x[0] = 1.0;
    Signal X = numfft::fft(x);
    double err = 0.0;
    for (const auto& v : X) err = std::max(err, std::abs(v - Complex{1.0, 0.0}));
    report("impulse: fft([1,0,...]) == 1", err, kTol);
}

// 2. Constant: x[n] = c -> X[0] = N*c, X[k!=0] = 0.
void test_constant() {
    const std::size_t n = 128;
    const Complex c{0.7, -0.3};
    Signal x(n, c);
    Signal X = numfft::fft(x);
    double err = std::abs(X[0] - c * static_cast<double>(n));
    for (std::size_t k = 1; k < n; ++k) err = std::max(err, std::abs(X[k]));
    report("constant: DC = N*c, others = 0", err, kTol * n);
}

// 3. Single-frequency sine at bin k0 -> peaks at k0 and N-k0 with amp N/2.
void test_sine() {
    const std::size_t n = 256;
    const std::size_t k0 = 5;
    Signal x(n);
    for (std::size_t t = 0; t < n; ++t)
        x[t] = std::sin(2.0 * kPi * static_cast<double>(k0 * t) / n);
    Signal X = numfft::fft(x);
    double err = 0.0;
    for (std::size_t k = 0; k < n; ++k) {
        Complex expected{0.0, 0.0};
        if (k == k0) expected = {0.0, -static_cast<double>(n) / 2.0};
        if (k == n - k0) expected = {0.0, static_cast<double>(n) / 2.0};
        err = std::max(err, std::abs(X[k] - expected));
    }
    report("sine: peaks at k0 and N-k0", err, kTol * n);
}

// 4. Complex exponential at bin k0 -> single non-zero bin X[k0] = N*A.
void test_complex_exponential() {
    const std::size_t n = 64;
    const std::size_t k0 = 11;
    const Complex amp{0.5, 0.25};
    Signal x(n);
    for (std::size_t t = 0; t < n; ++t) {
        const double angle = 2.0 * kPi * static_cast<double>(k0 * t) / n;
        x[t] = amp * Complex{std::cos(angle), std::sin(angle)};
    }
    Signal X = numfft::fft(x);
    double err = 0.0;
    for (std::size_t k = 0; k < n; ++k) {
        Complex expected = (k == k0) ? amp * static_cast<double>(n)
                                     : Complex{0.0, 0.0};
        err = std::max(err, std::abs(X[k] - expected));
    }
    report("complex exp: single bin X[k0] = N*A", err, kTol * n);
}

// 5. Random complex signal: FFT must match direct DFT.
void test_random_vs_dft() {
    std::mt19937_64 rng(20261004);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    double worst = 0.0;
    for (std::size_t n : {2UL, 4UL, 8UL, 16UL, 64UL, 256UL, 1024UL}) {
        Signal x(n);
        for (auto& v : x) v = {dist(rng), dist(rng)};
        worst = std::max(worst, max_abs_diff(numfft::fft(x), numfft::dft(x)));
    }
    report("random complex: fft == dft", worst, 1e-8);
}

// 6. Round-trip: ifft(fft(x)) == x, also idft/dft pair.
void test_round_trip() {
    std::mt19937_64 rng(7);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    const std::size_t n = 512;
    Signal x(n);
    for (auto& v : x) v = {dist(rng), dist(rng)};
    const double err_fft = max_abs_diff(numfft::ifft(numfft::fft(x)), x);
    const double err_dft = max_abs_diff(numfft::idft(numfft::dft(x)), x);
    report("round-trip: ifft(fft(x)) == x", err_fft, 1e-10);
    report("round-trip: idft(dft(x)) == x", err_dft, 1e-10);
}

// 7. Illegal lengths (0, 3, 100, 2^10+1) must throw InvalidLengthError.
void test_invalid_length() {
    bool ok = true;
    for (std::size_t n : {0UL, 1UL * 3, 100UL, 1025UL}) {
        Signal x(n, {0.0, 0.0});
        try {
            numfft::fft(x);
            ok = false; // should have thrown
        } catch (const numfft::InvalidLengthError&) {
        }
        try {
            numfft::ifft(x);
            ok = false;
        } catch (const numfft::InvalidLengthError&) {
        }
    }
    report_bool("invalid length: throws InvalidLengthError", ok);
}

// 8. FFT vs direct DFT on structured signals of several sizes.
void test_fft_matches_dft_sizes() {
    double worst = 0.0;
    for (std::size_t n : {1UL, 2UL, 4UL, 32UL, 128UL, 512UL}) {
        Signal x(n);
        for (std::size_t t = 0; t < n; ++t)
            x[t] = {std::cos(0.3 * t) + 0.1 * t, std::sin(0.7 * t)};
        worst = std::max(worst, max_abs_diff(numfft::fft(x), numfft::dft(x)));
    }
    report("fft == dft across sizes 1..512", worst, 1e-8);
}

} // namespace

int main() {
    std::printf("=== numfft test suite ===\n");
    test_impulse();
    test_constant();
    test_sine();
    test_complex_exponential();
    test_random_vs_dft();
    test_round_trip();
    test_invalid_length();
    test_fft_matches_dft_sizes();
    std::printf("=== %s (%d failure%s) ===\n",
                g_failures == 0 ? "ALL TESTS PASSED" : "TESTS FAILED",
                g_failures, g_failures == 1 ? "" : "s");
    return g_failures == 0 ? 0 : 1;
}
