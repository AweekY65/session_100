#include "fft.h"

#include <cmath>
#include <iostream>
#include <random>
#include <string>

using numfft::Complex;
using numfft::Signal;

namespace {

int g_failures = 0;

void report(const std::string& name, double max_err, double tol) {
    bool ok = max_err <= tol;
    if (!ok) ++g_failures;
    std::cout << (ok ? "[PASS] " : "[FAIL] ") << name
              << "  max_err = " << max_err
              << "  (tol = " << tol << ")\n";
}

void report_true(const std::string& name, bool ok) {
    if (!ok) ++g_failures;
    std::cout << (ok ? "[PASS] " : "[FAIL] ") << name << "\n";
}

constexpr double kTwoPi = 2.0 * 3.14159265358979323846;
constexpr double kTol = 1e-9;

Signal make_impulse(std::size_t n) {
    Signal x(n, {0.0, 0.0});
    x[0] = {1.0, 0.0};
    return x;
}

Signal make_constant(std::size_t n, double c) {
    return Signal(n, {c, 0.0});
}

Signal make_sine(std::size_t n, std::size_t bin) {
    Signal x(n);
    for (std::size_t t = 0; t < n; ++t) {
        x[t] = std::sin(kTwoPi * static_cast<double>(bin * t) / static_cast<double>(n));
    }
    return x;
}

Signal make_complex(std::size_t n) {
    Signal x(n);
    for (std::size_t t = 0; t < n; ++t) {
        x[t] = {std::cos(0.3 * static_cast<double>(t)),
                std::sin(0.7 * static_cast<double>(t))};
    }
    return x;
}

Signal make_random(std::size_t n, unsigned seed) {
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    Signal x(n);
    for (auto& v : x) v = {dist(rng), dist(rng)};
    return x;
}

void test_impulse() {
    // delta[n] 的 DFT 是全 1 频谱
    std::size_t n = 64;
    Signal X = numfft::fft_forward(make_impulse(n));
    double err = 0.0;
    for (const auto& v : X) err = std::max(err, std::abs(v - Complex{1.0, 0.0}));
    report("impulse: flat unit spectrum", err, kTol);
}

void test_constant() {
    // 常数 c 的 DFT 只有 bin 0 为 N*c
    std::size_t n = 128;
    double c = 2.5;
    Signal X = numfft::fft_forward(make_constant(n, c));
    double err = std::abs(X[0] - Complex{n * c, 0.0});
    for (std::size_t k = 1; k < n; ++k) err = std::max(err, std::abs(X[k]));
    report("constant: DC-only spectrum", err, kTol);
}

void test_single_frequency_sine() {
    // sin(2*pi*k0*t/N) -> X[k0] = -i*N/2, X[N-k0] = +i*N/2
    std::size_t n = 256, k0 = 7;
    Signal X = numfft::fft_forward(make_sine(n, k0));
    double err = 0.0;
    for (std::size_t k = 0; k < n; ++k) {
        Complex expected{0.0, 0.0};
        if (k == k0) expected = {0.0, -static_cast<double>(n) / 2.0};
        if (k == n - k0) expected = {0.0, static_cast<double>(n) / 2.0};
        err = std::max(err, std::abs(X[k] - expected));
    }
    report("single-frequency sine: two conjugate bins", err, kTol);
}

void test_complex_signal_vs_dft() {
    std::size_t n = 64;
    Signal x = make_complex(n);
    report("complex signal: FFT == direct DFT",
           numfft::max_abs_error(numfft::fft_forward(x), numfft::dft_forward(x)), kTol);
}

void test_random_signal_vs_dft() {
    std::size_t n = 128;
    Signal x = make_random(n, 42);
    report("random signal: FFT == direct DFT",
           numfft::max_abs_error(numfft::fft_forward(x), numfft::dft_forward(x)), kTol);
    report("random signal: inverse FFT == inverse DFT",
           numfft::max_abs_error(numfft::fft_inverse(numfft::dft_forward(x)),
                                 numfft::dft_inverse(numfft::dft_forward(x))), kTol);
}

void test_round_trip() {
    std::size_t n = 512;
    Signal x = make_random(n, 1234);
    Signal back = numfft::fft_inverse(numfft::fft_forward(x));
    report("round-trip: ifft(fft(x)) == x", numfft::max_abs_error(back, x), kTol);
    // 交叉验证：fft_inverse 与 dft_inverse 恢复结果一致
    Signal back2 = numfft::dft_inverse(numfft::fft_forward(x));
    report("round-trip: dft_inverse(fft(x)) == x", numfft::max_abs_error(back2, x), kTol);
}

void test_invalid_length() {
    bool threw = false;
    try {
        numfft::fft_forward(Signal(100, {1.0, 0.0})); // 100 不是 2 的幂
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    report_true("invalid length: N=100 rejected", threw);

    threw = false;
    try {
        numfft::fft_inverse(Signal(0));
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    report_true("invalid length: empty input rejected", threw);
}

void test_parseval() {
    // Parseval: sum|x|^2 = (1/N) * sum|X|^2
    std::size_t n = 256;
    Signal x = make_random(n, 7);
    Signal X = numfft::fft_forward(x);
    double t = 0.0, f = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        t += std::norm(x[i]);
        f += std::norm(X[i]);
    }
    report("Parseval identity", std::abs(t - f / static_cast<double>(n)), 1e-8);
}

} // namespace

int main() {
    test_impulse();
    test_constant();
    test_single_frequency_sine();
    test_complex_signal_vs_dft();
    test_random_signal_vs_dft();
    test_round_trip();
    test_invalid_length();
    test_parseval();

    std::cout << (g_failures == 0 ? "\nALL TESTS PASSED\n" : "\nSOME TESTS FAILED\n");
    return g_failures == 0 ? 0 : 1;
}
