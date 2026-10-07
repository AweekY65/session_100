#include "fft.h"

#include <cmath>

namespace numfft {

namespace {
constexpr double kPi = 3.14159265358979323846264338327950288;
constexpr double kTwoPi = 2.0 * kPi;
} // namespace

bool is_power_of_two(std::size_t n) {
    return n != 0 && (n & (n - 1)) == 0;
}

Signal dft(const Signal& x) {
    const std::size_t n = x.size();
    if (n == 0) throw InvalidLengthError("dft: empty input");
    Signal out(n);
    for (std::size_t k = 0; k < n; ++k) {
        Complex sum{0.0, 0.0};
        for (std::size_t t = 0; t < n; ++t) {
            const double angle = -kTwoPi * static_cast<double>(k * t % n) /
                                 static_cast<double>(n);
            sum += x[t] * Complex{std::cos(angle), std::sin(angle)};
        }
        out[k] = sum;
    }
    return out;
}

Signal idft(const Signal& X) {
    const std::size_t n = X.size();
    if (n == 0) throw InvalidLengthError("idft: empty input");
    Signal out(n);
    for (std::size_t t = 0; t < n; ++t) {
        Complex sum{0.0, 0.0};
        for (std::size_t k = 0; k < n; ++k) {
            const double angle = +kTwoPi * static_cast<double>(k * t % n) /
                                 static_cast<double>(n);
            sum += X[k] * Complex{std::cos(angle), std::sin(angle)};
        }
        out[t] = sum / static_cast<double>(n);
    }
    return out;
}

void fft_inplace(Signal& data) {
    const std::size_t n = data.size();
    if (!is_power_of_two(n)) {
        throw InvalidLengthError("fft: length must be a power of two");
    }
    if (n == 1) return;

    // Bit-reversal permutation, done in place with swaps (no extra buffer).
    for (std::size_t i = 1, j = 0; i < n; ++i) {
        std::size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(data[i], data[j]);
    }

    // Iterative butterfly stages: len = 2, 4, 8, ..., n.
    // Twiddle factors are computed by recurrence per stage to avoid
    // repeated trig calls; only O(log n) stages, no recursion, no copies.
    for (std::size_t len = 2; len <= n; len <<= 1) {
        const std::size_t half = len >> 1;
        const double stage_angle = -kTwoPi / static_cast<double>(len);
        const Complex w_step{std::cos(stage_angle), std::sin(stage_angle)};
        for (std::size_t base = 0; base < n; base += len) {
            Complex w{1.0, 0.0};
            for (std::size_t j = 0; j < half; ++j) {
                const Complex u = data[base + j];
                const Complex v = data[base + j + half] * w;
                data[base + j] = u + v;
                data[base + j + half] = u - v;
                w *= w_step;
            }
        }
    }
}

void ifft_inplace(Signal& data) {
    const std::size_t n = data.size();
    if (!is_power_of_two(n)) {
        throw InvalidLengthError("ifft: length must be a power of two");
    }
    // Inverse via conjugation trick: ifft(X) = conj(fft(conj(X))) / N.
    // This reuses the exact same butterfly network, guaranteeing that
    // forward and inverse scaling are consistent (1/N applied only here).
    for (auto& v : data) v = std::conj(v);
    fft_inplace(data);
    const double inv_n = 1.0 / static_cast<double>(n);
    for (auto& v : data) v = std::conj(v) * inv_n;
}

Signal fft(const Signal& x) {
    Signal out = x;
    fft_inplace(out);
    return out;
}

Signal ifft(const Signal& X) {
    Signal out = X;
    ifft_inplace(out);
    return out;
}

} // namespace numfft
