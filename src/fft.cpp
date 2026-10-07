#include "fft.h"

#include <cmath>
#include <stdexcept>

namespace numfft {

namespace {
constexpr double kTwoPi = 2.0 * 3.14159265358979323846;
}

bool is_power_of_two(std::size_t n) {
    return n != 0 && (n & (n - 1)) == 0;
}

double max_abs_error(const Signal& a, const Signal& b) {
    if (a.size() != b.size()) {
        throw std::invalid_argument("max_abs_error: size mismatch");
    }
    double err = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        err = std::max(err, std::abs(a[i] - b[i]));
    }
    return err;
}

Signal dft_forward(const Signal& x) {
    const std::size_t n = x.size();
    Signal out(n);
    for (std::size_t k = 0; k < n; ++k) {
        Complex sum{0.0, 0.0};
        for (std::size_t t = 0; t < n; ++t) {
            double angle = -kTwoPi * static_cast<double>(k * t) / static_cast<double>(n);
            sum += x[t] * Complex{std::cos(angle), std::sin(angle)};
        }
        out[k] = sum;
    }
    return out;
}

Signal dft_inverse(const Signal& X) {
    const std::size_t n = X.size();
    Signal out(n);
    for (std::size_t t = 0; t < n; ++t) {
        Complex sum{0.0, 0.0};
        for (std::size_t k = 0; k < n; ++k) {
            double angle = kTwoPi * static_cast<double>(k * t) / static_cast<double>(n);
            sum += X[k] * Complex{std::cos(angle), std::sin(angle)};
        }
        out[t] = sum / static_cast<double>(n);
    }
    return out;
}

namespace {

// 迭代式 radix-2 FFT 核心：bit-reversal 重排 + 逐级蝶形，全部原位完成。
// inverse=true 时使用共轭旋转因子并在末尾除以 N。
void fft_inplace(Signal& a, bool inverse) {
    const std::size_t n = a.size();

    // bit-reversal 置换（原位交换）
    for (std::size_t i = 1, j = 0; i < n; ++i) {
        std::size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) {
            j ^= bit;
        }
        j ^= bit;
        if (i < j) {
            std::swap(a[i], a[j]);
        }
    }

    // 逐级蝶形：len = 2, 4, 8, ..., n
    for (std::size_t len = 2; len <= n; len <<= 1) {
        double angle = kTwoPi / static_cast<double>(len) * (inverse ? 1.0 : -1.0);
        Complex wlen{std::cos(angle), std::sin(angle)};
        for (std::size_t i = 0; i < n; i += len) {
            Complex w{1.0, 0.0};
            for (std::size_t j = 0; j < len / 2; ++j) {
                Complex u = a[i + j];
                Complex v = a[i + j + len / 2] * w;
                a[i + j] = u + v;
                a[i + j + len / 2] = u - v;
                w *= wlen;
            }
        }
    }

    if (inverse) {
        double inv_n = 1.0 / static_cast<double>(n);
        for (auto& v : a) {
            v *= inv_n;
        }
    }
}

void check_size(const Signal& x) {
    if (x.empty()) {
        throw std::invalid_argument("fft: empty input");
    }
    if (!is_power_of_two(x.size())) {
        throw std::invalid_argument(
            "fft: length must be a power of two (got " + std::to_string(x.size()) + ")");
    }
}

} // namespace

Signal fft_forward(const Signal& x) {
    check_size(x);
    Signal a = x; // 单次拷贝，之后全部原位计算
    fft_inplace(a, false);
    return a;
}

Signal fft_inverse(const Signal& X) {
    check_size(X);
    Signal a = X;
    fft_inplace(a, true);
    return a;
}

} // namespace numfft
