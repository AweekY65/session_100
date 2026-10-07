#pragma once

#include <complex>
#include <vector>
#include <cstddef>

namespace numfft {

using Complex = std::complex<double>;
using Signal = std::vector<Complex>;

// 频率 bin 顺序约定（forward 输出 X[k], k = 0..N-1）：
//   k = 0          : 直流分量
//   1 <= k <= N/2  : 正频率 f_k = k / (N * dt)
//   N/2 < k <= N-1 : 负频率 f_k = (k - N) / (N * dt)
// 归一化约定：
//   forward: X[k] = sum_n x[n] * exp(-2*pi*i*k*n/N)   （不缩放）
//   inverse: x[n] = (1/N) * sum_k X[k] * exp(+2*pi*i*k*n/N) （1/N 缩放）
// 因此 inverse(forward(x)) == x。

// 直接 O(N^2) DFT，作为参考实现。任意长度均可。
Signal dft_forward(const Signal& x);
Signal dft_inverse(const Signal& X);

// 迭代 radix-2 Cooley-Tukey FFT（原位计算，无递归复制）。
// 输入长度必须是 2 的幂，否则抛出 std::invalid_argument。
Signal fft_forward(const Signal& x);
Signal fft_inverse(const Signal& X);

// 工具
bool is_power_of_two(std::size_t n);
double max_abs_error(const Signal& a, const Signal& b);

} // namespace numfft
