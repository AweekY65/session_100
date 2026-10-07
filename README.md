# numfft — 本地 DFT / FFT 数值库

纯本地 C++17 实现的离散傅里叶变换库，不依赖 FFTW、MKL 或任何外部 FFT 实现；
所有信号与频谱仅存在于内存或本地文件中。

## 功能

- `dft_forward` / `dft_inverse`：直接 O(N²) DFT，作为参考实现，支持任意长度。
- `fft_forward` / `fft_inverse`：迭代式 radix-2 Cooley-Tukey FFT，复数输入，
  长度为 2 的幂时与直接 DFT 在 1e-9 容差内一致。

## 蝶形运算

FFT 采用迭代（非递归）实现，避免递归中的大量数组复制：

1. **bit-reversal 置换**：对输入索引做二进制位反转并原位交换，将递归分治的
   重排步骤一次性完成。
2. **逐级蝶形**：子变换长度 `len` 从 2 开始逐级翻倍到 N。每级使用旋转因子
   `wlen = exp(±2πi/len)`，对每对元素执行蝶形：
   ```
   u = a[j];  v = a[j + len/2] * w;
   a[j]          = u + v;
   a[j + len/2]  = u - v;
   ```
   旋转因子 `w` 在每级内累乘更新（`w *= wlen`），全部计算原位（in-place）完成，
   除输出缓冲的一次拷贝外无额外内存分配。

## 频率索引约定

forward 输出 `X[k]`（`k = 0..N-1`）按标准 DFT 顺序排列：

| k 的范围        | 含义     | 物理频率（采样间隔 dt） |
| --------------- | -------- | ----------------------- |
| `k = 0`         | 直流分量 | 0                       |
| `1 <= k <= N/2` | 正频率   | `f_k = k / (N·dt)`      |
| `N/2 < k < N`   | 负频率   | `f_k = (k - N) / (N·dt)`|

## 归一化规则

- **forward**：`X[k] = Σ_n x[n]·exp(-2πi·k·n/N)`（不缩放）
- **inverse**：`x[n] = (1/N)·Σ_k X[k]·exp(+2πi·k·n/N)`（1/N 缩放全部放在逆变换）

因此 `fft_inverse(fft_forward(x)) == x`（数值容差内），两个方向缩放规则一致、无矛盾。

## 非 2 的幂输入

`fft_forward` / `fft_inverse` **明确拒绝**长度不是 2 的幂的输入（以及空输入），
抛出 `std::invalid_argument`，行为稳定可预期。需要处理任意长度时请使用
`dft_forward` / `dft_inverse`，或在调用 FFT 前自行补零（padding）。

## 构建与测试

```sh
make test    # 编译并运行测试（终端直接输出每项最大误差）
make bench   # 运行基准，输出 FFT/DFT 耗时及归一化复杂度指标
```

或使用 CMake：

```sh
cmake -B build && cmake --build build
cd build && ctest
```

测试覆盖：impulse（平坦单位频谱）、constant（仅直流）、单频正弦（共轭双谱线）、
复数信号与随机信号（FFT 对比直接 DFT）、round-trip（`ifft(fft(x)) == x`）、
非法长度拒绝、Parseval 恒等式。每项测试输出 `max_err`，全部在容差内通过。

基准输出中 `FFT_ns/NlogN` 与 `DFT_ns/N^2` 两列随 N 增大保持恒定，
分别验证 O(N log N) 与 O(N²) 的复杂度趋势。实测（N = 256 … 262144）：

```
N       FFT_ms      DFT_ms
256     4.5e-03     1.6e+00
4096    9.8e-02     4.0e+02
16384   4.7e-01     6.4e+03
262144  9.4e+00     -
```

## 目录结构

```
include/fft.h      公共接口与约定说明
src/fft.cpp        DFT / FFT / IFFT 实现
tests/test_fft.cpp 自动化测试（输出最大误差）
bench/bench.cpp    运行时间统计
```
