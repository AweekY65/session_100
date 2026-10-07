# numfft — 本地 DFT / FFT 数值库

一个纯本地实现的离散傅里叶变换数值库（C++17，仅依赖标准库）。
不调用 FFTW、MKL、NumPy 或任何云计算服务；所有信号、频谱与测试
数据只存在于本地内存或本地文件中。

## 目录结构

```
src/fft.h, src/fft.cpp   核心库：dft / idft / fft / ifft
tests/test_fft.cpp       自动化测试（终端直接执行，输出最大误差）
bench/bench_fft.cpp      运行时间统计（验证复杂度趋势）
Makefile                 构建与测试入口
```

## 算法说明

### 直接 DFT（参考实现）

`dft()` 按定义直接计算，复杂度 O(N²)，作为正确性参考：

```
X[k] = sum_{n=0}^{N-1} x[n] * exp(-2*pi*i*k*n/N)
```

### Radix-2 Cooley-Tukey FFT

`fft()` 使用**迭代式**按时间抽取（DIT）radix-2 Cooley-Tukey 算法，
复杂度 O(N log N)：

1. **位反转置换（bit-reversal permutation）**：原地交换，把输入重排为
   蝶形网络所需的顺序，不分配额外数组。
2. **蝶形运算（butterfly）**：对 len = 2, 4, 8, ..., N 共 log2(N) 级，
   每级把两个半段合并。核心蝶形公式为：

   ```
   A = u + w * v
   B = u - w * v
   ```
   其中 `w = exp(-2*pi*i*j/len)` 是旋转因子（twiddle factor）。
   旋转因子在每级内通过递推 `w *= w_step` 更新，避免重复三角函数调用。

整个变换**原地（in-place）**进行，无递归、无逐层复制，避免了递归
实现中常见的大量子数组拷贝。

### 逆变换

`ifft()` 通过共轭技巧复用同一蝶形网络：

```
ifft(X) = conj( fft( conj(X) ) ) / N
```

保证正/逆变换使用完全一致的约定，不会出现互相矛盾的缩放。

## 归一化与频率索引约定

- **正变换不缩放**：`X[k] = sum_n x[n] * exp(-2*pi*i*k*n/N)`
- **逆变换缩放 1/N**：`x[n] = (1/N) * sum_k X[k] * exp(+2*pi*i*k*n/N)`
- 因此 `ifft(fft(x)) == x`（数值容差内）。

**频率 bin 顺序**：输出 `X[k]`（k = 0 .. N-1）对应数字频率 `k/N`
周期/采样（采样率 fs 时为 `k*fs/N` Hz）：

- `k = 0`：直流（DC）分量；
- `k = 1 .. N/2 - 1`：正频率；
- `k = N/2`（N 为偶数）：奈奎斯特频率；
- `k = N/2 + 1 .. N-1`：负频率，即 `X[k]` 对应频率 `k - N`。

这是与 FFTW / NumPy 相同的标准排列（未做 fftshift）。

## 非 2 的幂长度

`fft()` / `ifft()` **明确拒绝**长度不是 2 的幂（或为 0）的输入，
抛出 `numfft::InvalidLengthError`（继承自 `std::invalid_argument`）。
该行为稳定且确定：不做隐式 padding，调用方如需 padding 应自行
补零到下一个 2 的幂。`dft()` / `idft()` 作为参考实现支持任意
非零长度。

## 构建与测试

```sh
make          # 编译库、测试与基准程序（输出到 build/）
make test     # 运行测试套件，逐项打印 max_err
make bench    # 运行基准，打印 DFT/FFT 耗时及归一化指标
make clean
```

### 测试覆盖

`make test` 直接在终端执行，每项测试输出最大绝对误差：

- 单位冲激（impulse）：频谱恒为 1
- 常数信号：仅 DC bin 为 N*c，其余为 0
- 单频正弦：峰值出现在 k0 与 N-k0
- 复指数信号：单 bin 响应 X[k0] = N*A
- 随机复数信号：FFT 与直接 DFT 对比（多种长度）
- round-trip：`ifft(fft(x)) == x`、`idft(dft(x)) == x`
- 非法长度（0、3、100、1025）：必须抛出 `InvalidLengthError`
- 多种长度（1..512）下 FFT 与 DFT 一致性

任一测试失败时进程以非零状态码退出。

### 复杂度验证

`make bench` 对 N = 64 .. 2^20 计时，并打印 `DFT/N^2` 与
`FFT/(N*log2 N)`。两者随 N 增长近似保持恒定，即分别验证了
O(N²) 与 O(N log N) 的复杂度趋势（DFT 在 N > 2^14 时跳过一次
全量计时以控制总耗时）。

## 使用示例

```cpp
#include "src/fft.h"

numfft::Signal x = {{1,0}, {2,0}, {3,0}, {4,0}};
auto X = numfft::fft(x);        // 正变换
auto y = numfft::ifft(X);       // 逆变换，y ≈ x
// 原地版本：
numfft::fft_inplace(x);
numfft::ifft_inplace(x);
```
