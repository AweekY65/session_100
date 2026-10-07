#pragma once

#include <complex>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace numfft {

using Complex = std::complex<double>;
using Signal = std::vector<Complex>;

// Thrown when the input length is not a power of two (FFT) or is zero.
class InvalidLengthError : public std::invalid_argument {
public:
    explicit InvalidLengthError(const char* msg) : std::invalid_argument(msg) {}
};

// Returns true iff n is a positive power of two.
bool is_power_of_two(std::size_t n);

// Reference O(N^2) discrete Fourier transform.
//
// Convention (shared by every function in this library):
//   forward:  X[k] = sum_{n=0}^{N-1} x[n] * exp(-2*pi*i*k*n/N)   (no scaling)
//   inverse:  x[n] = (1/N) * sum_{k=0}^{N-1} X[k] * exp(+2*pi*i*k*n/N)
//
// Frequency-bin order: X[k] corresponds to digital frequency k/N cycles/sample
// (or k*fs/N Hz for sample rate fs). Bins k > N/2 are the negative frequencies
// k - N. Bin 0 is DC; for even N, bin N/2 is the Nyquist bin.
Signal dft(const Signal& x);
Signal idft(const Signal& X);

// Iterative, in-place radix-2 Cooley-Tukey FFT (decimation in time).
// Same normalization convention as dft()/idft().
// Requires x.size() to be a power of two; throws InvalidLengthError otherwise.
Signal fft(const Signal& x);
Signal ifft(const Signal& X);

// In-place variants: transform `data` directly, no copy of the input.
void fft_inplace(Signal& data);
void ifft_inplace(Signal& data);

} // namespace numfft
