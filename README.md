# cpp-ap-filter-tools
- Created by and copywrited by me
- This code is made availiable under the GPLv3 license here: www.gnu.org/licenses/gpl-3.0.en.html
## Warning, not fully tested, use at your own risk
- These functions seem correct but have not been sufficiently tested
## List of functions
- lp2lp - lowpass to lowpass
- lp2hp - lowpass to highpass
- lp2bp - lowpass to bandpass
- bt_freq_warp - warp a frequency for bilinear transforms
- bilin - bilinear transform
## List of classes
- filtab
	- container class with two vectors a and b
- FastBilin
	- a class that caches results from previous bilinear transforms to make repeat transforms faster
- FastBPF
	- a class that caches results from previous lp2bp transforms to make repeat transforms faster
## General Notes
### analog filters are in the form
		b[0]s^N + b[0]s^(N - 1) + ... + b[N - 1]s + b[N]
	H(s) = -------------------------------------------------
		a[0]s^N + a[0]s^(N - 1) + ... + a[N - 1]s + a[N]
### digital filters are in the form
		b[0] + b[1]z^-1 + ... + b[N - 1]z^-(N - 1) + b[N]z^-N
	H(z) = ------------------------------------------------------
		a[0] + a[1]z^-1 + ... + a[N - 1]z^-(N - 1) + a[N]z^-N
### The a and b vectors in the filtab class are expected to be the same size
- a and b are assumed to be of the same order, with the unneeded coefficients being set to zero
### This header file has no functions to generate analog prototypes
- You have to provide the prototypes yourself
- I want to implement the creation of butterworth, chebyshev, and elliptical prototypes at some point
### Be mindful of template type inference, especially for bt_freq_warp
- It's best to explicitly declare the type of the template
- T is only intended to be either double or float
### The functions lp2lp, lp2hp, lp2bp, and bilin have two versions
- One in which the user is expected to provide the filtab of the prototype and the filtab of the transform
- If the filtab of the transform is correctly sized then no memory will be allocated (lp2bp and bilin excluded unless using the class versions) ### making the transform very fast
- Another in which the user only provides the prototype and the function will create the filtab of the transform
	- Because memory is allocated this version is slower than the previous
### The functions lp2hp, lp2bp, and bilin do not support inplace transformations, however lp2lp does
### lp2lp works as a general frequency scaler
### FastBilin and FastBPF are fast
- The bilin and lp2bpf located in these classes are about 9-10x faster than the ones located outside of their classes for repeat transforms of the same N
- These classes are fast because they cache results of previous transforms
- By default, the cache will erase once more than 10 transforms are stored
	- This is to prevent memory leaks, but there is certainly a better way to do this
### The poly_exp function is very bad
- There are probably some very fast algorithms to calculate a polynomial raised to a power, but with the FastBilin and FastBPF classes there is little motivation to implement a faster algorithm
### These functions have not been properly tested yet, but I think they might work
- Use at your own risk
