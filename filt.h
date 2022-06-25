/* 
Created by and copywrited by Dante Viloria - dviloria1919@gmail.com
This code is made availiable under the GPLv3 license here: www.gnu.org/licenses/gpl-3.0.en.html

General Notes:
analog filters are in the form
           b[0]s^N + b[0]s^(N - 1) + ... + b[N - 1]s + b[N]
	H(s) = -------------------------------------------------
           a[0]s^N + a[0]s^(N - 1) + ... + a[N - 1]s + a[N]

digital filters are in the form
		   b[0] + b[1]z^-1 + ... + b[N - 1]z^-(N - 1) + b[N]z^-N
	H(z) = -----------------------------------------------------
		   a[0] + a[1]z^-1 + ... + a[N - 1]z^-(N - 1) + a[N]z^-N

The a and b vectors in the filtab class are expected to be the same size
	a and b are assumed to be of the same order, with the unneeded coefficients being set to zero

This header file has no functions to generate analog prototypes
	You have to provide the prototypes yourself
	I want to implement the creation of butterworth, chebyshev, and elliptical prototypes at some point

Be mindful of template type inference, especially for the frequency warping function
	It's best to explicitly declare the type of the template
	T is only intended to be either double or float

The functions lp2lp, lp2hp, lp2bp, and bilin have two versions
	One in which the user is expected to provide the filtab of the prototype and the filtab of the transform
		If the filtab of the transform is correctly sized then no memory will be allocated (lp2bp and bilin excluded unless using the class versions) making the transform very fast
	Another in which the user only provides the prototype and the function will create the filtab of the transform
		Because memory is allocated this version is slower than the previous

The functions lp2hp, lp2bp, and bilin do not support inplace transformations, however lp2lp does

lp2lp works as a general frequency scaler

FastBilin and FastBPF are fast
	The bilin and lp2bpf located in these classes are about 9-10x faster than the ones located outside of their classes for repeat transforms of the same N
	These classes are fast because they cache results of previous transforms
	By default, the cache will erase once more than 10 transforms are stored
		This is to prevent memory leaks, but there is certainly a better way to do this

The poly_exp function is very bad
	There are probably some very fast algorithms to calculate a polynomial raised to a power, 
	but with the FastBilin and FastBPF classes there is little motivation to implement a faster algorithm

These functions have not been properly tested yet, but I think they might work
	Use at your own risk
*/
#include <iostream>
#include <vector>
#include <unordered_map>
#include <ctgmath>
namespace filt {
	// multiply two polynomials
	template <typename T>
	std::vector<T> poly_mult(const std::vector<T>& poly1, const std::vector<T>& poly2) {
		int a_max = poly2.size();
		int b_max = poly1.size();
		std::vector<T> c(a_max + b_max - 1);
		for (auto i = 0; i < a_max; i++) { 
			for (auto j = 0; j < b_max; j++) {
				c[i + j] += poly2[i] * poly1[j];
			} 
		}
		return c;
	}
	/* 
	only for integer powers >= 0
	this function is very inefficient and has a lot of room for improvement
	*/
	template <typename T>
	std::vector<T> poly_exp(const std::vector<T>& poly, int pow) {
		int a_max = poly.size();
		std::vector<T> poly_raised;
		if (pow <= 0) {
			poly_raised = { 1 };
		}
		else {
			poly_raised = poly;
			for (auto i = 0; i < pow - 1; i++) {
				poly_raised = poly_mult(poly_raised, poly);
			}
		}
		return poly_raised;
	}
	/*
	if a vector vec = {1} and length = 5 is given
	return {0, 0, 0, 0, 1}
	*/
	template <typename T>
	std::vector<T> pad_vec(const std::vector<T>& vec, int length) {
		if (vec.size() == length) {
			return vec;
		}
		std::vector<T> vec_padded(length);
		int offset = length - vec.size();
		for (auto i = offset; i < length; i++) {
			vec_padded[i] = vec[i - offset];
		}
		return vec_padded;
	}
	/*
	container for a standard filter in a b form
	*/
	template <typename T>
	class filtab {
	public:
		std::vector<T> a;
		std::vector<T> b;
		filtab() {}
		filtab(std::vector<T> b, std::vector<T> a) {
			this->b = b;
			this->a = a;
		}
		filtab(int N) {
			this->resize(N);
		}
		const int getN() const {
			return this->a.size() - 1;
		}
		/*
		resize the a and b vectors to accomodate a filter of order N
		return 1 if the vectors were resized
		return 0 if the vectors were not resized
		*/
		int resize(int N) {
			if (this->getN() == N) {
				return 0;
			}
			else {
				this->a.resize(N + 1);
				this->b.resize(N + 1);
				for (auto i = 0; i < N + 1; i++) {
					this->a[i] = 0;
					this->b[i] = 0;
				}
				return 1;
			}
		}
		/*
		divide all coefficients by a[0]
		*/
		filtab<T>& normalize() {
			T k = this->a[0];
			for (auto i = 0; i < this->a.size(); i++) {
				this->a[i] /= k;
				this->b[i] /= k;
			}
			return *this;
		}
		void print() {
			using namespace std;
			cout << "b: ";
			for (auto i = 0; i < b.size(); i++) {
				cout << b[i];
				if (i < b.size() - 1) {
					cout << ", ";
				}
			}
			cout << endl << "a: ";
			for (auto i = 0; i < a.size(); i++) {
				cout << a[i];
				if (i < a.size() - 1) {
					cout << ", ";
				}
			}

			cout << endl << endl;
		}
	};
	/* 
	transforms an analog lpf with a cutoff of 1 rad/s to a lpf with a cutoff of wo rad/s
	*/
	template <typename T>
	filtab<T> lp2lp(const filtab<T>& lpap, T wo) {
		filtab<T> lp;
		return lp2lp(lpap, lp, wo);
	}
	template <typename T>
	filtab<T>& lp2lp(const filtab<T>& lpap, filtab<T>& lp, T wo) {
		int N = lpap.getN();
		lp.resize(N);
		T k = 0;
		for (auto i = 0; i < N + 1; i++) {
			k = pow(wo, i);
			lp.a[i] = k * lpap.a[i];
			lp.b[i] = k * lpap.b[i];
		}
		lp.normalize();
		return lp;
	}
	/* 
	transforms an analog lpf with a cutoff of 1 rad/s to a hpf with a cutoff of wo rad/s
	*/
	template <typename T>
	filtab<T> lp2hp(const filtab<T>& lpap, T wo) {
		filtab<T> hp;
		return lp2hp(lpap, hp, wo);
	}
	template <typename T>
	filtab<T>& lp2hp(const filtab<T>& lpap, filtab<T>& hp, T wo) {
		int N = lpap.getN();
		hp.resize(N);
		T k = 0;
		for (auto i = 0; i < N + 1; i++) {
			k = pow(wo, -i);
			hp.a[N - i] = k * lpap.a[i];
			hp.b[N - i] = k * lpap.b[i];
		}
		hp.normalize();
		return hp;
	}
	/* 
	transforms an analog lpf with a cutoff of 1 rad/s to a bpf with a geometric center of wo rad/s and quality factor Q
	*/
	template <typename T>
	filtab<T> lp2bp(const filtab<T>& lpap, T wo, T Q) {
		filtab<T> bp;
		return lp2bp(lpap, bp, wo, Q);
	}
	template <typename T>
	filtab<T>& lp2bp(const filtab<T>& lpap, filtab<T>& bp, T wo, T Q) {
		int N = lpap.getN();
		bp.resize(2 * N);
		std::vector<T> l1{ (T)wo, 0 }, l2{ (T)Q, 0, (T)(Q * wo * wo) }, part;
		for (auto i = 0; i < N + 1; i++) {
			part = pad_vec(poly_mult(poly_exp(l1, i), poly_exp(l2, N - i)), 2 * N + 1);
			for (auto j = 0; j < 2 * N + 1; j++) {
				bp.a[j] += lpap.a[i] * part[j];
				bp.b[j] += lpap.b[i] * part[j];
			}
		}
		bp.normalize();
		return bp;
	}
	/* 
	Warps an analog filter's frequency such that an analog filter with cutoff frequency bt_freq_warp(wo, fs) will produce a digital filter with cutoff wo when a bilinear transform is performed with a sampling frequency of fs.
	*/
	template <typename T>
	T bt_freq_warp(T wo, T fs) {
		const T pi = 3.14159265358979323846;
		return 2 * fs * tan(pi * wo / fs);
	}
	/*
	converts an analog filter to a digital filter
	*/
	template <typename T>
	filtab<T> bilin(const filtab<T>& analog, T fs) {
		filtab<T> digital;
		return bilin(analog, digital, fs);
	}
	template <typename T>
	filtab<T>& bilin(const filtab<T>& analog, filtab<T>& digital, T fs) {
		int N = analog.getN();
		digital.resize(N);
		T k1 = 2 * fs, k2 = 0;
		std::vector<T> l1{ 1, 1 }, l2{ 1, -1 }, part;
		for (auto i = 0; i < N + 1; i++) {
			part = poly_mult(poly_exp(l1, i), poly_exp(l2, N - i));
			k2 = pow(k1, -i);
			for (auto j = 0; j < N + 1; j++) {
				digital.a[j] += k2 * analog.a[i] * part[j];
				digital.b[j] += k2 * analog.b[i] * part[j];
			}
		}
		digital.normalize();
		return digital;
	}
	
	/* 
	A class that has a member function that computes bilinear transforms. 
	This class caches data from previous transforms to make repeated transforms of the same order more efficient. 
	A max number of cached transforms can be set and if that number is reached the cache is erased. 
	A better system should be implemented but this should prevent memory leaks. 
	*/
	template <typename T>
	class FastBilin {
	private:
		int no_of_cached_transforms = 0;
		const int max_no_of_cached_transforms = 10;
		std::vector<T> l1{ 1, 1 };
		std::vector<T> l2{ 1, -1 };
		std::unordered_map<int, std::vector<std::vector<T>>> cached_transforms;

		/* 
		Compute a cached transform if necessary or possible
		return 0 if a cached transform exists for N
		return 1 if a cached transform was generated
		return -1 if a transform does not exist and cannot be generated
		*/
		int generate_transform(int N) {
			if (this->cached_transforms.count(N) > 0) {
				return 0;
			}
			else if (N <= 0) {
				return -1;
			}
			else {
				using namespace std;
				if (this->no_of_cached_transforms == this->max_no_of_cached_transforms) {
					this->cached_transforms.clear();
					this->no_of_cached_transforms = 0;
				}
				vector<vector<T>>& transform = this->cached_transforms[N];
				transform.resize(N + 1);
				for (auto i = 0; i < N + 1; i++) {
					transform[i] = poly_mult(poly_exp(this->l1, i), poly_exp(this->l2, N - i));
				}
				this->no_of_cached_transforms++;
				return 1;
			}
		}
	public:
		FastBilin() {
		}
		FastBilin(std::vector<int> list_of_Ns) {
			generate_transforms(list_of_Ns);
		}
		int get_no_cached_transforms() {
			return no_of_cached_transforms;
		}
		void generate_transforms(std::vector<int> list_of_Ns) {
			for (auto i = 0; i < list_of_Ns.size(); i++) {
				this->generate_transform(list_of_Ns[i]);
			}
		}
		/*
		convert an analog filter to a digital filter
		*/
		filtab<T> bilin(const filtab<T>& analog, T fs) {
			filtab<T> digital;
			return this->bilin(analog, digital, fs);
		}
		filtab<T>& bilin(const filtab<T>& analog, filtab<T>& digital, T fs) {
			using namespace std;
			int N = analog.getN();
			digital.resize(N);
			this->generate_transform(N);
			vector<vector<T>>& transform = cached_transforms[N];
			T k1 = 2 * fs, k2 = 0;
			for (auto i = 0; i < N + 1; i++) {
				k2 = pow(k1, -i);
				for (auto j = 0; j < N + 1; j++) {
					digital.a[j] += k2 * analog.a[i] * transform[i][j];
					digital.b[j] += k2 * analog.b[i] * transform[i][j];
				}
			}
			digital.normalize();
			return digital;
		}
	};
	/* 
	A class that has a member function transforms lowpass analog prototypes with wo of 1 rad / s to bandpass filters
	This class caches data from previous transforms to make repeated transforms of the same order more efficient. 
	The transform cache is based on the order of the analog prototype, not the bandpass filter.
	This means that a cached transform for N = 2 results in a BPF of N = 4;
	A max number of cached transforms can be set and if that number is reached the cache is erased. 
	A better system should be implemented but this should prevent memory leaks.
	*/
	template <typename T>
	class FastBPF {
	private:
		std::vector<T> l1{ 1, 0 };
		std::vector<T> l2{ 1, 0, 1 };
		std::unordered_map<int, std::vector<std::vector<T>>> cached_transforms;
		int no_of_cached_transforms = 0;
		const int max_no_of_cached_transforms = 10;

		/* 
		Compute a cached transform if necessary or possible.
		return 0 if a cached transform exists for N.
		return 1 if a cached transform was generated.
		return -1 if a transform does not exist and cannot be generated .
		*/
		int generate_transform(int N) {
			if (this->cached_transforms.count(N) > 0) {
				return 0;
			}
			else if (N <= 0) {
				return -1;
			}
			else {
				using namespace std;
				if (this->no_of_cached_transforms == this->max_no_of_cached_transforms) {
					this->cached_transforms.clear();
					this->no_of_cached_transforms = 0;
				}
				vector<vector<T>>& transform = this->cached_transforms[N];
				transform.resize(N + 1);
				for (auto i = 0; i < N + 1; i++) {
					transform[i] = pad_vec(poly_mult(poly_exp(this->l1, i), poly_exp(this->l2, N - i)), 2 * N + 1);
				}
				this->no_of_cached_transforms++;
				return 1;
			}
		}
	public:
		FastBPF() {
		}
		FastBPF(std::vector<int> list_of_Ns) {
			generate_transforms(list_of_Ns);
		}
		int get_no_cached_transforms() {
			return no_of_cached_transforms;
		}
		/* 
		note transforms are generated according to the N of the prototype filter
		i.e. a prototype with N = 2 gives a bandpass filter with N = 4 
		*/
		void generate_transforms(std::vector<int> list_of_Ns) {
			for (auto i = 0; i < list_of_Ns.size(); i++) {
				this->generate_transform(list_of_Ns[i]);
			}
		}
		/*
		convert an analog filter to a digital filter
		*/
		filtab<T> lp2bp(const filtab<T>& lpap, T wo, T Q) {
			filtab<T> bp;
			return this->lp2bp(lpap, bp, wo, Q);
		}
		filtab<T>& lp2bp(const filtab<T>& lpap, filtab<T>& bp, T wo, T Q) {
			using namespace std;
			int N = lpap.getN();
			bp.resize(2 * N);
			this->generate_transform(N);
			vector<vector<T>>& transform = cached_transforms[N];
			T k1 = 0, k2 = 0;
			for (auto i = 0; i < N + 1; i++) {
				k1 = pow(Q, N - i);
				for (auto j = 0; j < 2 * N + 1; j++) {
					k2 = pow(wo, i);
					bp.a[j] += k1 * lpap.a[i] * transform[i][j];
					bp.b[j] += k1 * lpap.b[i] * transform[i][j];
				}
			}
			lp2lp(bp, bp, wo);
			bp.normalize();
			return bp;
		}
	};
}