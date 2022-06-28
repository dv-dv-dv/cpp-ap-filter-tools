#include <iostream>
#include <vector>
#include <unordered_map>
#include <ctgmath>
namespace filt {
	template <typename T> 
	T choose(T n, T k) {
		T a1 = 1, a2 = 1, a3 = 1;
		for (auto i = 2; i <= n; i++) {
			a1 *= i;
			if (i <= k) {
				a2 *= i;
			}
			if (i <= n - k) {
				a3 *= i;
			}
		}
		return a1 / (a2 * a3);
	}
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
	transforms an analog lpf with a cutoff of 1 rad/s to a bpf with a lower cutoff of wl and an upper cutoff of wu
	*/
	template <typename T>
	filtab<T> lp2bp(const filtab<T>& lpap, T wl, T wu) {
		filtab<T> bp;
		return lp2bp(lpap, bp, wl, wu);
	}
	template <typename T>
	filtab<T>& lp2bp(const filtab<T>& lpap, filtab<T>& bp, T wl, T wu) {
		T wo = sqrt(wu * wl);
		T Q = wo / (wu - wl);
		int N = lpap.getN();
		T k = 0, k1, k2;
		bp.resize(2 * N);
		for (auto i = 0; i < N + 1; i++) {
			for (auto j = 0; j < i + 1; j++) {
				if ((i + j) % 2 == 0) {
					k = choose<T>(N - j, (i - j) / 2);
					k1 = pow(Q, -j);
					k2 = pow(wo, i);
					bp.a[i] += lpap.a[j] * k * k1 * k2;
					bp.b[i] += lpap.b[j] * k * k1 * k2;
					if (i != 2 * N - i){
						k2 = pow(wo, (2 * N - i));
						bp.a[2 * N - i] += lpap.a[j] * k * k1 * k2;
						bp.b[2 * N - i] += lpap.b[j] * k * k1 * k2;
					}
				}
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
		using namespace std;
		int N = analog.getN(), dim = N + 1, r1 = 0, r2 = 0, sign = 1;
		//test.assign((N + 1) * 2, 1);
		vector<T> test(2 * (N + 1), 1);
		digital.resize(N);
		T k = 2 * fs, k2 = 0;
		for (auto i = 0; i < dim; i++) {
			k2 = pow(k, -i);
			digital.a[0] += k2 * analog.a[i];
			digital.b[0] += k2 * analog.b[i];
		}
		for (auto i = 1; i < dim; i++) {
			r1 = i % 2;
			r2 = (i + 1) % 2;
			sign *= -1;
			test[r1 * dim] = choose(N, i) * sign;
			digital.a[i] += analog.a[0] * test[r1 * dim];
			digital.b[i] += analog.b[0] * test[r1 * dim];
			for (auto j = 1; j < dim; j++) {
				test[r1 * dim + j] = test[r1 * dim + j - 1] + test[r2 * dim + j - 1] + test[r2 * dim + j];
				k2 = pow(k, -j);
				digital.a[i] += k2 * analog.a[j] * test[r1 * dim + j];
				digital.b[i] += k2 * analog.b[j] * test[r1 * dim + j];
			}
		}
		digital.normalize();
		return digital;
	}
	/*
	bilinear transform specialization for filter order of 1
	*/
	template <typename T>
	filtab<T> bilin1(const filtab<T>& analog1, T fs) {
		filtab<T> digital1;
		return bilin1(analog1, digital1, fs);
	}
	template <typename T>
	filtab<T>& bilin1(const filtab<T>& analog1, filtab<T>& digital1, T fs) {
		digital1.resize(1);
		T k = 2 * fs;
		digital1.a[0] = analog1.a[0] * k + analog1.a[1];
		digital1.a[1] = -analog1.a[0] * k + analog1.a[1];

		digital1.b[0] = analog1.b[0] * k + analog1.b[1];
		digital1.b[1] = -analog1.b[0] * k + analog1.b[1];
		digital1.normalize();
		return digital1;
	}
	/*
	bilinear transform specialization for filter order of 2
	*/
	template <typename T>
	filtab<T> bilin2(const filtab<T>& analog2, T fs) {
		filtab<T> digital2;
		return bilin2(analog2, digital2, fs);
	}
	template <typename T>
	filtab<T>& bilin2(const filtab<T>& analog2, filtab<T>& digital2, T fs) {
		digital2.resize(2);
		T k = 2 * fs; T k2 = k * k;
		digital2.a[0] = analog2.a[0] * k2 + analog2.a[1] * k + analog2.a[2];
		digital2.a[1] = -2 * analog2.a[0] * k2 + 2 * analog2.a[2];
		digital2.a[2] = analog2.a[0] * k2 - analog2.a[1] * k + analog2.a[2];

		digital2.b[0] = analog2.b[0] * k2 + analog2.b[1] * k + analog2.b[2];
		digital2.b[1] = -2 * analog2.b[0] * k2 + 2 * analog2.b[2];
		digital2.b[2] = analog2.b[0] * k2 - analog2.b[1] * k + analog2.b[2];

		digital2.normalize();
		return digital2;
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
}