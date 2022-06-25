// filterdesign.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "filt.h"
int main()
{
    using namespace filt; using namespace std;
    filtab<double> buttap({ 0, 0, 1 }, { 1, 1.41421, 1 }); // create a second order butterworth LPAP
    filtab<double> lpf = lp2lp<double>(buttap, 1000);

    cout << "butterworth LPF coefficients with 1000 rad/s cutoff" << endl;
    lpf.print();

    filtab<double> hpf(buttap.getN());
    lp2hp<double>(buttap, hpf, 1250); // in this case lp2hp will not perform memory allocation, making this very fast

    cout << "butterworth HPF coefficients with 1250 rad/s cutoff" << endl;
    hpf.print();

    filtab<double> bpf1(2 * buttap.getN()), bpf2;
    FastBPF<double> fastbpf({ 1, 2, 3, 4, 5 }); // instantiate FastBPF class and precompute transforms for 1 <= N <= 5
    bpf1 = fastbpf.lp2bp(buttap, bpf1, 1000, 0.5); // because a transform for an AP of N = 2 was precalculated and bpf was instantiated with order 2N, this line did not allocate memory
    bpf2 = lp2bp(buttap, 1000.0, 0.5);

    cout << "butterworth BPF coefficients with 100 rad/s center and Q of 2" << endl;
    bpf1.print();
    bpf2.print();

    FastBilin<double> fb({ 1, 2, 3, 4, 5 }); // instantiate FastBilin class and precompute transforms for 1 <= N <= 5
    filtab<double> lpf_digital1 = fb.bilin(lp2lp<double>(buttap, bt_freq_warp<double>(1000, 44100)), 44100);
    filtab<double> lpf_digital2 = bilin<double>(lp2lp<double>(buttap, bt_freq_warp<double>(1000, 44100)), 44100);

    cout << "digital butterworth LPF coefficients with 1000 rad/s cutoff for fs = 44100" << endl;
    lpf_digital1.print();
    lpf_digital2.print();
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
