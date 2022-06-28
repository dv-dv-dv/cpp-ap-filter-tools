#include <iostream>
#include "filt.h"
int main()
{
    using namespace filt; using namespace std;
    filtab<double> buttap({ 0, 0, 1 }, { 1, 1.41421, 1 }); // create a second order butterworth LPAP

    // create a LPF with a cutoff of 1000 rad/s
    filtab<double> lpf = lp2lp(buttap, 1000.0);

    cout << "butterworth LPF coefficients with 1000 rad/s cutoff" << endl;
    lpf.print();

    // create a HPF with a cutoff of 1250 rad/s
    filtab<double> hpf(buttap.getN()); // create a filtab of the correct order
    lp2hp(buttap, hpf, 1250.0); // no memory is allocated/deallocated in this function 

    cout << "butterworth HPF coefficients with 1250 rad/s cutoff" << endl;
    hpf.print();

    // create a bpf with wl of 1000 and wu of 10000 rad/s
    filtab<double> bpf1 = lp2bp(buttap, 1000.0, 10000.0);
    cout << "coefficients of a butterworth BPF with lower cutoff of 1kHz and lower cutoff of 10kHz" << endl;
    bpf1.print();

    // create the same bpf but with the FastBPF class
    FastBPF<double> fastbpf({ 2 }); // instantiate FastBPF class and precompute transforms for N = 2
    filtab<double> bpf2(2 * buttap.getN()); // create a filtab of the correct order to store the bpf in
    bpf2 = fastbpf.lp2bp(buttap, bpf2, 1000.0, 10000.0); // no memory is allocated/deallocated in this function
    bpf2.print();

    // do bilinear transforms with frequency warping
    filtab<double> bpf_warped = lp2bp(buttap, bt_freq_warp(1000.0, 44100.0), bt_freq_warp(10000.0, 44100.0));
    filtab<double> lpf_warped = lp2lp(buttap, bt_freq_warp(1000.0, 44100.0));
    
    filtab<double> lpf_digital1 = bilin(lpf_warped, 44100.0);
    filtab<double> bpf_digital1 = bilin(bpf_warped, 44100.0);

    // do bilinear transform via FastBilin class 
    FastBilin<double> fb({ 2, 4 }); // instantiate the fast bilin class and precompute for N = 2, 4
    filtab<double> lpf_digital2 = fb.bilin(lpf_warped, 44100.0);
    filtab<double> bpf_digital2 = fb.bilin(bpf_warped, 44100.0);
    
    cout << "digital butterworth LPF coefficients with 1000kHz cutoff for fs = 44100Hz" << endl;
    lpf_digital1.print();
    lpf_digital2.print();

    cout << "coefficients of a digital butterworth BPF with lower cutoff of 1kHz and lower cutoff of 10kHz for fs = 44100Hz" << endl;
    bpf_digital1.print();
    bpf_digital2.print();
}
