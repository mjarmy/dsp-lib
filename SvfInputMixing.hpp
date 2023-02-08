#pragma once

#include <cmath>

//------------------------------------------------------------------------------
//
// SvfInputMixing is a State Variable Filter that provides a configurable mix of
// lowpass, bandpass and highpass filtering. It can also be used as a peaking
// filter or notch filter.
//
// SvfInputMixing is a C++ port of the algorithm described in the paper "Input
// mixing linear trapezoidal State Variable Filter (SVF) in state increment
// form", by Andy Simper of Cytomic.
//
// https://en.wikipedia.org/wiki/State_variable_filter
// https://cytomic.com/files/dsp/SvfInputMixing.pdf
//
//------------------------------------------------------------------------------

class SvfInputMixing {

  public:

    SvfInputMixing() {}

    ~SvfInputMixing() {}

    // init() changes the cutoff frequency, resonance, and/or sample rate.
    //
    // If you are changing the sample rate, you should call clear() before
    // calling init().
    //
    // cutoff: cutoff frequency -- In Hz. The range is [16.0, sampleRate/2.0].
    // res: resonance, aka Q -- The range is [0.0, 0.999].
    // sampleRate -- samples per second, e.g. 44100.
    //
    void init(double cutoff, double res, double sampleRate) {

        double k = 2 - 2 * res;
        double w = M_PI * cutoff / sampleRate;
        double s1 = std::sin(w);
        double s2 = std::sin(2 * w);
        double nrm = 1 / (2 + (2 - k) * s2);
        double s1n = 2 * s1 * s1 * nrm;
        double s2n = s2 * nrm;

        g0m0 = k * s1n + s2n;
        g0m1 = -s1n;
        g0m2 = -s2n;
        g1 = -s1n;
        g2 = -s2n;
        g3m0 = s1n;
        g3m1 = s2n;
        g3m2 = -s1n - k * s2n;
        g4 = s2n;
        g5 = -s1n - k * s2n;
    }

    // setMix() changes the output mixture of lowpass, bandpass and highpass
    // filtering.
    //
    // For many use cases, the sum of the three parameters should be 1.0, e.g.:
    //
    //      lowpass only: [1, 0, 0]
    //      bandpass only: [0, 1, 0]
    //      highpass only: [0, 0, 1]
    //      even mix of lowpass and highpass: [0.5, 0.0, 0.5]
    //      mostly lowpass plus a little bandpass: [0.9, 0.1, 0]
    //      etc, etc.
    //
    // You can also create a notch or peaking filter:
    //
    //      notch: [1, 0, 1]
    //      peaking: [-1, 0, 1]
    //
    void setMix(double lowMix_, double bandMix_, double highMix_) {
        lowMix = lowMix_;
        bandMix = bandMix_;
        highMix = highMix_;
    }

    // clear() resets the filter.
    void clear() {
        ic1eq = 0.0;
        ic2eq = 0.0;
    }

    // tick() filters the input based on the current settings.
    double tick(double input) {

        double vlow = lowMix * input;
        double vband = bandMix * input;
        double vhigh = highMix * input;

        double t1 =
            vlow * g0m0 + vband * g0m1 + vhigh * g0m2 + g1 * ic1eq + g2 * ic2eq;
        double t2 =
            vlow * g3m0 + vband * g3m1 + vhigh * g3m2 + g4 * ic1eq + g5 * ic2eq;
        double vc2 = t2 + ic2eq;
        ic1eq = ic1eq + 2.0 * t1;
        ic2eq = ic2eq + 2.0 * t2;

        return vhigh + vc2;
    }

  private:

    double g0m0 = 0.0;
    double g0m1 = 0.0;
    double g0m2 = 0.0;
    double g1 = 0.0;
    double g2 = 0.0;
    double g3m0 = 0.0;
    double g3m1 = 0.0;
    double g3m2 = 0.0;
    double g4 = 0.0;
    double g5 = 0.0;

    double lowMix = 0.0;
    double bandMix = 0.0;
    double highMix = 0.0;

    double ic1eq = 0.0;
    double ic2eq = 0.0;
};
