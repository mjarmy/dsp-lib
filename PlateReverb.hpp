/*
MIT License

Copyright (c) 2023 Mike Jarmy

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#pragma once

#include <array>
#include <cmath>
#include <memory>

//-----------------------------------------------------------------------------
// PlateReverb is an algorithmic reverb that simulates a plate reverb.
//
// PlateReverb is an implementation of the classic plate reverb algorithm
// described by Jon Dattorro.
//
// Dattorro, J. 1997. "Effect Design Part 1: Reverberators and Other Filters."
// Journal of the Audio Engineering Society, Vol. 45, No. 9
//
// https://en.wikipedia.org/wiki/Reverb_effect#Plate_reverb
// https://ccrma.stanford.edu/~dattorro/EffectDesignPart1.pdf
//
// Parameters:
//
//    mix:        Dry/wet mix.
//    predelay:   Delay before reverb.
//    lowpass:    Apply a lowpass filter before reverb.
//    decay:      How quickly the reverb decays.
//    size:       The size of our imaginary plate.
//    damping:    How quickly high frequencies decay during reverb.
//
//-----------------------------------------------------------------------------

class PlateReverb {

  public:

    static constexpr double kMaxPredelay = 0.1; // seconds
    static constexpr double kMaxSize = 2.0;

    PlateReverb() {}
    ~PlateReverb() {}

    void setSampleRate(double sampleRate_) {
        sampleRate = sampleRate_;

        // Ratio of our sample rate to the sample rate that is used in
        // Dattorro's paper.
        double r = sampleRate / 29761.0;

        // Predelay
        predelayLine.reset(new DelayLine(std::ceil(sampleRate * kMaxPredelay)));

        // Lowpass filters
        lowpass.setSampleRate(sampleRate);
        leftTank.damping.setSampleRate(sampleRate);
        rightTank.damping.setSampleRate(sampleRate);

        // Diffusers
        diffusers[0].reset(new DelayAllpass(std::ceil(142 * r), 0.75));
        diffusers[1].reset(new DelayAllpass(std::ceil(107 * r), 0.75));
        diffusers[2].reset(new DelayAllpass(std::ceil(379 * r), 0.625));
        diffusers[3].reset(new DelayAllpass(std::ceil(277 * r), 0.625));

        // Tanks
        double maxModDepth = 8.0 * kMaxSize * r;
        leftTank.resetDelayLines(
            std::ceil(kMaxSize * 672 * r), -0.7, // apf1
            maxModDepth,
            std::ceil(kMaxSize * 4453 * r),      // del1
            std::ceil(kMaxSize * 1800 * r), 0.5, // apf2
            std::ceil(kMaxSize * 3720 * r)       // del2
        );
        rightTank.resetDelayLines(
            std::ceil(kMaxSize * 908 * r), -0.7, // apf1
            maxModDepth,
            std::ceil(kMaxSize * 4217 * r),      // del1
            std::ceil(kMaxSize * 2656 * r), 0.5, // apf2
            std::ceil(kMaxSize * 3163 * r)       // del2
        );

        leftTank.lfo.setSampleRate(sampleRate);
        rightTank.lfo.setSampleRate(sampleRate);
        leftTank.lfo.setFrequency(1.0);
        rightTank.lfo.setFrequency(0.95);

        // Tap points
        baseLeftTaps = {
            266 * r,  // rightTank.del1
            2974 * r, // rightTank.del1
            1913 * r, // rightTank.apf2
            1996 * r, // rightTank.del2
            1990 * r, // leftTank.del1
            187 * r,  // leftTank.apf2
            1066 * r, // leftTank.del2
        };
        baseRightTaps = {
            353 * r,  // leftTank.del1
            3627 * r, // leftTank.del1
            1228 * r, // leftTank.apf2
            2673 * r, // leftTank.del2
            2111 * r, // rightTank.del1
            335 * r,  // rightTank.apf2
            121 * r,  // rightTank.del2
        };
    }

    void setMix(double m /* [0, 1] */) { mix = clamp(m, 0.0, 1.0); }

    void setPredelay(double pd /* in seconds, [0, 0.1] */) {
        predelay = clamp(pd, 0.0, kMaxPredelay) * sampleRate;
    }

    void setLowpass(double cutoff /* Hz */) {
        cutoff = clamp(cutoff, 16.0, 20000.0);
        lowpass.setCutoff(cutoff);
    }

    void setDecay(double dr /* [0, 1) */) {
        decayRate = clamp(dr, 0.0, 0.9999999);
        leftTank.setDecay(decayRate);
        rightTank.setDecay(decayRate);
    }

    // The size parameter scales the delay time for all of the delay lines and
    // APFs in each tank, and for all of the tap points.
    //
    // Note that there is no size parameter in Dattorro's paper; it is an
    // extension to the original algorithm.
    void setSize(double sz /* [0, 2] */) {

        double sizeRatio = clamp(sz, 0.0, kMaxSize) / kMaxSize;

        // Scale the tank delays and APFs in each tank
        leftTank.setSizeRatio(sizeRatio);
        rightTank.setSizeRatio(sizeRatio);

        // Scale the taps
        for (int i = 0; i < kNumTaps; i++) {
            leftTaps[i] = baseLeftTaps[i] * sizeRatio;
            rightTaps[i] = baseRightTaps[i] * sizeRatio;
        }
    }

    void setDamping(double cutoff /* Hz */) {
        cutoff = clamp(cutoff, 16.0, 20000.0);

        leftTank.damping.setCutoff(cutoff);
        rightTank.damping.setCutoff(cutoff);
    }

    void process(
        double dryLeft, double dryRight, double* leftOut, double* rightOut) {

        double sum = dryLeft + dryRight;

        // predelay
        sum = predelayLine->tapAndPush(predelay, sum);

        // input lowpass
        sum = lowpass.process(sum);

        // diffusers
        sum = diffusers[0]->process(sum, diffusers[0]->getSize());
        sum = diffusers[1]->process(sum, diffusers[1]->getSize());
        sum = diffusers[2]->process(sum, diffusers[2]->getSize());
        sum = diffusers[3]->process(sum, diffusers[3]->getSize());

        // tanks
        double leftIn = sum + rightTank.out * decayRate;
        double rightIn = sum + leftTank.out * decayRate;
        leftTank.process(leftIn);
        rightTank.process(rightIn);

        // tap for output
        double wetLeft = rightTank.del1->tap(leftTaps[0])   //  266
                         + rightTank.del1->tap(leftTaps[1]) // 2974
                         - rightTank.apf2->tap(leftTaps[2]) // 1913
                         + rightTank.del2->tap(leftTaps[3]) // 1996
                         - leftTank.del1->tap(leftTaps[4])  // 1990
                         - leftTank.apf2->tap(leftTaps[5])  //  187
                         - leftTank.del2->tap(leftTaps[6]); // 1066

        double wetRight = leftTank.del1->tap(rightTaps[0])     //  353
                          + leftTank.del1->tap(rightTaps[1])   // 3627
                          - leftTank.apf2->tap(rightTaps[2])   // 1228
                          + leftTank.del2->tap(rightTaps[3])   // 2673
                          - rightTank.del1->tap(rightTaps[4])  // 2111
                          - rightTank.apf2->tap(rightTaps[5])  //  335
                          - rightTank.del2->tap(rightTaps[6]); //  121

        // mix
        *leftOut = dryLeft * (1 - mix) + wetLeft * mix;
        *rightOut = dryRight * (1 - mix) + wetRight * mix;
    }

  private:

    //--------------------------------------------------------------
    // OnePoleFilter is a one-pole low pass filter.
    //--------------------------------------------------------------

    class OnePoleFilter {

      public:

        OnePoleFilter() {}
        ~OnePoleFilter() {}

        void setSampleRate(double sampleRate_) {
            invSampleRate = 1 / sampleRate_;
            recalc();
        }

        void setCutoff(double cutoff_ /* Hz */) {
            cutoff = cutoff_;
            recalc();
        }

        double process(double x) {
            z = x * a + z * b;
            return z;
        }

      private:

        double invSampleRate = 1;
        double cutoff = 0;

        double a = 0;
        double b = 0;
        double z = 0;

        void recalc() {
            b = std::exp(-2 * M_PI * cutoff * invSampleRate);
            a = 1 - b;
        }
    };

    //--------------------------------------------------------------
    // DelayLine
    //--------------------------------------------------------------

    class DelayLine {

      public:

        DelayLine(uint64_t size_) : size(size_) {

            // For speed, create a bigger buffer than we really need.
            uint64_t bufferSize = ceilPowerOfTwo(size);
            buffer.reset(new double[bufferSize]);
            std::memset(&buffer[0], 0, bufferSize * sizeof(double));

            mask = bufferSize - 1;

            writeIdx = 0;
        }

        ~DelayLine() {}

        inline void push(double val) {
            buffer[writeIdx++] = val;
            writeIdx &= mask;
        }

        inline double tap(double delay /* samples */) {
            // We always want to be able to properly handle any delay value that
            // gets passed in here, without going past the original size.
            assert(delay <= size);

            int64_t d = delay;
            double frac = 1 - (delay - d);

            int64_t readIdx = (writeIdx - 1) - d;
            double a = buffer[(readIdx - 1) & mask];
            double b = buffer[readIdx & mask];

            return a + (b - a) * frac;
        }

        // This does "read-before-write".
        inline double tapAndPush(double delay, double val) {
            double out = tap(delay);
            push(val);
            return out;
        }

        inline uint64_t getSize() { return size; }

      private:

        const uint64_t size;

        std::unique_ptr<double[]> buffer;
        uint64_t mask;

        uint64_t writeIdx;

        static uint64_t ceilPowerOfTwo(uint64_t n) {
            return (uint64_t)std::pow(2, std::ceil(std::log(n) / std::log(2)));
        }
    };

    //------------------------------------------
    // DelayAllpass
    //------------------------------------------

    class DelayAllpass {

      public:

        DelayAllpass(uint64_t size_, double gain_)
            : delayLine(size_), gain(gain_) {}

        ~DelayAllpass() {}

        inline double process(double x, double delay) {
            double wd = delayLine.tap(delay);
            double w = x + gain * wd;
            double y = -gain * w + wd;
            delayLine.push(w);
            return y;
        }

        inline void setGain(double gain_) { gain = gain_; }

        inline double tap(double delay) { return delayLine.tap(delay); }

        inline uint64_t getSize() { return delayLine.getSize(); }

      private:

        DelayLine delayLine;
        double gain;
    };

    //--------------------------------------------------------------
    // Lfo
    //--------------------------------------------------------------

    class Lfo {

      public:

        Lfo() {}
        ~Lfo() {}

        void setSampleRate(double sampleRate_) {
            invSampleRate = 1 / sampleRate_;
            recalc();
        }

        void setFrequency(double freq_) {
            freq = freq_;
            recalc();
        }

        inline double process() {
            double out = -fastSin(phase);

            phase += phaseInc;
            if (phase > M_PI) {
                phase = -M_PI;
            }

            return out;
        }

      private:

        double invSampleRate = 1;
        double freq = 0;

        double phaseInc = 0;
        double phase = -M_PI;

        void recalc() {
            phaseInc = freq * invSampleRate;
            phaseInc *= 2 * M_PI;
        }

        // Parabolic sine approximation.
        // Range is [-pi, pi].
        //
        // https://web.archive.org/web/20100613230051/http://www.devmaster.net/forums/showthread.php?t=5784
        inline double fastSin(double x) {
            static constexpr double B = 4 / M_PI;
            static constexpr double C = -4 / (M_PI * M_PI);
            static constexpr double P = 0.225;

            double y = B * x + C * x * std::abs(x);

            // Extra precision.
            y = P * (y * std::abs(y) - y) + y;

            return y;
        }
    };

    //------------------------------------------
    // Tank
    //------------------------------------------

    class Tank {

      public:

        Tank() {}
        ~Tank() {}

        void resetDelayLines(
            uint64_t apf1Size_, double apf1Gain_, // First APF
            double maxModDepth_,
            uint64_t delay1Size_,                 // First delay
            uint64_t apf2Size_, double apf2Gain_, // Second APF
            uint64_t delay2Size_                  // Second delay
        ) {
            apf1Size = apf1Size_;
            maxModDepth = maxModDepth_;
            double maxApf1Size = apf1Size + maxModDepth + 1;
            apf1.reset(new DelayAllpass(maxApf1Size, apf1Gain_));

            del1.reset(new DelayLine(delay1Size_));
            apf2.reset(new DelayAllpass(apf2Size_, apf2Gain_));
            del2.reset(new DelayLine(delay2Size_));

            // We've changed the various delay line sizes and associated values,
            // so update the sizeRatio values too.
            recalcSizeRatio();
        }

        void setDecay(double decayRate_) {
            decayRate = decayRate_;
            apf2->setGain(clamp(decayRate + 0.15, 0.25, 0.5));
        }

        void setSizeRatio(double sizeRatio_) {
            sizeRatio = sizeRatio_;
            recalcSizeRatio();
        }

        void process(double val) {

            // APF1: "Controls density of tail."
            val = apf1->process(val, apf1Delay + lfo.process() * modDepth);
            val = del1->tapAndPush(del1Delay, val);

            val = damping.process(val);
            val *= decayRate;

            // APF2: "Decorrelates tank signals."
            val = apf2->process(val, apf2Delay);
            val = del2->tapAndPush(del2Delay, val);

            out = val;
        }

        double out = 0.0;

        std::unique_ptr<DelayAllpass> apf1 = nullptr;
        std::unique_ptr<DelayAllpass> apf2 = nullptr;
        std::unique_ptr<DelayLine> del1 = nullptr;
        std::unique_ptr<DelayLine> del2 = nullptr;
        OnePoleFilter damping;
        Lfo lfo;

      private:

        uint64_t apf1Size;
        double maxModDepth = 0;
        double modDepth = 0;

        double apf1Delay = 0;
        double apf2Delay = 0;
        double del1Delay = 0;
        double del2Delay = 0;

        double decayRate = 0;
        double sizeRatio = 0;

        void recalcSizeRatio() {

            apf1Delay = apf1Size * sizeRatio;
            modDepth = maxModDepth * sizeRatio;

            apf2Delay = apf2->getSize() * sizeRatio;
            del1Delay = del1->getSize() * sizeRatio;
            del2Delay = del2->getSize() * sizeRatio;
        }
    };

    //--------------------------------------------------------------
    //--------------------------------------------------------------
    //--------------------------------------------------------------

    double sampleRate = 1.0;

    double mix = 0.0;
    double predelay = 0.0;
    double decayRate = 0.0;

    std::unique_ptr<DelayLine> predelayLine = nullptr;
    OnePoleFilter lowpass;
    std::array<std::unique_ptr<DelayAllpass>, 4> diffusers;

    Tank leftTank;
    Tank rightTank;

    static const int kNumTaps = 7;
    std::array<double, kNumTaps> baseLeftTaps = {};
    std::array<double, kNumTaps> baseRightTaps = {};
    std::array<double, kNumTaps> leftTaps = {};
    std::array<double, kNumTaps> rightTaps = {};

    static inline double clamp(double val, double low, double high) {
        return std::min(std::max(val, low), high);
    }
};
