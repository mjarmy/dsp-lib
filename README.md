# dsp-lib

`dsp-lib` is a collection of utilities related to audio DSP

### PlateReverb

PlateReverb is an algorithmic reverb that simulates 
a [plate reverb](https://en.wikipedia.org/wiki/Reverb_effect#Plate_reverb).

PlateReverb is an implementation of the classic plate reverb algorithm
described by Jon Dattorro.

Dattorro, J. 1997. ["Effect Design Part 1: Reverberators and Other Filters."](https://ccrma.stanford.edu/~dattorro/EffectDesignPart1.pdf)
Journal of the Audio Engineering Society, Vol. 45, No. 9

### SvfInputMixing

`SvfInputMixing` is a [State Variable Filter](https://en.wikipedia.org/wiki/State_variable_filter) 
that provides a configurable mix of lowpass, bandpass and highpass filtering. It
can also be used as a peaking filter or notch filter.

`SvfInputMixing` is a C++ port of the algorithm described in the [paper](https://cytomic.com/files/dsp/SvfInputMixing.pdf) 
"Input mixing linear trapezoidal State Variable Filter (SVF) in state increment
form", by Andy Simper of Cytomic.
