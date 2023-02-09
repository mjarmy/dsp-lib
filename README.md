# dsp-lib

`dsp-lib` is a collection of utilities related to audio DSP

### SvfInputMixing

`SvfInputMixing` is a [State Variable Filter](https://en.wikipedia.org/wiki/State_variable_filter) 
that provides a configurable mix of lowpass, bandpass and highpass filtering. It
can also be used as a peaking filter or notch filter.

`SvfInputMixing` is a C++ port of the algorithm described in the [paper](https://cytomic.com/files/dsp/SvfInputMixing.pdf) 
"Input mixing linear trapezoidal State Variable Filter (SVF) in state increment
form", by Andy Simper of Cytomic.
