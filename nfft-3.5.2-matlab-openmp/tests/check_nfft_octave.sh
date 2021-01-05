#!/bin/sh
echo "Running Octave nfft tests..."
"/usr/bin/octave-cli" --eval "try; addpath('/LOCAL/tovo/nfft-3.5.2/matlab/tests','/LOCAL/tovo/nfft-3.5.2/matlab/nfft'); nfftUnitTestsRunAndExit; catch; disp('Error running nfftUnitTestsRunAndExit'); end; exit(1);"
