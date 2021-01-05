#!/bin/sh
echo "Running Octave nfsft tests..."
"/usr/bin/octave-cli" --eval "try; addpath('/LOCAL/tovo/nfft-3.5.2/matlab/tests','/LOCAL/tovo/nfft-3.5.2/matlab/nfsft'); nfsftUnitTestsRunAndExit; catch; disp('Error running nfsftUnitTestsRunAndExit'); end; exit(1);"
