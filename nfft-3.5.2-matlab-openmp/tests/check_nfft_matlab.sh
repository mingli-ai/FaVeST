#!/bin/sh
echo "Running MATLAB nfft tests..."
"/matlab" -wait -nodesktop -nosplash -r "try; diary('check_nfft_matlab.output'); addpath('/LOCAL/tovo/nfft-3.5.2/matlab/tests','/LOCAL/tovo/nfft-3.5.2/matlab/nfft'); nfftUnitTestsRunAndExit; catch; disp('Error running nfftUnitTestsRunAndExit'); end; exit(1);"
