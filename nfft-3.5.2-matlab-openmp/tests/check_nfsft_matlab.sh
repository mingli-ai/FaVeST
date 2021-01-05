#!/bin/sh
echo "Running MATLAB nfsft tests..."
"/matlab" -wait -nodesktop -nosplash -r "try; diary('check_nfsft_matlab.output'); addpath('/LOCAL/tovo/nfft-3.5.2/matlab/tests','/LOCAL/tovo/nfft-3.5.2/matlab/nfsft'); nfsftUnitTestsRunAndExit; catch; disp('Error running nfsftUnitTestsRunAndExit'); end; exit(1);"
