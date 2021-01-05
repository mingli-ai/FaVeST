#!/bin/sh
echo "Running Octave nfsoft tests..."
"/usr/bin/octave-cli" --eval "try; addpath('/LOCAL/tovo/nfft-3.5.2/matlab/tests','/LOCAL/tovo/nfft-3.5.2/matlab/nfsft','/LOCAL/tovo/nfft-3.5.2/matlab/nfsoft'); perform_exhaustive_tests_flag=0; nfsoftUnitTestsRunAndExit; catch; disp('Error running nfsoftUnitTestsRunAndExit'); end; exit(1);"
