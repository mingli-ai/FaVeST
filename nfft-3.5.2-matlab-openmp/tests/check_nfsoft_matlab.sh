#!/bin/sh
echo "Running MATLAB nfsoft tests..."
"/matlab" -wait -nodesktop -nosplash -r "try; diary('check_nfsoft_matlab.output'); addpath('/LOCAL/tovo/nfft-3.5.2/matlab/tests','/LOCAL/tovo/nfft-3.5.2/matlab/nfsft','/LOCAL/tovo/nfft-3.5.2/matlab/nfsoft'); perform_exhaustive_tests_flag=0; nfsoftUnitTestsRunAndExit; catch; disp('Error running nfsoftUnitTestsRunAndExit'); end; exit(1);"
