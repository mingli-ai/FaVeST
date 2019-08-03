% Download NFFT Package and Run Demo
clear;clc;
str = computer;
switch str
    case 'PCWIN'
        fprintf('Downloading NFFT package wtih multicore support OpenMP and 32bit ...\n');    
        url = 'https://www-user.tu-chemnitz.de/~potts/nfft/windows/nfft-3.5.0-mexw32-openmp.zip';
        untar(url);
        fprintf('Succefully installed NFFT package!\n'); 
        fprintf('==========================================================================\n'); 
        fprintf('****** Demo for Reconstruction by FaVeST ****** \n');
        run Demo.m
    case 'PCWIN64'
         fprintf('Downloading NFFT package wtih multicore support OpenMP and 64bit ...\n');     
         url = 'https://www-user.tu-chemnitz.de/~potts/nfft/windows/nfft-3.5.0-mexw64-openmp.zip';
         unzip(url);
         fprintf('Succefully installed NFFT package!\n'); 
         fprintf('==========================================================================\n'); 
         fprintf('****** Demo for Reconstruction by FaVeST ****** \n');
         run Demo.m
    case 'GLNX86'
        warning('Unmatched Operating System, please download the correct version of NFFT package from the the website:\n');
        fprintf('https://www-user.tu-chemnitz.de/~potts/nfft/download.php');
    case 'GLNXA64'
         fprintf('Downloading NFFT package wtih multicore support OpenMP and 64bit ...\n'); 
         url = 'https://www-user.tu-chemnitz.de/~potts/nfft/linux/nfft-3.5.0-mexa64-octave-5.1-openmp.tar.gz';
         unzip(url);
         fprintf('Succefully installed NFFT package!\n'); 
         fprintf('==========================================================================\n'); 
         fprintf('****** Demo for Reconstruction by FaVeST ****** \n');
         run Demo.m
    case 'MACI'
         warning('Unmatched Operating System, please download the correct version of NFFT package from the the website:\n');
         fprintf('https://www-user.tu-chemnitz.de/~potts/nfft/download.php');
    case 'MACI64'
         fprintf('Downloading NFFT package wtih multicore support OpenMP and 64bit ...\n'); 
         url = 'https://www-user.tu-chemnitz.de/~potts/nfft/macos/nfft-3.5.0-mexmaci64-openmp.zip';
         unzip(url);
         fprintf('Succefully installed NFFT package!\n'); 
         fprintf('==========================================================================\n'); 
         fprintf('****** Demo for Reconstruction by FaVeST ****** \n');
         run Demo.m
    otherwise
        warning('Unmatched Operating System, please download the right version of NFFT package in the the website:\n');
        fprintf('https://www-user.tu-chemnitz.de/~potts/nfft/download.php');
end
