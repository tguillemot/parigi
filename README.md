# About
  - Author    :: Thierry Guillemot <thierry.guillemot.work@gmail.com>,
                 Julie Delon <julie.delon@parisdescartes.fr> and Agnès
                 Desolneux <agnes.desolneux@cmla.ens-cachan.fr>
  - Copyright :: (C) 2015 IPOL Image Processing On Line http://www.ipol.im/
  - Licence   :: GPL V3+

  Copying and distribution of this file, with or without modification,
  are permitted in any medium without royalty provided the copyright
  notice and this notice are preserved.  This file is offered as-is,
  without any warranty.
# Overview
  This source code provides an implementation of the "Parigi": a
  patch-based approach to remove impulse-gaussian noise image, as
  described in IPOL.

  The 'bin/parigi' program reads an PNG image, add impulse-gaussian
  noise and then apply the parigi denoising algorithm. The algorithm
  estimates the probability of impulse noise and extract similar
  patches. From them, it estimates for each pixels an histogram of all
  possible values. The denoised image is defined from the best mode of
  each histograms.

  Only 8bit GRAY PNG images are handled.
# Requirement
  The code is written in UTF8 C++, and should compile on any system with
  an UTF8 C++ compiler.

  The libpng and FFTW header and libraries are required on the system
  for compilation and execution. On Linux, just use your package
  manager to install it:
```
sudo apt-get install libpng
sudo apt-get install libfftw3-dev
```


  For more information, see http://www.libpng.org/pub/png/libpng.html
  and http://www.fftw.org/.

# Compilation
  Simply use the provided makefile, with the command 'make'.  The
  makefile will produce a program called : 'bin/parigi'.

  It is possible to compile the program using OpenMP with the command
  'make OMP=1'.

  The 'parigi' program is used to corrupted an original image with an
  impulse-gaussian noise and the apply the parigi restoration.

# Usage
## Parigi
### Description
  The 'parigi' program is used to corrupted an original image with an
  impulse-gaussian noise and the apply the parigi restoration.
  It takes 8 parameters with other optional:
```
bin/parigi input.png noisy.png denoised.png p_gaussian_noise p_impulse_noise two_step mixed_image verbose
```
  - input.png  :: input image.
  - noisy.png  :: resulting noisy image corrupted with an impulse-gaussian noise.
  - denoised.png :: resulting parigi denoised image.
  - p_gaussian_noise :: noise gaussian standard deviation [0. 255.0].
  - p_impulse_noise :: probability of the impulse noise [0.0, 1.0].
  - two_step :: apply parigi algorithm two times (allows to obtain best visual results) [0 or 1].
  - mixed_image :: apply a criterium to mix the noisy and the denoised images (allow to obtain best visual results) [0 or 1].
  - verbose :: activate the verbose mode

### Usage examples
  If you want to compute the parigi restoration of the file
  'input/barbara.png' corrupted with a gaussian noise of standard
  deviation of 10 and an impulse noise with a probability of 40% by
  applying two steps with the mixed criterium without the verbose
  mode and saving the noisy image into 'noisy.png' and the denoised
  image into 'denoised.png', you can use the command:
```
bin/parigi input/barbara.png noisy.png denoised.png 10 0.4 1 1 0
```

# Bugs reports
You can report any bug with the github interface:
https://github.com/tguillemot/parigi
