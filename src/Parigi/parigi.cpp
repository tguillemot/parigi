/*
 * Copyright (c) 2015, Thierry Guillemot
 * <thierry.guillemot.work@gmail.com> All rights reserved.
 *
 * This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

#include "Parigi/parigi.h"

#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <fftw3.h>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

int LaunchParigi(
    const char  *input_filename,
    const char  *noisy_filename,
    const char  *denoised_filename,
    const float &p_gaussian_noise,
    const float &p_impulse_noise,
    const bool  &two_iteration,
    const bool  &mixed_image,
    const bool  &verbose
    )
{
  //! Parameters
  const size_t maxIter                = two_iteration ? 2 : 1;
  const size_t half_size_patch        = 3;
  const size_t half_size_neighborhood = 5;


  size_t nbThreads=1;
#ifdef _OPENMP
  nbThreads = omp_get_max_threads();
#endif

  //! Load images
  vector<float> im_input;
  ImageSize imSize;
  if (loadImage(input_filename, im_input, imSize, verbose) != EXIT_SUCCESS)
  {
    cerr << "Failed reading image " << input_filename << endl;
    return EXIT_FAILURE;
  }

  //! Compute the noisy image by adding an impulse noise
  vector<float> im_noisy;
  addGaussianNoise(im_input,  im_noisy, p_gaussian_noise,  verbose);
  addImpulseNoise(im_noisy, im_noisy, p_impulse_noise, verbose);
  if (saveImage(noisy_filename, im_noisy, imSize, 0.f, 255.f) != EXIT_SUCCESS)
  {
    cerr << "Failed writing image " << noisy_filename << endl;
    return EXIT_FAILURE;
  }

  //! Apply Nstep
  vector<float> u_i = im_noisy;
  vector<float> roadmask;
  float p_estimated;
  for(size_t iter=0; iter<maxIter; ++iter)
  {
    cout << "Applying iteration " << iter+1 << "..." << endl;

    //! Compute the noise probability
    estimp_road(u_i, imSize, &roadmask, &p_estimated);

    //! Estimate the number of patches necessary at each pixel for
    //! impulse denoising
    const int number_patch = estim_npatchs(p_estimated);

    //! Compute the beta coefficient beta in a lookup table
    const vector<float> coeff = ComputeIncompleteBetaFunction((1-p_estimated)*(1-p_estimated), half_size_patch);

    //! Extraction of the sub images to be used by OpenMP
    vector<vector<float> > u_i_sub;
    ImageSize              imSize_sub;
    SubdiveImage(u_i, imSize, 2*nbThreads, 2*half_size_patch + half_size_neighborhood, &u_i_sub, &imSize_sub);

    //!Computation of the denoised image
    vector<vector<float> > hat_mode_i_sub(u_i_sub.size());
    vector<vector<float> > hat_sigma_i_sub(u_i_sub.size());
    #pragma omp parallel for schedule(dynamic, 2) shared(u_i_sub, imSize_sub)
    for(size_t k=0; k<u_i_sub.size(); ++k)
    {
      RestoreImpulseNoise(u_i_sub[k], imSize_sub,
                            half_size_patch, half_size_neighborhood,
                            p_estimated, number_patch, coeff,
                            &hat_mode_i_sub[k], &hat_sigma_i_sub[k]);
    }

    //! Recompose the final image from the sub images after
    vector<float> sigma_i = RecomposeImageWithBoundaries(hat_sigma_i_sub, imSize, 2*nbThreads, half_size_patch);
    vector<float> mode_i = RecomposeImageWithBoundaries(hat_mode_i_sub, imSize, 2*nbThreads, half_size_patch);

    if(mixed_image) ComputeMixedCriterium(u_i, mode_i, sigma_i, &u_i);
    else u_i = mode_i;
  }

  if (saveImage(denoised_filename, u_i, imSize, 0.f, 255.f) != EXIT_SUCCESS)
  {
    cerr << "Failed writing image " << denoised_filename << endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

//! Private function to compare the pair
template<typename A, typename B>
bool pairCompare(pair<A, B>& firstElem, const std::pair<A, B>& secondElem)
{
  return firstElem.first < secondElem.first;
}

int RestoreImpulseNoise(
    const vector<float>  &im_noisy_bound,
    const ImageSize      &imSize_bound,
    const float          &half_size_patch,
    const float          &half_size_neighborhood,
    const float          &p_estimated,
    const int            &number_patch,
    const vector<float>  &coeff,
    vector<float>        *modes,
    vector<float>        *sigma
    )
{
  //! Constant used by the function
  const int length_patch = 2*half_size_patch + 1;
  const int length_neigh = 2*half_size_neighborhood + 1;
  const int dim_patch    = length_patch * length_patch;
  const int dim_neigh    = length_neigh * length_neigh;
  const int offset       = half_size_neighborhood + half_size_patch;
  //! Image size less the patch length
  ImageSize imSize_neigh;
  imSize_neigh.height    = imSize_bound.height - 2*half_size_patch;
  imSize_neigh.width     = imSize_bound.width  - 2*half_size_patch;
  imSize_neigh.nChannels = imSize_bound.nChannels;
  imSize_neigh.wh        = imSize_neigh.height*imSize_neigh.width;
  imSize_neigh.whc       = imSize_neigh.wh*imSize_neigh.nChannels;
  //! Image size of the restore pixels
  ImageSize imSize;
  imSize.height          = imSize_bound.height - 2*offset;
  imSize.width           = imSize_bound.width  - 2*offset;
  imSize.nChannels       = imSize_bound.nChannels;
  imSize.wh              = imSize.height*imSize.width;
  imSize.whc             = imSize.wh*imSize.nChannels;

  //! Patches computation
  const vector<float> patches = ExtractCenteredPatches(im_noisy_bound, imSize_bound, half_size_patch);

   // Distance computation between patches
  vector<pair<float, const float*> > dist_patch = ComputeNearestPatches(coeff, patches, imSize_neigh, imSize, dim_patch, half_size_neighborhood, number_patch);

  // Histogram computation from these patches
  const size_t size_histogram = 1<<8;
  vector<float> histograms = ComputeHistograms(dist_patch, imSize, half_size_patch, size_histogram, number_patch, dim_neigh);

  //! Computation of the main mode of each histograms
  const float sigma_max = 25.f;
  ComputeModeFFTW(histograms, imSize.wh, size_histogram, p_estimated, sigma_max, modes, sigma);

  return EXIT_SUCCESS;
}

int ComputeMixedCriterium(
    const vector<float> &im_noisy,
    const vector<float> &modes,
    const vector<float> &sigma,
    vector<float>       *im_denoised)
{
  //! Restore only for noisy points
  im_denoised->resize(im_noisy.size());

  for(size_t i=0; i<im_noisy.size(); ++i)
  {
      if(abs(modes[i] - im_noisy[i])>sigma[i])
      {
        im_denoised->at(i) = modes[i];
      }
      else
      {
        im_denoised->at(i) = im_noisy[i];
      }
  }

  return EXIT_SUCCESS;
}

vector<pair<float, const float*> > ComputeNearestPatches(
    const vector<float> &coeff,
    const vector<float> &patches,
    const ImageSize     &imSize_neigh,
    const ImageSize     &imSize,
    const int           &dim_patch,
    const int           &half_size_neighborhood,
    const size_t        &number_patch
    )
{
  const int length_neigh = 2*half_size_neighborhood + 1;
  const int dim_neigh    = length_neigh * length_neigh;

  vector<pair<float, const float*> > nearest_patch(imSize.wh*dim_neigh);

  vector<float> dist(dim_patch);
  int c=0;
  for(int w=0; w<imSize.width; ++w)
  {
    for(int h=0; h<imSize.height; ++h)
    {
      const float* patch_ref = patches.data() + imSize_neigh.at(h+half_size_neighborhood, w+half_size_neighborhood, c)*dim_patch;
      for(int j=0; j<length_neigh; ++j)
      {
        for(int i=0; i<length_neigh; ++i)
        {
          const float* patch_comp = patches.data() + imSize_neigh.at(h+i, w+j, c)*dim_patch;

          for(int k=0; k<dim_patch; ++k)
          {
            dist[k] = (patch_ref[k]-patch_comp[k])*(patch_ref[k]-patch_comp[k]);
          }

          sort(dist.begin(),dist.end());

          float d = .0f;
          for(int k=0; k<dim_patch; ++k)
          {
            d += coeff[k]*dist[k];
          }
          int pos_dist = i + j*length_neigh + imSize.at(h, w, c)*dim_neigh;
          nearest_patch[pos_dist] = make_pair(d, patch_comp);
        }
      }

      partial_sort(nearest_patch.begin() + imSize.at(h,w,c)*dim_neigh,
                   nearest_patch.begin() + imSize.at(h,w,c)*dim_neigh + number_patch,
                   nearest_patch.begin() + (imSize.at(h,w,c)+1)*dim_neigh,
                   pairCompare<float, const float*>);
    }
  }

  return nearest_patch;
}

vector<float> ComputeHistograms(
    const vector<pair<float, const float*> > &dist_patch,
    const ImageSize                          &imSize,
    const int                                &half_length_patch,
    const size_t                             &size_histogram,
    const size_t                             &number_patch,
    const size_t                             &dim_neigh
    )
{
  const size_t length_patch = 2*half_length_patch+1;
  const float  minimum      = 0;
  const float  maximum      = 255;
  const float  alpha        = (size_histogram-1)/(maximum-minimum);

  vector<float> histograms(imSize.wh*size_histogram, .0f);
  int c=0;
  for(int w=0; w<imSize.width; ++w)
  {
    for(int h=0; h<imSize.height; ++h)
    {
      for(size_t k=0; k<number_patch; ++k)
      {
        const float* patch = dist_patch[k + imSize.at(h, w, c)*dim_neigh].second;
        for(int j=-half_length_patch; j<=half_length_patch; ++j)
        {
          for(int i=-half_length_patch; i<=half_length_patch; ++i)
          {
            if( (h+i<0) || (h+i>=imSize.height) || (w+j<0) || (w+j>=imSize.width) ) continue;
            int pos_hist = floor(alpha*(patch[i+half_length_patch + (j+half_length_patch)*length_patch] - minimum));
            ++histograms[pos_hist + imSize.at(h+i, w+j, c)*size_histogram];
          }
        }
      }
    }
  }

  return histograms;
}

int ComputeModeFFTW(
    const vector<float> &histograms,
    const size_t        &size,
    const size_t        &histogram_size,
    const float         &p,
    const int           &sigma_max,
    vector<float>       *modes,
    vector<float>       *scale
    )
{
  //! Output allocation
  vector<float> max_value(size);
  modes->resize(size);
  scale->resize(size);

  //! Initialisation of FFT
  fftw_plan     fft;
  vector<double>  data(histogram_size);
  fftw_complex *fft_data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * histogram_size);
#pragma omp critical
  fft = fftw_plan_dft_r2c_1d(histogram_size, data.data(), fft_data, FFTW_ESTIMATE);

  //! Compute Histogram FFT
  fftw_complex *fft_histograms = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * size*histogram_size);
  for(size_t k=0; k<size; ++k)
  {
    for(size_t l=0; l<histogram_size; ++l)
    {
      data[l] = histograms[l+k*histogram_size];
    }
    fftw_execute_dft_r2c(fft, data.data(), fft_data);
    for(size_t l=0; l<histogram_size; ++l)
    {
      fft_histograms[l+k*histogram_size][0] = fft_data[l][0];
      fft_histograms[l+k*histogram_size][1] = fft_data[l][1];
    }
  }

  //! Compute kernel FFTW
  fftw_complex *fft_kernels = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * sigma_max*histogram_size);
  for(int sigma=1; sigma<=sigma_max; ++sigma)
  {
    const float inv_sigma   = 1.f / (2*sigma*sigma);
    const int   half_kernel = 3*sigma;

    //! Kernel computation
    float inv_cst = 0;
    for(int k=0; k<=half_kernel; ++k)
    {
      data[k]  = exp(-inv_sigma*k*k);
      if(k) inv_cst += data[k];
    }
    fill(data.begin()+half_kernel+1, data.end()-half_kernel, .0f);
    inv_cst = 1.f / (2*inv_cst + data[0]);
    for(int k=0; k<=half_kernel; ++k)
    {
      data[k] = log( (1-p)*data[k]*inv_cst + p/256 ) - log(p/256);
      if(k) data[256-k] = data[k];
    }

    //! Kernel FFT computation
    fftw_execute_dft_r2c(fft, data.data(), fft_data);
    for(size_t l=0; l<histogram_size; ++l)
    {
      fft_kernels[l+(sigma-1)*histogram_size][0] = fft_data[l][0];
      fft_kernels[l+(sigma-1)*histogram_size][1] = fft_data[l][1];
    }
  }

  const float minimum = 0;
  const float maximum = 255;
  const float alpha   = (maximum - minimum) / (histogram_size-1);

  fftw_complex *fft_convolve   = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * histogram_size);
  fftw_plan ifft;
#pragma omp critical
  ifft = fftw_plan_dft_c2r_1d(histogram_size, fft_convolve, data.data(), FFTW_ESTIMATE);
  for(size_t k=0; k<size; ++k)
  {
    for(int sigma=0; sigma<sigma_max; ++sigma)
    {
      //! Multiplication of the two FFT
      for(size_t l=0; l<histogram_size; ++l)
      {
        fft_convolve[l][0] = fft_histograms[l+k*histogram_size][0] * fft_kernels[l+sigma*histogram_size][0]
                           - fft_histograms[l+k*histogram_size][1] * fft_kernels[l+sigma*histogram_size][1];
        fft_convolve[l][1] = fft_histograms[l+k*histogram_size][1] * fft_kernels[l+sigma*histogram_size][0]
                           + fft_histograms[l+k*histogram_size][0] * fft_kernels[l+sigma*histogram_size][1];
      }

      //! iFFT computation
      fftw_execute_dft_c2r(ifft, fft_convolve, data.data());

      //! Maximum computation
      for(size_t l=0; l<histogram_size; ++l)
      {
        if( (!sigma && !l) || data[l]>max_value[k] )
        {
          max_value[k]  = data[l];
          modes->at(k)  = alpha*l+minimum;
          scale->at(k)  = sigma+1;
        }
      }
    }
  }

  //! Memory free
  fftw_destroy_plan(fft);
  fftw_destroy_plan(ifft);
  fftw_free(fft_kernels);
  fftw_free(fft_histograms);
  fftw_free(fft_convolve);

  return EXIT_SUCCESS;
}

int estimp_road(
    const vector<float> &im_noisy,
    const ImageSize &imSize,
    vector<float> *roadmask,
    float *p_estimated
    )
{
  //! Parameters
  const int    half_length = 1;
  const size_t Nroadbis    = 4;
  const float  threshold   = 70;

  //! Constant computation
  const int length = 2*half_length + 1;
  const int size   = length*length;

  //! Vector initialisation
  *roadmask = vector<float>(imSize.wh, .0f);
  *p_estimated = .0f;

  //! Extraction des patches
  vector<float> patches(size);
  const int c = 0;
  for(int h=1; h<imSize.height-1-1; ++h)
  {
    for(int w=1; w<imSize.width-1-1; ++w)
    {
      //! Compute patch
      for(int j=-half_length; j<=half_length; ++j)
      {
        for(int i=-half_length; i<=half_length; ++i)
        {
          size_t pos   = i+half_length + (j+half_length)*length;
          patches[pos] = abs(im_noisy[imSize.symmetrise_at(h+j, w+i, c)] - im_noisy[imSize.at(h,w,c)]);
        }
      }

      //! Sort vector
      sort(patches.begin(), patches.end());

      //! Estimation of the differences
      float sum = patches[1];
      for(size_t i=2 ; i<=Nroadbis; ++i)
      {
        sum += patches[i];
      }
      if(sum > threshold)
      {
        roadmask->at(imSize.at(h, w, c)) = 255.f;
        ++(*p_estimated);
      }
    }
  }

  *p_estimated /= (imSize.width-2)*(imSize.height-2);
  cout << "Estimated probability of noise: " << *p_estimated << endl;

  return EXIT_SUCCESS;
}

int estim_npatchs(const float &p)
{
  const int pos = floor(10.0001f*p);
  const int n_values[8] = {7, 9, 12, 16, 22, 35, 60, 90};

  if(pos <= 0)
  {
    return n_values[0];
  }
  else if(pos >= 8)
  {
    return n_values[7];
  }

  int value = 10*(p-0.1f*pos)*(n_values[pos]-n_values[pos-1]);
  value = value > 0 ? value : 0;
  return value + n_values[pos-1];
}

vector<float> ComputeIncompleteBetaFunction(
    const float &x,
    const int   &patch_half_length
    )
{
  const int size = (2*patch_half_length+1)*(2*patch_half_length+1);
  vector<float> coeff(size);

  //! Computation of the x^k/k!
  coeff[0] = x;
  for(int k=1; k<size; ++k)
  {
    coeff[k] = coeff[k-1]*x;
  }

  //! Binomial coefficients multiplication
  double alpha = 1;
  for(int k=size-1; k>0; --k)
  {
    alpha *= (k+1)*(1-x)/(size-k);
    coeff[k-1] *= alpha;
  }

  //! Binomial sum
  for(int k=size-1; k>0; --k)
  {
    coeff[k-1] += coeff[k];
  }

  // We return the vector of the I_x(k,n-k+1)
  return coeff;
}
