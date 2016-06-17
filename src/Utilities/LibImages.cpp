/*
 * Copyright (c) 2013, Marc Lebrun <marc.lebrun.ik@gmail.com>
 * Copyright (c) 2015, Thierry Guillemot <thierry.guillemot.work@gmail.com>
 * All rights reserved.
 *
 * This program is free software: you can use, modify and/or
 * redistribute it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later
 * version. You should have received a copy of this license along
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file LibImages.cpp
 * @brief Usefull functions on images
 *
 * @author Marc Lebrun <marc.lebrun.ik@gmail.com>
 **/

#include "Utilities/LibImages.h"

#include <iostream>
#include <sstream>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>

#include "io_png.h"
#include "mt19937ar.h"
#include "utils.h"

using namespace std;

int loadImage(
    const char*    p_name,
    vector<float> &o_im,
    ImageSize     &o_imSize,
    const bool     p_verbose
    )
{
  if (p_verbose)
  {
    cout << endl << "Read input image '" << p_name << "'...";
  }

  size_t w, h, c;
  float *imTmp = read_png_f32(p_name, &w, &h, &c);
  if (!imTmp)
  {
    cout << "error :: " << p_name << " not found or not a correct png image" << endl;
    return EXIT_FAILURE;
  }

  if (p_verbose)
  {
    cout << "done." << endl;
  }

  //! test if image is really a color image and exclude the alpha channel
  if (c == 2)
  {
    cout << "The gray image has a alpha channel. We will remove it." << endl;
    c = 1;
  }
  else if(c >= 3)
  {
    cout << "Parigi does not work with color images." << endl;
    exit(-1);
  }

  if (p_verbose)
  {
    cout << "image size :" << endl;
    cout << " - width          = " << w << endl;
    cout << " - height         = " << h << endl;
    cout << " - nb of channels = " << c << endl;
  }

  o_imSize.width      = w;
  o_imSize.height     = h;
  o_imSize.nChannels  = c;
  o_imSize.wh         = w * h;
  o_imSize.whc        = w * h * c;
  o_im.resize(w * h * c);
  for(size_t _c=0; _c<c; ++_c)
    {
        for(size_t _h=0; _h<h; ++_h)
        {
            for(size_t _w=0; _w<w; ++_w)
            {
              o_im[o_imSize.at(_h,_w,_c)] = imTmp[_c*w*h + _h*w + _w];
            }
        }
    }

  delete[] imTmp;

  return EXIT_SUCCESS;
}

int saveImage(
    const char          *p_name,
    const vector<float> &i_im,
    const ImageSize     &p_imSize,
    const float         &p_min,
    const float         &p_max
    )
{
  float* imTmp = new float[p_imSize.whc];

  for(int _c=0; _c<p_imSize.nChannels; ++_c)
  {
    for(int _h=0; _h<p_imSize.height; ++_h)
    {
      for(int _w=0; _w<p_imSize.width; ++_w)
      {
        imTmp[_c*p_imSize.wh + _h*p_imSize.width + _w] = i_im[p_imSize.at(_h,_w,_c) ];
      }
    }
  }

  //! Check for boundary problems
  for (int k = 0; k < p_imSize.whc; k++)
  {
    imTmp[k] = imTmp[k] < p_min ? p_min : (imTmp[k] > p_max ? p_max : imTmp[k]);
  }

  if (write_png_f32(p_name, imTmp, p_imSize.width, p_imSize.height, p_imSize.nChannels) != 0)
  {
    cout << "... failed to save png image :'" << p_name << "'" << endl;
    return EXIT_FAILURE;
  }

  delete[] imTmp;

  return EXIT_SUCCESS;
}

void addImpulseNoise(
    const vector<float> &i_im,
    vector<float>       &o_imNoisy,
    const float          p,
    const bool           p_verbose
    )
{
  if (p_verbose)
  {
    cout << "Add noise [probability = " << p << "] ...";
  }

  o_imNoisy = i_im;
  mt_init_genrand((unsigned long int) time (NULL) + (unsigned long int) getpid());

  for (unsigned k = 0; k < i_im.size(); ++k)
  {
    const double noise = mt_genrand_res53();
    const double value = mt_genrand_res53();
    o_imNoisy[k] = noise>p ? i_im[k] : 256*value;
  }

  if (p_verbose)
  {
    cout << "done." << endl;
  }
}

void addGaussianNoise(
    const vector<float> &i_im,
    vector<float>       &o_imNoisy,
    const float               p_sigma,
    const bool                p_verbose
    )
{
  if (p_verbose)
  {
    cout << "Add noise [sigma = " << p_sigma << "] ...";
  }

  o_imNoisy = i_im;
  mt_init_genrand((unsigned long int) time (NULL) + (unsigned long int) getpid());

  for (unsigned k = 0; k < i_im.size(); ++k)
  {
    const double a = mt_genrand_res53();
    const double b = mt_genrand_res53();

    o_imNoisy[k] += p_sigma * (float) (sqrtl(-2.0l * log(a)) * cos(2.0l * M_PI * b));
    o_imNoisy[k] = Crop(o_imNoisy[k], 0.f, 255.f);
  }

  if (p_verbose)
  {
    cout << "done." << endl;
  }
}

int computePsnr(
    const vector<float> &i_im1,
    const vector<float> &i_im2,
    float               &o_psnr,
    float               &o_rmse,
    const char*          p_imageName,
    const bool           p_verbose
    )
{
  if (i_im1.size() != i_im2.size())
  {
    cout << "Can't compute PSNR & RMSE: images have different sizes: " << endl;
    cout << "i_im1 : " << i_im1.size() << endl;
    cout << "i_im2 : " << i_im2.size() << endl;
    return EXIT_FAILURE;
  }

  float sum = 0.f;
  for (unsigned k = 0; k < i_im1.size(); k++)
    sum += (i_im1[k] - i_im2[k]) * (i_im1[k] - i_im2[k]);

  o_rmse = sqrtf(sum / (float) i_im1.size());
  o_psnr = 20.f * log10f(255.f / o_rmse);

  if (p_verbose)
  {
    cout << p_imageName << endl;
    cout << "PSNR = " << o_psnr << endl;
    cout << "RMSE = " << o_rmse << endl;
  }

  return EXIT_SUCCESS;
}

vector<float> ComputeNormalizedCumulativeHistogram(
    const float  *data_ptr,
    const size_t  size,
    const float   min,
    const float   max,
    const size_t  histogram_size
    )
{
  vector<float> histogram(histogram_size, 0);
  size_t size_max = histogram_size - 1;

  //! Computation of the histogram
  for(size_t pos=0; pos<size; ++pos)
  {
    ++histogram[(size_t) ( size_max*(data_ptr[pos]- min)/(max - min))];
  }

  //! Normalization of the cumulative histogram
  for(size_t pos=0; pos<histogram_size-1; ++pos)
  {
    histogram[pos+1] += histogram[pos];
    histogram[pos]   /= size;
  }
  histogram[histogram_size-1] /= size;

  return histogram;
}

vector<float> ComputeNormalizedHistogram(
    const float  *data_ptr,
    const size_t  size,
    const float   min,
    const float   max,
    const size_t  histogram_size
    )
{
  vector<float> histogram(histogram_size, 0);
  size_t size_max = histogram_size - 1;

  //! Computation of the histogram
  for(size_t pos=0; pos<size; ++pos)
  {
    ++histogram[(size_t) ( size_max*(data_ptr[pos]- min)/(max - min))];
  }

  //! Normalization of the cumulative histogram
  for(size_t pos=0; pos<histogram_size; ++pos)
  {
    histogram[pos]   /= size;
  }

  return histogram;
}

vector<float> ExtractPatches(
    const vector<float> &im_input,
    const ImageSize     &imSize,
    const int            length
    )
{
  const int dim_patch    = length*length*imSize.nChannels;
  const int im_height    = imSize.height-length;
  const int im_width     = imSize.width-length;
  const int number_patch = im_height*im_width;

  vector<float> patches(number_patch*dim_patch);

  for(int c=0; c<imSize.nChannels; ++c)
  {
    for(int w=0; w<im_width; ++w)
    {
      for(int h=0; h<im_height; ++h)
      {
        for(int j=0; j<length; ++j)
        {
          for(int i=0; i<length; ++i)
          {
            size_t pos_patch = i + j*length + c*length*length;
            size_t pos_pixel = h + w*im_height;
            patches[pos_pixel*dim_patch + pos_patch] = im_input[imSize.at(h+i, w+j, c)];
          }
        }
      }
    }
  }

  return patches;
}

vector<float> ExtractCenteredPatches(
    const vector<float> &im_input,
    const ImageSize     &imSize,
    const int            l
    )
{
  const int length       = 2*l+1;
  const int dim_patch    = length*length;
  const int im_height    = imSize.height-2*l;
  const int im_width     = imSize.width-2*l;
  const int im_wh        = im_height*im_width;
  const int number_patch = im_height*im_width*imSize.nChannels;
  vector<float> patches(number_patch*dim_patch);

  int c=0;
  for(int w=l; w<imSize.width-l; ++w)
  {
    for(int h=l; h<imSize.height-l; ++h)
    {
      for(int j=-l; j<=l; ++j)
      {
        for(int i=-l; i<=l; ++i)
        {
          size_t pos_patch = i+l + (j+l)*length;
          size_t pos_pixel = (h-l) + (w-l)*im_height + c*im_wh;
          patches[pos_pixel*dim_patch + pos_patch] = im_input[imSize.at(h+i, w+j, c)];
        }
      }
    }
  }

  return patches;
}

int ExtractSubImageWithBoundary(
    const vector<float> &im_input,
    const ImageSize     &imSize_input,
    const ImageSize     &imSize_output,
    const size_t        &h_offset,
    const size_t        &w_offset,
    const size_t        &boundary_length,
    vector<float>       *im_output
    )
{
  im_output->resize(imSize_output.whc);

  for(int c=0; c<imSize_output.nChannels; ++c)
  {
    for(int w=0; w<imSize_output.width; ++w)
    {
      for(int h=0; h<imSize_output.height; ++h)
      {
        im_output->at(imSize_output.at(h,w,c)) = im_input[imSize_input.symmetrise_at(h-boundary_length+h_offset, w-boundary_length+w_offset, c)];
      }
    }
  }

  return EXIT_SUCCESS;
}

int EstimateSubImagesDivision(
    const size_t &number_sub_images,
    size_t *a,
    size_t *b
                          )
{
  if(number_sub_images == 1)
  {
    *a = 1;
    *b = 1;
    return EXIT_SUCCESS;
  }
  else if(number_sub_images%2) //! Check if the number ispower of 2
  {
    cerr << "The number of division must be a power of 2" << endl;
    return EXIT_FAILURE;
  }

  *b = 2;
  while (number_sub_images % *b > 0)
  {
    ++(*b);
  }
  *a = number_sub_images / *b;

  if (*b > *a) {
    *a = *b;
    *b = number_sub_images / *a;
  }

  return EXIT_SUCCESS;
}

int SubdiveImage(
    const vector<float>    &im_input,
    const ImageSize        &imSize_input,
    const size_t           &number_sub_images,
    const size_t           &boundary_length,
    vector<vector<float> > *im_subdivide,
    ImageSize              *imSize_subdivide
    )
{
  //! Estimate the number of division on each direction
  size_t h_div, w_div;
  EstimateSubImagesDivision(number_sub_images, &h_div, &w_div);
  //! We ensure that the biggest division is done over the maximum size
  im_subdivide->resize(h_div*w_div);
  if(imSize_input.width > imSize_input.height)
  {
    size_t swap_value = h_div;
    h_div = w_div;
    w_div = swap_value;
  }
  const unsigned hTmp = ceil(float(imSize_input.height) / float(h_div));
  const unsigned wTmp = ceil(float(imSize_input.width)  / float(w_div));
  //! Compute the size of the sub images
  imSize_subdivide->height    = hTmp + 2* boundary_length;
  imSize_subdivide->width     = wTmp + 2* boundary_length;
  imSize_subdivide->nChannels = imSize_input.nChannels;
  imSize_subdivide->wh        = imSize_subdivide->height * imSize_subdivide->width;
  imSize_subdivide->whc       = imSize_subdivide->wh * imSize_subdivide->nChannels;
  //! Computation of the subimages with the boundaries
  for(size_t y=0; y<w_div; ++y)
  {
    for(size_t x=0; x<h_div; ++x)
    {
      ExtractSubImageWithBoundary(im_input, imSize_input, *imSize_subdivide, x*hTmp, y*wTmp, boundary_length, &((*im_subdivide)[x+y*h_div]));
    }
  }

  return EXIT_SUCCESS;
}

int PartialRecomposeImageWithBoundaries(
    const vector<float> &im_subdivide,
    const ImageSize     &imSize_subdivide,
    const ImageSize     &imSize_output,
    const int           &h_offset,
    const int           &w_offset,
    const int           &boundary_length,
    const int           &boundary_w_length,
    const int           &boundary_h_length,
    vector<float>       *im_output
    )
{
  //! Compute the image
  for(int c=0; c<imSize_output.nChannels; ++c)
  {
    for(int w=0; w<imSize_subdivide.width-boundary_w_length; ++w) //X
    {
      if(w+w_offset >= imSize_output.width) continue;
      for(int h=0; h<imSize_subdivide.height-boundary_h_length; ++h) //X
      {
        if(h+h_offset >= imSize_output.height) continue;
        im_output->at(imSize_output.at(h+h_offset, w+w_offset, c)) = im_subdivide[imSize_subdivide.at(h+boundary_length, w+boundary_length, c)];
      }
    }
  }

  return EXIT_SUCCESS;
}

vector<float> RecomposeImageWithBoundaries(
    const vector<vector<float> > &im_subdivide,
    const ImageSize              &imSize_input,
    const size_t                 &number_sub_images,
    const size_t                 &boundary_length
    )
{
  vector<float> recomposed_image(imSize_input.whc);

  //! Estimate the number of division on each direction
  size_t h_div, w_div;
  EstimateSubImagesDivision(number_sub_images, &h_div, &w_div);
  //! We ensure that the biggest division is done over the maximum size
  if(imSize_input.width > imSize_input.height)
  {
    size_t swap_value = h_div;
    h_div = w_div;
    w_div = swap_value;
  }

  const unsigned hTmp = ceil(float(imSize_input.height) / float(h_div));
  const unsigned wTmp = ceil(float(imSize_input.width)  / float(w_div));

  //! Compute the size of the sub images
  ImageSize imSize_subdivide;
  imSize_subdivide.height    = hTmp + 2* boundary_length;
  imSize_subdivide.width     = wTmp + 2* boundary_length;
  imSize_subdivide.nChannels = imSize_input.nChannels;
  imSize_subdivide.wh        = imSize_subdivide.height * imSize_subdivide.width;
  imSize_subdivide.whc       = imSize_subdivide.wh * imSize_subdivide.nChannels;
  //! Computation of the subimages with the boundaries
  for(size_t y=0; y<w_div; ++y)
  {
    int boundary_w_length = (imSize_input.width%2) && (y==w_div-1) ? 2*boundary_length + 1 : 2*boundary_length;
    for(size_t x=0; x<h_div; ++x)
    {
      int boundary_h_length = (imSize_input.height%2) && (x==h_div-1) ? 2*boundary_length + 1 : 2*boundary_length;
      PartialRecomposeImageWithBoundaries(im_subdivide[x+y*h_div], imSize_subdivide, imSize_input, x*hTmp, y*wTmp, boundary_length, boundary_w_length, boundary_h_length, &recomposed_image);
    }
  }

  return recomposed_image;
}

vector<float> RecomposeImage(
    const vector<vector<float> > &im_subdivide,
    const ImageSize              &imSize_input,
    const size_t                 &number_sub_images
    )
{
  return RecomposeImageWithBoundaries(im_subdivide, imSize_input, number_sub_images, 0);
}
