/*
 * Copyright (c) 2015, Thierry Guillemot <thierry.guillemot.work@gmail.com>
 * All rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Parigi/parigi.h"
#include "Utilities/LibImages.h"

using namespace std;

/**
 * @brief print usage function
 */
void PrintUsage()
{
  cerr << "usage : image_input image_noisy image_denoised "
       << "p_gaussian_noise p_impulse_noise "
       << "two_step mixed_image verbose" << endl;
}

/**
 * @brief main function call
 */
int main(int argc, char *argv[])
{
  if(argc!=9)
  {
    PrintUsage();
    return EXIT_FAILURE;
  }

  const char* image_input       = argv[1];
  const char* image_noisy       = argv[2];
  const char* image_denoised    = argv[3];
  const float p_gaussian_noise  = atof(argv[4]);
  const float p_impulse_noise = atof(argv[5]);
  const bool  two_step          = static_cast<bool>(atof(argv[6]));
  const bool  mixed_image       = static_cast<bool>(atof(argv[7]));
  const bool  verbose           = static_cast<bool>(atof(argv[8]));
  
  LaunchParigi(image_input, image_noisy, image_denoised,
               p_gaussian_noise, p_impulse_noise,
               two_step, mixed_image, verbose);
  
  return EXIT_SUCCESS;
}
