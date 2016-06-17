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

#ifndef PARIGI_PARIGI_H
#define PARIGI_PARIGI_H

#include "Utilities/LibImages.h"

/**
 * @brief This function generates the noisy image from an input image
 *        and apply the parigi denoising with one/two iterations and
 *        apply/not apply the mixed criterium.
 *
 * @param input_filename: input image filename
 * @param noisy_filename: output noisy image filename
 * @param denoised_filename: output denoised image filename
 * @param p_gaussian_noise: sigma value of the gaussian noise added
 *        to the original image
 * @param p_impulse_noise: probability of the impulse noise added
 *        to the original image
 * @param two_iteration: 'true' for two iterations of parigi,
                         'false' for one iteration
 * @param mixed_image: 'true' to mix the noised and the denoised images
 *                     'false' otherwise
 * @return: EXIT_SUCCESS if the code is executed without problem
            EXIT_FAILURE otherwise
 */
int LaunchParigi(
    const char  *input_filename,
    const char  *noisy_filename,
    const char  *denoised_filename,
    const float &p_gaussian_noise,
    const float &p_impulse_noise,
    const bool  &two_iteration,
    const bool  &mixed_image,
    const bool  &verbose
    );

/**
 * @brief Restore an image corrupted by an impulse and a
 *        gaussian noise.
 *
 * @param im_noisy_bound: input noisy image.
 * @param imSize_bound: size of the noisy image
 * @param half_size_patch: patch half size
 * @param half_size_neighborhood: neighborhood half size
 * @param p_estimated: impulse noise probability
 * @param number_patch: patches number needed by the restoration
 * @param coeff: binomial coefficients need by the restoration
 *
 * @return mode: pointer to the resulting denoised image
 * @return sigma: pointer to the resulting sigma value of the mode
 */

int RestoreImpulseNoise(
    const std::vector<float> &im_noisy_bound,
    const ImageSize          &imSize_bound,
    const float              &half_size_patch,
    const float              &half_size_neighborhood,
    const float              &p_estimated,
    const int                &number_patch,
    const std::vector<float> &coeff,
    std::vector<float>       *modes,
    std::vector<float>       *sigma
    );

/**
 * @brief Compute the mixed restore image from the noisy image, the
 *        denoised image and the sigma value of the modes.
 *
 * @param im_noisy: noisy image
 * @param modes: denoised image
 * @param sigma: sigma values of the modes
 *
 * @return mixed denoised image
 **/
int ComputeMixedCriterium(
    const std::vector<float> &im_noisy,
    const std::vector<float> &modes,
    const std::vector<float> &sigma,
    std::vector<float>       *im_denoised);

/**
 * @brief Compute the number_patch nearest neighbors and return a
 * vector of the pointer of the patches and their distances in the
 * ascending order.
 *
 * @param coeff: binomial coefficients
 * @param patches: vector containing all the extracted patches of the
 *        noisy image
 * @param imSize_neigh: size of the noisy image
 * @param imSize_neigh: size of the noisy image without the
 *        neighborhood boundary.
 * @param dim_patch: patch dimensionality
 * @param half_size_neighborhood: neighborhood half size
 * @param number_patch: number_patch needed by the restoration
 *
 * @return vector containing a pointer to the number_patch nearest
 *         patches and their distances
 **/
std::vector<std::pair<float, const float*> > ComputeNearestPatches(
    const std::vector<float> &coeff,
    const std::vector<float> &patches,
    const ImageSize          &imSize_neigh,
    const ImageSize          &imSize,
    const int                &dim_patch,
    const int                &half_size_neighborhood,
    const size_t             &number_patch
    );

/**
 * @brief Compute the histograms from the nearest neighbors
 *
 * @param dist_patch: pair of the nearest patches and their distances
 * @param imSize: image size
 * @param half_lenght_patch: patch half length
 * @param size_histogram: histogram dimension
 * @param number_patch: number of patches used by the restoration
 * @param dim_neigh: dimension of the neighborhood
 *
 * @return vector containing all the compute histograms
 **/
std::vector<float> ComputeHistograms(
    const std::vector<std::pair<float, const float*> > &dist_patch,
    const ImageSize                                    &imSize,
    const int                                          &half_length_patch,
    const size_t                                       &size_histogram,
    const size_t                                       &number_patch,
    const size_t                                       &dim_neigh
    );

/**
 * @brief Compute the best mode and scale for each histogram
 *
 * @params histograms: vector containing all histograms extracted from
 *         the noisy image
 * @param size: number of histograms contained in the vector histograms
 * @param histogram_size: dimensionality of the histogram
 * @param p: impulse noise probabilty
 * @param sigma_max: max sigma value used for the convolution
 *
 * @return modes: pointer to the vector containing the mode values of
 *                each histograms
 * @return scale: pointer to the vector containing the scale values of
 *                each histogram modes
 **/
int ComputeModeFFTW(
    const std::vector<float> &histograms,
    const size_t             &size,
    const size_t             &histogram_size,
    const float              &p,
    const int                &sigma_max,
    std::vector<float>       *modes,
    std::vector<float>       *scale
    );

/**
 * @brief Estimate the probability and the mask of the corrupted pixel
 * of the image
 *
 * @param im_noisy: noisy image
 * @param imSize: size of the noisy image
 *
 * @return roadmask: pointer to the mask image of the corrupted pixels
 * @return p_estimage: pointer to the probability of the impulse noise
 */
int estimp_road(
    const std::vector<float> &im_noisy,
    const ImageSize          &imSize,
    std::vector<float>       *roadmask,
    float                    *p_estimated
    );

/**
 * @brief Return the number of patches needed by the parigi
 * restoration. These values are the one we have defined in our IPOL
 * paper.
 *
 * @param p: probability of the impulse noise.
 *
 * @return Number of patches needed by the restoration.
 **/
int estim_npatchs(const float &p);

/**
 * @brief Compute the binomial needed by the parigi restoration
 *
 * @param x: probability of the impulse noise
 * @param patch_half_length: patch half length used by the restoration
 *
 * @return lookup table of the binomial coefficient needed by the
 * parigi restoration
 */
std::vector<float> ComputeIncompleteBetaFunction(
    const float &x,
    const int   &patch_half_length
    );

#endif // PARIGI_PARIGI_H
