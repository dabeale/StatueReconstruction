#ifndef SEGMENTATIONMODEL_H
#define SEGMENTATIONMODEL_H

/* ----------------------------------------------------------------------
 * Copyright (C) 2016 Daniel Beale. All rights reserved.
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 * ---------------------------------------------------------------------- */

#include <vector>
#include <limits>
#include <array>
#include <map>

#include "ImageInfo.h"
#include "InfiniteGMM.h"
#include "maxflow.h"
#include "mesh.h"
#include "Matrix.h"

namespace SegmentationGUI
{

/**
 * @brief The SegmentationModelParameters class
 * A class containing the model parameters. This includes the parameters for the Markov random field
 * and also the infinite Gaussian mixture model.
 */
class SegmentationModelParameters
{
public:
    /**
     * @brief The SegType enum
     * An enumeration of Markov random field algorithms
     */
    enum SegType
    {
        GraphCut,
        ICM
    };

    /**
     * @brief SegmentationModelParameters
     * the default constructor
     */
    SegmentationModelParameters() :
        m_Alpha({100.0,100.0}), m_NSample(5), m_edgeWeight(10.0), m_terminalWeight(0.1), m_type(GraphCut) {}

    std::array<double, 2> m_Alpha;      ///< The dirichlet prior hyper parameter Alpha
    uint32_t m_NSample;                 ///< The number of samples to take when evaluating edge costs (dirichlet prior)
    double m_edgeWeight;                ///< Edge weight for graph cut
    double m_terminalWeight;            ///< Terminal weight for graph cut
    SegType m_type;                     ///< The segmentation type
};

/**
 * @brief The DepthModel struct
 * A collection of parameters for the depth model which is computed from the user annotations
 * and depth maps. It represents the a Gaussian model for the probability of depth, given the depth
 * map, and a binomial model for the probability that the foreground or background is defined.
 *
 * @todo This depth model works better when there is more noise in the point cloud. For point clouds
 * which exhibit next to no noise, the covariances are so small that it cannot infer what
 * is behind the object in question. The solution would be to compute the covariances in the
 * user strokes on each image.
 */
struct DepthModel
{
    typedef std::array<double, 2> MeanVar;              ///< A mean and variance type
    typedef std::array<double, 2> BinaryProbability;    ///< A binary probability type

    MeanVar m_foregroundDepth;                          ///< The Gaussian model for the foreground depths
    MeanVar m_backgroundDepth;                          ///< The Gaussian models for the background depths

    BinaryProbability m_foregroundUndefinedProbability; ///< The probability that the foreground is undefined
    BinaryProbability m_backgroundUndefinedProbability; ///< The probability that the background is undefined
};


/**
 * @brief The SegmentationModel class
 *
 * The segmentation model combines depth and colour information in order to disambiguate the foreground and
 * background. The parameters of the model are learned from the annotations that are drawn on a few images
 * by the user, allowing the algorithm to infer foreground and background in unlabelled images.
 *
 * There are 3 collections of parameters in the model,
 *  1. The parameters for the depth model is broken in to two,
 *     - A mean and variance for the foreground and background depth,
 *     - The binomial probability that the foreground and background is underfined
 *  2. A Gaussian mixture model representing the probability that a pixel takes on a RGB value.
 *     This model is trained using a Dirichlet process model, with a Monte Carlo algorithm to estimate the parameters.
 *
 * The graphical model for a single pixel in the segmentation algorithm is as follows,
 *
 *         ( depth model )       ( colour model )
 *                    \           /
 *                     \         /
 *                   ( ij_th pixel class )
 * i.e. the depth model and colour model are independent.
 *
 * The pixels are arranged in to a Markov random field with 4 neighbour connectivity.
 */
class SegmentationModel
{
public:
    /**
     * @brief SegmentationModel
     * The explicit constructor for the Segmentation model, which takes the collection of images with
     * their information, the 3D point cloud and a collection of camera matrices.
     * @param iminfos The images
     * @param mesh The 3D point cloud
     * @param ProjectionMatrices The projection matrices, in the same order as for iminfos.
     */
    SegmentationModel(std::vector< ImageInfo >& iminfos, const Mesh& mesh,
                      const std::vector<Math::Matrix>& ProjectionMatrices);

    /**
     * @brief CreateData
     * Create a matrix of foreground colours and a matrix of background colours by uniformly sampling
     * the user annotations. the number of samples is governed by the member variable m_NSamples.
     */
    void CreateData(  );

    /**
     * @brief RunSegmentation
     * Compute all of the parameters for the model given the user input, and then infer the classes of each pixel
     * given the parameters.
     * @param variance The variance hyperparameter for the IGMM
     * @param alpha The Alpha Dirichlet parameter for the IGMM
     */
    void RunSegmentation( double variance, double alpha);

    /**
     * @brief ComputeDepthModels
     * Compute the depth models for the point cloud given the calibrated views.
     */
    void ComputeDepthModels();

    /**
     * @brief ComputePointProbabilities
     * Compute the probabilities that each point in the point cloud belongs to the fg / bg, given the
     * probabilities in the image space.
     */
    void ComputePointProbabilities();

    /**
     * @brief EvaluateMRF
     * Evaluate the class of a pixel given the model parameters and the Markov random field assumption.
     * @param msg The interface for printing messages to stdout.
     */
    void EvaluateMRF(Stream::Message &msg);

    /**
     * @brief GetParameters
     * Return the chosen parameters for the model
     * @return
     */
    inline SegmentationModelParameters& GetParameters(){return m_params;}
private:
    std::vector< ImageInfo >& m_imageinfos;                 ///< The images
    const Mesh& m_3dmesh;                                   ///< The 3D mesh
    const std::vector<Math::Matrix>& m_ProjectionMatrices;  ///< The projection matrices

    Math::Matrix m_ForegroundColours; ///< A matrix of foreground points
    Math::Matrix m_BackgroundColours; ///< A matrix of selected background colours (this matrix is sampled)

    std::map<int, DepthModel> m_depthModel;                         ///< The depth models for each camera that has been draw on
    std::vector<DepthModel::BinaryProbability> m_pointProbability;  ///< The probability that each point belongs to the foreground / background repectively.
    DepthModel::BinaryProbability m_fgprobabilityUndefined;         ///< The probabilities that a pixel does not have a depth value associated with it
    DepthModel::BinaryProbability m_bgprobabilityUndefined;         ///< The probabilities that a pixel does not have a depth value associated with it


    uint32_t m_NSamples;                    ///< The number of colour samples to take
    Stream::Message m_msg;                  ///< The message streamer
    bool m_depthModelsComputed;             ///< True when the depth models have been computed
    std::vector<Math::Matrix> m_fprob;      ///< the foreground probabilities for each image
    std::vector<Math::Matrix> m_bprob;      ///< the background probabilities for each image

    SegmentationModelParameters m_params;   ///< The model parameters
};
}
#endif
