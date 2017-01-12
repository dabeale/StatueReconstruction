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

class SegmentationModelParameters
{
public:
    enum SegType
    {
        GraphCut,
        ICM
    };

    SegmentationModelParameters() :
        m_Alpha({100.0,100.0}), m_NSample(5), m_edgeWeight(10.0), m_terminalWeight(0.1), m_type(GraphCut) {}

    std::array<double, 2> m_Alpha; ///< The dirichlet prior
    uint32_t m_NSample; ///< The number of samples to take when evaluating edge costs (dirichlet prior)
    double m_edgeWeight; ///< Edge weight for graph cut
    double m_terminalWeight; ///< Terminal weight for graph cut
    SegType m_type;
};

struct DepthModel
{
    typedef std::array<double, 2> MeanVar;
    typedef std::array<double, 2> BinaryProbability;

    MeanVar m_foregroundDepth; ///< The Gaussian model for the foreground depths
    MeanVar m_backgroundDepth; ///< The Gaussian models for the background depths

    BinaryProbability m_foregroundUndefinedProbability; ///< The probability that the foreground is undefined
    BinaryProbability m_backgroundUndefinedProbability; ///< The probability that the background is undefined
}; ///< The parameters for the Gaussian depth model

class SegmentationModel
{
public:
    SegmentationModel(std::vector< ImageInfo >& iminfos, const Mesh& mesh,
                      const std::vector<Math::Matrix>& ProjectionMatrices);

    void CreateData(  ); ///< Sample data from the images

    void RunSegmentation( double variance, double alpha); ///< Run the segmentation on an QImage

    void ComputeDepthModels(); ///< Compute the depth models from a single viewpoint

    void ComputePointProbabilities(); ///< Compute the probabilities that each point in the point cloud belongs to the fg / bg

    void EvaluateMRF(Stream::Message &msg);

    inline SegmentationModelParameters& GetParameters(){return m_params;}
private:
    std::vector< ImageInfo >& m_imageinfos; ///< The image infos
    const Mesh& m_3dmesh;
    const std::vector<Math::Matrix>& m_ProjectionMatrices;

    Math::Matrix m_ForegroundColours; ///< A matrix of foreground points
    Math::Matrix m_BackgroundColours; ///< A matrix of selected background colours (this matrix is sampled)

    std::map<int, DepthModel> m_depthModel; ///< The depth models for each camera that has been draw on
    std::vector<DepthModel::BinaryProbability> m_pointProbability; ///< The probability that each point belongs to the foreground / background repectively.
    DepthModel::BinaryProbability m_fgprobabilityUndefined; ///< The probabilities that a pixel does not have a depth value associated with it
    DepthModel::BinaryProbability m_bgprobabilityUndefined; ///< The probabilities that a pixel does not have a depth value associated with it


    uint32_t m_NSamples; ///< The number of colour samples to take

    Stream::Message m_msg; ///< The message streamer

    bool m_depthModelsComputed; ///< True when the depth models have been computed

    std::vector<Math::Matrix> m_fprob; ///< the foreground probabilities for each image
    std::vector<Math::Matrix> m_bprob; ///< the background probabilities for each image

    SegmentationModelParameters m_params;
};

#endif
