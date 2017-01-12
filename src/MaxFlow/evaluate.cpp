#include "maxflow.h"

namespace MaxFlow {
    Math::Matrix EvaluateFGBGFromProbabilities( const Math::Matrix& fgprobs, const Math::Matrix& bgprobs, HyperParameters& hypes )
    {
        hypes.msg.print(2) << "Running graph cut" << std::endl;
        const uint32_t M = fgprobs.Rows();
        const uint32_t N = fgprobs.Cols();

        // create graph
        hypes.msg.print(2) << "Filling graph" << std::endl;
        typedef Graph<double,double,double> GraphType;
        GraphType *g = new GraphType( M*N, 4*N*M);
        g->add_node(M*N);

        for( uint32_t i=0; i<M;++i)
        for( uint32_t j=0; j<N;++j)
        {
            uint32_t index = j*M + i;
            Arr2d expected0 {hypes.N*fgprobs(i,j), hypes.N*bgprobs(i,j)};
            double sum = hypes.N + hypes.Alpha[0] + hypes.Alpha[1];

            if( i < M - 1 )
            {
                Arr2d expected {hypes.N*fgprobs(i+1,j), hypes.N*bgprobs(i+1,j)};
                Arr2d predicted { (expected[0] + hypes.Alpha[0])/sum, (expected[1] + hypes.Alpha[1])/sum };
                double cprob = hypes.edgeWeight*std::exp(LogMultinomial( expected0, predicted ));
                g->add_edge(index, index+1, cprob , 0.0f);
            }

            if( i > 0 )
            {
                Arr2d expected {hypes.N*fgprobs(i-1,j), hypes.N*bgprobs(i-1,j)};
                Arr2d predicted { (expected[0] + hypes.Alpha[0])/sum, (expected[1] + hypes.Alpha[1])/sum };
                double cprob = hypes.edgeWeight*std::exp(LogMultinomial( expected0, predicted ));
                g->add_edge(index, index-1, cprob , 0.0f);
            }

            if( j < N - 1 )
            {
                Arr2d expected {hypes.N*fgprobs(i,j+1), hypes.N*bgprobs(i,j+1)};
                Arr2d predicted { (expected[0] + hypes.Alpha[0])/sum, (expected[1] + hypes.Alpha[1])/sum };
                double cprob = hypes.edgeWeight*std::exp(LogMultinomial( expected0, predicted ));
                g->add_edge(index, index+M, cprob , 0.0f);
            }

            if( j > 0 )
            {
                Arr2d expected {hypes.N*fgprobs(i,j-1), hypes.N*bgprobs(i,j-1)};
                Arr2d predicted { (expected[0] + hypes.Alpha[0])/sum, (expected[1] + hypes.Alpha[1])/sum };
                double cprob = hypes.edgeWeight*std::exp(LogMultinomial( expected0, predicted ));
                g->add_edge(index, index-M, cprob , 0.0f);
            }

            g->add_tweights(index, hypes.terminalWeight*(fgprobs(i,j)), 0.0);
            g->add_tweights(index, 0.0, hypes.terminalWeight*(bgprobs(i,j)));
        }

        hypes.msg.print(2) << "Computing flow" << std::endl;
        double totalflow = g->maxflow();
        hypes.msg.print(2) << "Done. The total flow is : " << totalflow << std::endl;

        hypes.msg.print(1) << "Labelling" << std::endl;
        Math::Matrix labels(M,N);
        for( uint32_t i=0; i<M;++i)
        for( uint32_t j=0; j<N;++j)
        {
            uint32_t index = j*M + i;
            labels(i,j) = static_cast<double>(g->what_segment(index));
        }

        return labels;
    }
}
