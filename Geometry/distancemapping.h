#ifndef DISTANCEMAPPING_H
#define DISTANCEMAPPING_H

#include <vector>
#include <iostream>
#include <cstdint>

namespace distmap
{
    class SaitosAlgorithm
    {
    public:
        double beta = 1.; //Aspect ratio in z-direction
        float bgthreshold = 0.;

        SaitosAlgorithm(){}
        SaitosAlgorithm(double beta_):beta(beta_){}

        void SetBackgroundThreshold(float value) {bgthreshold = value;}

    public: //Interface functions
        float* SEDM_foreground(uint8_t* input_stack, int inshape[3]);
        float* SEDM_background(uint8_t* input_stack, int inshape[3]);
        float* sedm2edm(float* sedm_stack, int shape[3]);

    private: //Implementation
        void FirstTransformation_fg(uint8_t* input_stack, float* output_stack, int shape[3]);
        void FirstTransformation_bg(uint8_t* input_stack, float* output_stack, int shape[3]);
        void SecondTransformation(float* output_stack, int shape[3]);
        void ThirdTransformation(float* output_stack, int shape[3]);
    };
}

#endif // DISTANCEMAPPING.H

