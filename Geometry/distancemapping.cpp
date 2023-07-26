#include "distancemapping.h"
#include <cmath>
#include <omp.h>

namespace distmap
{

/********************** Interface *****************************/
float* SaitosAlgorithm::SEDM_foreground(uint8_t* input_stack, int inshape[3])
{
    long long int nslice = inshape[0]*inshape[1];
    long long int nstack = inshape[2]*nslice;
    float* sedm_stack = (float*) calloc(nstack, sizeof(*sedm_stack));
    FirstTransformation_fg(input_stack, sedm_stack, inshape);
    SecondTransformation(sedm_stack, inshape);
    ThirdTransformation(sedm_stack, inshape);
    return sedm_stack;
}
float* SaitosAlgorithm::SEDM_background(uint8_t* input_stack, int inshape[3])
{
    long long int nslice = inshape[0]*inshape[1];
    long long int nstack = inshape[2]*nslice;
    float* sedm_stack = (float*) calloc(nstack, sizeof(*sedm_stack));
    FirstTransformation_bg(input_stack, sedm_stack, inshape);
    SecondTransformation(sedm_stack, inshape);
    ThirdTransformation(sedm_stack, inshape);
    return sedm_stack;
}
float* SaitosAlgorithm::sedm2edm(float *sedm_stack, int shape[3])
{
	long long int nslice = shape[0]*shape[1];
	long long int nstack = shape[2]*nslice;

	float* edm_stack = (float*) calloc(nstack, sizeof(*sedm_stack));

	#pragma omp parallel for
    for(uint64_t pos = 0; pos < nstack; pos++)
        edm_stack[pos] = sqrt((double) sedm_stack[pos]);
    return edm_stack;
}

/********************** Implementation *****************************/
void SaitosAlgorithm::FirstTransformation_fg(uint8_t* input_stack, float* output_stack, int shape[3])
{
    long long int nslice = shape[0]*shape[1];
    long long int nstack = shape[2]*nslice;

    uint64_t n_lines = nstack/shape[0]; //number of lines in first dimension to be processed

    #pragma omp parallel
    {
        #pragma omp for
        for (uint64_t n = 0; n < n_lines; n++)
        {
            //Forward scan
            uint64_t df = shape[0]-1;
            for(uint64_t xpos = n*shape[0]; xpos < (n+1)*shape[0]; xpos++)
            {
                if(input_stack[xpos] > bgthreshold)
                    df++;
                else
                    df = 0;

                output_stack[xpos] = df*df;
            }

            //Backward scan
            uint64_t db = shape[0]-1;
            for(uint64_t xpos = ((n+1)*shape[0]-1); xpos >= n*shape[0]; xpos--)
            {


                if(input_stack[xpos] != 0)
                    db++;
                else
                    db = 0;

                output_stack[xpos] = std::min(output_stack[xpos], (float) db*db);

                if (xpos == 0)
                    break;
            }
        }
    }
    return;
}
void SaitosAlgorithm::FirstTransformation_bg(uint8_t* input_stack, float* output_stack, int shape[3])
{
    long long int nslice = shape[0]*shape[1];
    long long int nstack = shape[2]*nslice;

    uint64_t n_lines = nstack/shape[0]; //number of lines in first dimension to be processed

    #pragma omp parallel
    {
        #pragma omp for
        for (uint64_t n = 0; n < n_lines; n++)
        {
            //Forward scan
            uint64_t df = shape[0]-1;
            for(uint64_t xpos = n*shape[0]; xpos < (n+1)*shape[0]; xpos++)
            {
                if(input_stack[xpos] <= bgthreshold)
                    df++;
                else
                    df = 0;

                output_stack[xpos] = df*df;
            }

            //Backward scan
            uint64_t db = shape[0]-1;
            for(uint64_t xpos = ((n+1)*shape[0]-1); xpos >= n*shape[0]; xpos--)
            {


                if(input_stack[xpos] <= bgthreshold)
                    db++;
                else
                    db = 0;

                output_stack[xpos] = std::min(output_stack[xpos], (float) db*db);

                if (xpos == 0)
                    break;
            }
        }
    }
    return;
}
void SaitosAlgorithm::SecondTransformation(float* output_stack, int shape[3])
{
    int nx = shape[0];
	long long int nslice = shape[0]*shape[1];
	long long int nstack = shape[2]*nslice;

	uint64_t* scan_positions = (uint64_t*) calloc(nstack, sizeof(*scan_positions));

    //std::vector<uint64_t> scan_positions(nstack, 0);

    //scan positions
    long long int pos0 = 0;
    for (int z = 0; z < shape[2]; z++)
    {
        for (int x = 0; x < shape[0]; x++)
        {
            for (long long int ypos = z*nslice+x; ypos < (z+1)*nslice; ypos += nx)
            {
                scan_positions[pos0] = ypos;
                pos0++;
            }
        }
    }

    uint64_t n_lines = nstack/shape[1]; //number of columns in second dimension to be processed

    #pragma omp parallel
    {
        #pragma omp for
        for (long long int n = 0; n < n_lines; n++)
        {
            long long int pos = n*shape[1];
            int r_max, r_start, r_end;
            uint64_t r_temp;
            std::vector<float> workArray(shape[1],0);

            for (int y = 0; y < shape[1]; y++)
            {
                workArray[y] = output_stack[scan_positions[pos]];
                pos++;
            }
            pos -= shape[1];

            for (int y = 0; y < shape[1]; y++)
            {
                uint64_t dsecond = workArray[y];
                if (dsecond != 0)
                {
                    r_max = sqrt(dsecond)+1;
                    r_start = std::min(r_max, y);
                    r_end = std::min(r_max, shape[1]-y);

                    for(int r = -r_start; r < r_end; r++)
                    {
                        r_temp = workArray[y+r]+(r*r);
                        if (r_temp < dsecond)
                            dsecond = r_temp;
                    }
                }

                output_stack[scan_positions[pos]] = dsecond;
                pos++;
            }
        }
    }

    free(scan_positions);
    return;
}
void SaitosAlgorithm::ThirdTransformation(float* output_stack, int shape[3])
{
    double beta2 = beta*beta;
    long long int n_slice = shape[0]*shape[1];
    long long int nstack = shape[2]*n_slice;

    #pragma omp parallel
    {
        #pragma omp for
        for (uint64_t idx = 0; idx < n_slice; idx++)
        {
            int r_max, r_start, r_end;
            uint64_t r_temp;
            std::vector<float> workArray(shape[2],0);

            for (uint64_t z = 0; z < shape[2]; z++)
            {
                workArray[z] = output_stack[idx+z*n_slice];
            }
            for (int z = 0; z < shape[2]; z++)
            {
                uint64_t dthird = workArray[z];
                if (dthird != 0)
                {
                    r_max = sqrt(dthird)+1;
                    r_start = std::min(r_max, z);
                    r_end = std::min(r_max, shape[2]-z);

                    for(int r = -r_start; r < r_end; r++)
                    {
                        r_temp = workArray[z+r]+(beta2*r*r);
                        if (r_temp < dthird)
                            dthird = r_temp;
                    }
                }
                output_stack[idx+z*n_slice] = dthird;
            }
        }
    }
    return;
}
}

