#ifndef AUXILIARY_H
#define AUXILIARY_H

#include <iostream>
#include <string.h>
#include <cstdint>
#include <vector>

namespace aux
{
    /*String-Manipulation
    *********************************************************/
    std::string zfill_int2string(int inint, const unsigned int &zfill);

    /*Numpy-like
    *********************************************************/
    std::vector<float> linspace(float startval, float endval, uint64_t bins);
    std::vector<double> linspace(double startval, double endval, uint64_t bins);

    void normalize_greyscales(float* image, int shape[3]);
    void idx2xyz(const int64_t &pos, const int shape[3], int &outx, int &outy, int &outz);

    float* upscale_bonescrewlabels(float* labels_bone, float* label_screw, int shape[3], bool bone_preference = false);

    std::vector<std::vector<float>> get_convhull_vertices(float* image, int shape[3], float foreground_label, uint8_t* output);

}

#endif // AUXILIARY_H
