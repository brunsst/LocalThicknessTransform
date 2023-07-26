#include <iostream>
#include <vector>
#include <algorithm>
#include <omp.h>

namespace aux
{
    /*String-Manipulation
    *********************************************************/
    std::string zfill_int2string(int inint, const unsigned int &zfill)
    {
        std::string outstring = std::to_string(inint);
        while(outstring.length() < zfill)
            outstring = "0" + outstring;
        return outstring;
    }

    /*Numpy-like
    *********************************************************/
    std::vector<float> linspace(float startval, float endval, uint64_t bins)
    {
        std::vector<float> linspaced(bins);
        float delta = (endval-startval)/(bins-1);
        for(uint64_t i = 0; i < (bins-1); i++)
        {
            linspaced[i] = startval + delta * i;
        }
        linspaced[bins-1] = endval;
        return linspaced;
    }
    std::vector<double> linspace(double startval, double endval, uint64_t bins)
	{
		std::vector<double> linspaced(bins);
		double delta = (endval-startval)/(bins-1);
		for(uint64_t i = 0; i < (bins-1); i++)
		{
			linspaced[i] = startval + delta * i;
		}
		linspaced[bins-1] = endval;
		return linspaced;
	}
    /*
     ************************************************************/

    void normalize_greyscales(float* image, int shape[3])
    {
    	long long int nslice = shape[0]*shape[1];
    	long long int nstack = shape[2]*nslice;

    	float maxval = -1e21;
    	float minval = 1e21;

		#pragma omp parallel for reduction(max: maxval), reduction(min: minval)
    	for (long long int idx = 0; idx < nstack; idx++)
    	{
    		float val = image[idx];
    		if (val > maxval) maxval = val;
    		if (val < minval) minval = val;
    	}
		#pragma omp parallel for
    	for (long long int idx = 0; idx < nstack; idx++)
    		image[idx] = (image[idx]-minval)/(maxval-minval);

    	return;
    }
    void idx2xyz(const int64_t &pos, const int shape[3], int &outx, int &outy, int &outz)
    {
    	outz = pos/(shape[0]*shape[1]);
    	int remainder = pos-(outz*shape[0]*shape[1]);
    	outy = remainder/shape[0];
    	outx = remainder-(outy*shape[0]);
    	return;
    }
    float* upscale_bonescrewlabels(float* labels_bone, float* labels_screw, int shape[3], bool bone_preference)
    {
    	int nx = shape[0]; int ny = shape[1]; int nz = shape[2];
    	long long int nslice = nx*ny;
    	long long int nstack = nslice*nz;

    	int nx2 = 2*shape[0]; int ny2 = 2*shape[1]; int nz2 = 2*shape[2];
    	long long int nslice2 = nx2*ny2;
    	long long int nstack2 = nslice2*nz2;

    	float* output = (float*) calloc(nstack2, sizeof(*output));

		#pragma omp parallel for
    	for (long long int idx = 0; idx < nstack; idx++)
    	{
    		int z = idx/nslice;
    		int y = (idx-z*nslice)/nx;
    		int x = idx-z*nslice-y*nx;

    		float label = 0;
    		if (labels_screw[idx] != 0) label = 2;
    		else if (labels_bone[idx] != 0) label = 1;
    		if (bone_preference && labels_bone[idx] != 0) label = 1;

    		output[(2*z)*nslice2 + (2*y)*nx2 + (2*x)] = label;
    		output[(2*z)*nslice2 + (2*y)*nx2 + (2*x+1)] = label;
    		output[(2*z)*nslice2 + (2*y+1)*nx2 + (2*x)] = label;
    		output[(2*z+1)*nslice2 + (2*y)*nx2 + (2*x)] = label;

    		output[(2*z)*nslice2 + (2*y+1)*nx2 + (2*x+1)] = label;
    		output[(2*z+1)*nslice2 + (2*y)*nx2 + (2*x+1)] = label;
    		output[(2*z+1)*nslice2 + (2*y+1)*nx2 + (2*x)] = label;
    		output[(2*z+1)*nslice2 + (2*y+1)*nx2 + (2*x+1)] = label;
    	}

    	return output;
    }

    std::vector<std::vector<float>> get_convhull_vertices(float* image, int shape[3], float foreground_label, uint8_t* output)
    {
        std::vector<std::vector<float>> vertices;

        int nx = shape[0]; int ny = shape[1]; int nz = shape[2];
        long long int nslice = nx*ny;
        long long int nstack = nz*nslice;

        for (long long int idx = 0; idx < nstack; idx++)
        {
        int z = idx/nslice;
        int y = (idx-z*nslice)/nx;
        int x = idx-z*nslice-y*nx;

        float value = image[idx];
        if (value != foreground_label) continue; //not foreground

        output[idx] = 1;

        if (z == 0 || z == nz-1 || y == 0 || y == ny-1 || x == 0 || x == nx-1){
            vertices.push_back({(float) x,(float) y,(float) z});
            continue;
        } //image boundary = on bou8ndary

        //6-connected
        if (image[idx-1] != foreground_label || image[idx+1] != foreground_label) {vertices.push_back({(float) x,(float) y,(float) z}); continue;}
        if (image[idx-nx] != foreground_label || image[idx+nx] != foreground_label) {vertices.push_back({(float) x,(float) y,(float) z}); continue;}
        if (image[idx-nslice] != foreground_label || image[idx+nslice] != foreground_label) {vertices.push_back({(float) x,(float) y,(float) z}); continue;}
        }

        return vertices;
    }
}
