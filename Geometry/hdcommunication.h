#ifndef HDCOMMUNICATION_H
#define HDCOMMUNICATION_H

#include <vector>
#include <string.h>
#include <iostream>
#include <tiffio.h>
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <errno.h>
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <sstream>
#include <omp.h>

namespace hdcom
{
    class HdCommunication
    {
    public:
    	short last_bps = 0;

    	/*************** Reading and saving unknown dimensionality ***************/
    	std::vector<std::string> GetFilelist_And_ImageSequenceDimensions(std::string path, int outshape[3], bool &is_rgb);
    	float* GetTif_unknowndim_32bit(std::string inpath, int outshape[3], bool dspprogress=false);
    	uint8_t* GetTif_unknowndim_8bit(std::string path, int outshape[3], bool dspprogress=false, bool rescale=false);
    	void SaveTif_unknowndim_32bit(float *image, int imgshape[3], std::string path, std::string name, std::string subdir3D = "", int slice_nr = -1);

    	/*************** Read ImageJ 3D-tif and/or read without libtiff ***************/
    	float* Custom3DTifReader(std::string path, int outshape[3], bool verbose = false);

        /*************** Reading a single greyscale image as 1D vector ***************/
        uint8_t* Get2DTifImage_8bitPointer(std::string file, int outshape[2], bool rescale = false);
        float* Get2DTifImage_32bitPointer(std::string file, int outshape[2]);
        void GetRGBTif_8bitChannels(std::string file, int outshape[2], uint8_t* &R, uint8_t* &G, uint8_t* &B);
        void GetRGBTif_32bitChannels(std::string file, int outshape[2], float* &R, float* &G, float* &B);
        void Insert2DTifImage_32bitPointer(std::string file, int outshape[2], float *imgpointer, long long int pos0);
        void Insert2DTifImage_8bitPointer(std::string file, int outshape[2], uint8_t *imgpointer, long long int pos0);

        /**************** Reading a greyscale image sequence (with overloading) ****************/
        float* Get3DTifSequence_32bitPointer(std::vector<std::string> &filelist, int outshape[3], bool dspprogress);
        float* Get3DTifSequence_32bitPointer(std::vector<std::string> &filelist, int outshape[3], int firstslice, int lastslice); //lastslice is inclusive
        uint8_t* Get3DTifSequence_8bitPointer(std::vector<std::string> &filelist, int outshape[3], bool dspprogress);

        /********************** Saving a single greyscale image **********************/
        template <typename T> void Save2DTifImage_8bit(T* image, int imgshape[2], std::string path, std::string name, long long int firstpos = 0, bool rescale = false);
        template <typename T> void Save2DTifImage_32bit(T* image, int imgshape[2], std::string path, std::string name, long long int firstpos = 0);
        void Save2DTifImage_as16bit(float *image, int imgshape[2], std::string path, std::string name, int64_t pos=0);
        void Save2DTifImage_RGB(uint8_t *image, int imgshape[2], std::string path, std::string name);

        /********************** Saving a greyscale sequence **********************/
        template <typename T> void SaveTifSequence_8bit(T* image, int imgshape[3], std::string path, std::string name, bool dspprogress, bool rescale=false);
        template <typename T> void SaveTifSequence_32bit(T *image, int imgshape[3], std::string path, std::string name, bool dspprogress);
        void SaveTifSequence_as16bit(int z0, float *image, int imgshape[3], std::string path, std::string name, bool dspprogress);

        /********************** Miscellaneous **********************/
        void Save3DVector_vtk(float *u, int shape[3], std::string path, std::string name, std::string header="");
        void save2inr(std::string outpath, std::string filename, uint8_t *image, int shape[3]);

        //Get a vector with tif-files in a directory:
        int GetFilelist(std::string const directory, std::vector<std::string> &outfiles);
        std::vector<std::string> GetFilelist(std::string path, int outshape[3]);
        void makedir(const std::string path);

    private:
        bool hasEnding(std::string const &str1, std::string const &str2);

    };
}

#endif // HDCOMMUNICATION_H
