#include <iostream>
#include <omp.h>
#include <fstream>

#include "Geometry/hdcommunication.h"
#include "Geometry/localthicknesstransform.h"
#include "Geometry/distancemapping.h"

using namespace std;

int main(int argc, char* argv[])
{
    //////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////
    string inpath = "";
    string outpath = "";

    uint8_t active_color = 255;
    int min_radius = 1;
    float vxlsize = 1.0;
    //////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////////////
    if ("extract command line arguments)")
    {
        for (uint16_t i = 1; i < argc; i++)
        {
            if (std::string(argv[i]) == "-i"){
                i++;
                inpath = std::string(argv[i]);
            }
            else if (std::string(argv[i]) == "-o"){
                i++;
                outpath = std::string(argv[i]);
            }
            else if (std::string(argv[i]) == "-color"){
                i++;
                active_color = std::stoi(argv[i]);
            }
            else if (std::string(argv[i]) == "-r_min"){
                i++;
                min_radius = std::stoi(argv[i]);
            }
            else if (std::string(argv[i]) == "-vxl"){
                i++;
                vxlsize = std::stod(argv[i]);
            }
        }
    }
    //////////////////////////////////////////////////////////////////////////////////////

    hdcom::HdCommunication hdcom;
    int shape[3];
    std::vector<std::string> filelist = hdcom.GetFilelist(inpath, shape);
    uint8_t* all_labels = hdcom.Get3DTifSequence_8bitPointer(filelist, shape, true);

    int nx = shape[0]; int ny = shape[1]; int nz = shape[2];
    long long int nslice = nx*ny; long long int nstack = nz*nslice;
    uint8_t* active_labels = (uint8_t*) calloc(nstack, sizeof(*active_labels));

    int first_z = nz-1;
    int last_z = 0;

    //mask out everything but pits
    #pragma omp parallel for reduction(min: first_z), reduction(max: last_z)
    for(long long int idx = 0; idx < nstack; idx++)
    {
        int z = idx/nslice;
        if (all_labels[idx] == active_color)
        {
            if (z < first_z) first_z = z;
            if (z > last_z) last_z = z;
        }

        if(all_labels[idx] == active_color) active_labels[idx] = active_color;
    }
    cout << "z-range: " << first_z << " " << last_z+1 << endl;

    cout << "creating distance map...\r";
    cout.flush();
    distmap::SaitosAlgorithm distmap;
    float* sedm = distmap.SEDM_foreground(active_labels, shape);
    cout << endl;

    cout << "local thickness transform..." << endl;
    locthick::LocalThicknessTransform locthick;
    locthick.minradius = min_radius;

    float* valid_sphere_centers = (float*) calloc(nstack,sizeof(*valid_sphere_centers));
    float* squared_locthick = (float*) calloc(nstack,sizeof(*squared_locthick));
    locthick.Run(sedm, shape, valid_sphere_centers, squared_locthick);

    float* locthick_map = (float*) calloc(nstack, sizeof(*locthick_map));

    #pragma omp parallel for
    for(long long int idx = 0; idx < nstack; idx++)
    {
        float value = sqrt(squared_locthick[idx]);
        if (value == 1)
            locthick_map[idx] = vxlsize;
        else
            locthick_map[idx] = 2*vxlsize*value; //locthick as diameter and in um
    }
    //////////////////////////////////////////////////////////////////////////////////////

    //Create output
    /////////////////////////////////////////////////////////////////////////////////////
    hdcom.makedir(outpath);
    hdcom.SaveTifSequence_32bit(locthick_map, shape, outpath, "locthick", true);

    return 0;
}
