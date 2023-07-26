#ifndef LOCALTHICKNESSTRANSFORM_H
#define LOCALTHICKNESSTRANSFORM_H

#include <vector>
#include <iostream>
#include <cstdint>
#include <omp.h>
#include <algorithm>
#include <cmath>

namespace locthick
{
    typedef float sedm_type;

    struct IndexSorting
    {
        uint64_t idx;
        uint64_t value;

        IndexSorting(uint64_t idx_, uint64_t value_):idx(idx_),value(value_){}

        //sort the results with decreasing correlation:
        bool operator< (const IndexSorting &rhs) const
        {
            return value > rhs.value;
        }
    };

    class LocalThicknessTransform
    {
    public:
        bool only_localthicknessmap = true;
        bool use_sortedapproach = true; //slower but needs less spheres
        bool use_backupapproach = true;
        double minradius = 2;

    private:
        int shape[3]; //shape of the image stack
        uint64_t n_row; //shift in index for the next voxel in dim1
        uint64_t n_slice; //shift in index for the next voxel in dim2

    public:
        void Run(sedm_type* &sedm, int datashape[3], sedm_type* &out_valid_sphere_centers, sedm_type* &out_squared_localthickness_map)
        {
            //Set up some parameters that will be used frequently
            shape[0] = datashape[0];
            shape[1] = datashape[1];
            shape[2] = datashape[2];
            n_row = shape[0];
            n_slice = shape[0]*shape[1];
            long long int nstack = n_slice*shape[2];

            sedm_type* backup;
            uint64_t voxelstoanalyze0 = positionstoevaluate(sedm);

            //Step 1
            std::cout << "    eliminating by direct neighbour...";
            std::cout.flush();
            out_valid_sphere_centers = eliminatebydirectneighbour(sedm);
            uint64_t voxelstoanalyze1 = positionstoevaluate(out_valid_sphere_centers);
            std::cout << "done! (" << round(voxelstoanalyze1/1000.)/1000. << "M/" << round(voxelstoanalyze0/1000.)/1000. << "M possible centers)" << std::endl;

            //Step1.5
            if(minradius > 0)
            {
                std::cout << "    thresholding centers at r=" << minradius << "...";
                std::cout.flush();
                uint64_t minradius_squared = minradius*minradius;
                for(uint64_t i = 0; i < nstack; i++)
                {
                    if(out_valid_sphere_centers[i] < minradius_squared)
                        out_valid_sphere_centers[i] = 0;
                }
                uint64_t voxelstoanalyze_tmp = positionstoevaluate(out_valid_sphere_centers);
                std::cout << "    done! (" << round(voxelstoanalyze_tmp/1000.)/1000. << "M/" << round(voxelstoanalyze1/1000.)/1000. << "M possible centers)" << std::endl;
                voxelstoanalyze1 = voxelstoanalyze_tmp;
            }

            //Step2
            uint64_t voxelstoanalyze3;
            if (use_sortedapproach == true)
            {
                std::cout << "    using sorting approach...         ";
                std::cout.flush();
                drawlocthickmap_sorted(out_valid_sphere_centers, out_squared_localthickness_map, voxelstoanalyze1);
                uint64_t voxelstoanalyze2 = positionstoevaluate(out_valid_sphere_centers);
                std::cout << "done! (" << round(voxelstoanalyze2/1000.)/1000. << "M/" << round(voxelstoanalyze1/1000.)/1000. << "M possible centers)" << std::endl;
                voxelstoanalyze3 = voxelstoanalyze2;
            }
            else
            {
                std::cout << "    eliminating with local maxima...  ";
                std::cout.flush();
                if(use_backupapproach == false)
                    drawlocthickmap_maximaonly(out_valid_sphere_centers, out_squared_localthickness_map);
                else
                    backup = drawlocthickmap_maximaonly_removeandbackup(out_valid_sphere_centers, out_squared_localthickness_map);
                uint64_t voxelstoanalyze2 = positionstoevaluate(out_valid_sphere_centers);
                std::cout << "done! (" << round(voxelstoanalyze2/1000.)/1000. << "M/" << round(voxelstoanalyze1/1000.)/1000. << "M possible centers)" << std::endl;

                std::cout << "    calculating max. inscr. spheres...";
                std::cout.flush();
                drawlocthickmap(out_valid_sphere_centers, out_squared_localthickness_map);
                voxelstoanalyze3 = positionstoevaluate(out_valid_sphere_centers);
                std::cout << "done! (" << round(voxelstoanalyze3/1000.)/1000. << "M/" << round(voxelstoanalyze2/1000.)/1000. << "M possible centers)" << std::endl;
            }

            if(only_localthicknessmap == true)
                return; //early exit if nothing but the map is needed
            if(use_backupapproach == true)
                restorefrombackup(out_valid_sphere_centers, backup);

            //Step4
            //rerunning eliminates redundancies
            std::cout << "    eliminating redundant centers...  ";
            std::cout.flush();
            drawlocthickmap(out_valid_sphere_centers, out_squared_localthickness_map);
            uint64_t voxelstoanalyze4 = positionstoevaluate(out_valid_sphere_centers);
            std::cout << "done! (" << round(voxelstoanalyze4/1000.)/1000. << "M/" << round(voxelstoanalyze3/1000.)/1000. << "M possible centers)" << std::endl;
        }

    private:
        /********** Functions reducing the amount of spheres to be checked **********/
        //First step:
        //Iterate once over the images and eliminate all sedm values that are entirely covered
        //by an immediate neighbour in the 26 neighbourhood.
        //(not performed in parallel)
        sedm_type* eliminatebydirectneighbour(sedm_type* &sedm);

        //2nd step:
        //calculate maximal spheres only from centers that are a local maximum in their 6neighbourhood
        //this should remove a great chunk of possible centers!
        void drawlocthickmap_maximaonly(sedm_type* &valid, sedm_type* &locthickmap);
        sedm_type* drawlocthickmap_maximaonly_removeandbackup(sedm_type* &valid, sedm_type* &locthickmap);
        void restorefrombackup(sedm_type* &valid, sedm_type* &backup);

        //3rd and 4th step:
        //do the previous step with all remaining centers!
        sedm_type* drawlocthickmap(sedm_type* &valid, sedm_type* &locthickmap);
        sedm_type* drawlocthickmap_sorted(sedm_type* &valid, sedm_type* &locthickmap, uint64_t n_validcenters);

        /***************************** Helper functions *****************************/
        uint64_t positionstoevaluate(sedm_type* &valid); //quick check how many spheres are still potentially part of the transform

        //Draws a sphere and checks if it contributes to the map:
        bool evaluateposition(sedm_type* &valid, sedm_type* &locthickmap, const int &x, const int &y, const int &z,
                              const uint64_t &origin, const uint64_t &squared_radius);

    };
}

#endif // LOCALTHICKNESSTRANSFORM_H
