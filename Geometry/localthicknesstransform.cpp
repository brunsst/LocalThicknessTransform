#include "localthicknesstransform.h"
#include "auxiliary.h"
#include <omp.h>

namespace locthick
{
    sedm_type* LocalThicknessTransform::eliminatebydirectneighbour(sedm_type* &sedm)
    {
        double value;
        uint64_t n_slice = shape[0]*shape[1];
        uint64_t n_row = shape[0];
        long long int nstack = shape[2]*n_slice;

        sedm_type* valid = (sedm_type*) calloc(nstack,sizeof(*valid));

        #pragma omp parallel for
        for (long long int idx = 0; idx < nstack; idx++)
            valid[idx] = 1.;

        int radius = 1;
        double sqrt2 = sqrt((double) radius*(double) radius + (double) radius*(double) radius);
        double sqrt3 = sqrt((double) radius*(double) radius + (double) radius*(double) radius + (double) radius*(double) radius);
        uint64_t cutoff;

        for (int z = 0; z < shape[2]; z ++)
        {
            uint64_t idx = z*n_slice;

            for (int y = 0; y < shape[1]; y ++)
            {
                for (int x = 0; x < shape[0]; x++)
                {
                        if(sedm[idx] == 0)
                        {
                            idx++;
                            continue;
                        }
                        value = sqrt(sedm[idx]);
                        cutoff = (uint64_t) round((value+radius)*(value+radius));

                        if((x < shape[0]-radius) && (cutoff <= sedm[idx+radius]))         {valid[idx] = 0; idx++; continue;}
                        if((x > radius-1) &&        (cutoff <= sedm[idx-radius]))         {valid[idx] = 0; idx++; continue;}
                        if((y < shape[1]-radius) && (cutoff <= sedm[idx+radius*n_row]))   {valid[idx] = 0; idx++; continue;}
                        if((y > radius-1) &&        (cutoff <= sedm[idx-radius*n_row]))   {valid[idx] = 0; idx++; continue;}
                        if((z < shape[2]-radius) && (cutoff <= sedm[idx+radius*n_slice])) {valid[idx] = 0; idx++; continue;}
                        if((z > radius-1) &&        (cutoff <= sedm[idx-radius*n_slice])) {valid[idx] = 0; idx++; continue;}

                        cutoff = (uint64_t) round((value+sqrt2)*(value+sqrt2));

                        if((x < shape[0]-radius) && (y < shape[1]-radius) && (cutoff <= sedm[idx+radius+radius*n_row])) {valid[idx] = 0; idx++; continue;}
                        if((x > radius-1) && (y < shape[1]-radius) &&        (cutoff <= sedm[idx-radius+radius*n_row])) {valid[idx] = 0; idx++; continue;}
                        if((x < shape[0]-radius) && (y > radius-1) &&        (cutoff <= sedm[idx+radius-radius*n_row])) {valid[idx] = 0; idx++; continue;}
                        if((x > radius-1) && (y > radius-1) &&               (cutoff <= sedm[idx-radius-radius*n_row])) {valid[idx] = 0; idx++; continue;}

                        if((x < shape[0]-radius) && (z < shape[2]-radius) && (cutoff <= sedm[idx+radius+radius*n_slice])) {valid[idx] = 0; idx++; continue;}
                        if((x > radius-1) && (z < shape[2]-radius) &&        (cutoff <= sedm[idx-radius+radius*n_slice])) {valid[idx] = 0; idx++; continue;}
                        if((x < shape[0]-radius) && (z > radius-1) &&        (cutoff <= sedm[idx+radius-radius*n_slice])) {valid[idx] = 0; idx++; continue;}
                        if((x > radius-1) && (z > radius-1) &&               (cutoff <= sedm[idx-radius-radius*n_slice])) {valid[idx] = 0; idx++; continue;}

                        if((z < shape[2]-radius) && (y < shape[1]-radius) && (cutoff <= sedm[idx+radius*n_slice+radius*n_row])) {valid[idx] = 0; idx++; continue;}
                        if((z > radius-1) && (y < shape[1]-radius) &&        (cutoff <= sedm[idx-radius*n_slice+radius*n_row])) {valid[idx] = 0; idx++; continue;}
                        if((z < shape[2]-radius) && (y > radius-1) &&        (cutoff <= sedm[idx+radius*n_slice-radius*n_row])) {valid[idx] = 0; idx++; continue;}
                        if((z > radius-1) && (y > radius-1) &&               (cutoff <= sedm[idx-radius*n_slice-radius*n_row])) {valid[idx] = 0; idx++; continue;}

                        cutoff = (uint64_t) round((value+sqrt3)*(value+sqrt3));

                        if((x < shape[0]-radius) && (y < shape[1]-radius) && (z < shape[2]-radius) && (cutoff <= sedm[idx+radius+radius*n_row+radius*n_slice])) {valid[idx] = 0; idx++; continue;}
                        if((x > radius-1) && (y < shape[1]-radius) && (z < shape[2]-radius) &&        (cutoff <= sedm[idx-radius+radius*n_row+radius*n_slice])) {valid[idx] = 0; idx++; continue;}
                        if((x < shape[0]-radius) && (y > radius-1) && (z < shape[2]-radius) &&        (cutoff <= sedm[idx+radius-radius*n_row+radius*n_slice])) {valid[idx] = 0; idx++; continue;}
                        if((x < shape[0]-radius) && (y < shape[1]-radius) && (z > radius-1) &&        (cutoff <= sedm[idx+radius+radius*n_row-radius*n_slice])) {valid[idx] = 0; idx++; continue;}
                        if((x > radius-1) && (y > radius-1) && (z < shape[2]-radius) &&               (cutoff <= sedm[idx-radius-radius*n_row+radius*n_slice])) {valid[idx] = 0; idx++; continue;}
                        if((x > radius-1) && (y < shape[1]-radius) && (z > radius-1) &&               (cutoff <= sedm[idx-radius+radius*n_row-radius*n_slice])) {valid[idx] = 0; idx++; continue;}
                        if((x < shape[0]-radius) && (y > radius-1) && (z > radius-1) &&               (cutoff <= sedm[idx+radius-radius*n_row-radius*n_slice])) {valid[idx] = 0; idx++; continue;}
                        if((x > radius-1) && (y > radius-1) && (z > radius-1) &&                      (cutoff <= sedm[idx-radius-radius*n_row-radius*n_slice])) {valid[idx] = 0; idx++; continue;}

                    idx++;
                }
            }
        }

        for (uint64_t i = 0; i < nstack; i++)
        {
            valid[i] = valid[i]*sedm[i];
        }
        return valid;
    }
    void LocalThicknessTransform::drawlocthickmap_maximaonly(sedm_type* &valid, sedm_type* &locthickmap)
    {
        uint64_t n_points_evaluated = 0;
        uint64_t removed = 0;

        long long int nslice = shape[0]*shape[1];
        long long int nstack = shape[2]*nslice;

        #pragma omp parallel
        {
            #pragma omp for
            for (uint64_t idx = 0; idx < nstack; idx++)
            {
                int x,y,z;
                bool change;

                if (valid[idx] != 0)
                {
                    aux::idx2xyz(idx,shape,x,y,z);

                    if ((x > 0) && valid[idx-1] > valid[idx]){continue;}
                    if ((y > 0) && valid[idx-n_row] > valid[idx]){continue;}
                    if ((z > 0) && valid[idx-n_slice] > valid[idx]){continue;}
                    if ((x < shape[0]-1) && valid[idx+1] > valid[idx]){continue;}
                    if ((y < shape[1]-1) && valid[idx+n_row] > valid[idx]){continue;}
                    if ((z < shape[2]-1) && valid[idx+n_slice] > valid[idx]){continue;}

                    change = evaluateposition(valid,locthickmap,x,y,z,idx,valid[idx]);
                    if (change == false)
                        valid[idx] = 0;
                }
            }
        }
        return;
    }
    sedm_type* LocalThicknessTransform::drawlocthickmap_maximaonly_removeandbackup(sedm_type* &valid, sedm_type* &locthickmap)
    {
        uint64_t n_points_evaluated = 0;
        uint64_t removed = 0;

        long long int nslice = shape[0]*shape[1];
        long long int nstack = shape[2]*nslice;
        sedm_type* backup = (sedm_type*) calloc(nstack,sizeof(*backup));

        #pragma omp parallel
        {
            #pragma omp for
            for (uint64_t idx = 0; idx < nstack; idx++)
            {
                int x,y,z;
                bool change;

                if (valid[idx] != 0)
                {
                    aux::idx2xyz(idx,shape,x,y,z);

                    if ((x > 0) && valid[idx-1] > valid[idx]){continue;}
                    if ((y > 0) && valid[idx-n_row] > valid[idx]){continue;}
                    if ((z > 0) && valid[idx-n_slice] > valid[idx]){continue;}
                    if ((x < shape[0]-1) && valid[idx+1] > valid[idx]){continue;}
                    if ((y < shape[1]-1) && valid[idx+n_row] > valid[idx]){continue;}
                    if ((z < shape[2]-1) && valid[idx+n_slice] > valid[idx]){continue;}

                    change = evaluateposition(valid,locthickmap,x,y,z,idx,valid[idx]);
                    if (change == false)
                        valid[idx] = 0;
                    else
                    {
                        backup[idx] = valid[idx];
                        valid[idx] = 0;
                    }
                }
            }
        }
        return backup;
    }
    void LocalThicknessTransform::restorefrombackup(sedm_type* &valid, sedm_type* &backup)
    {
        long long int nslice = shape[0]*shape[1];
        long long int nstack = shape[2]*nslice;

        for (int i = 0; i < nstack; i++)
        {
            if(backup[i] > 0)
                valid[i] = backup[i];
        }
    }
    sedm_type* LocalThicknessTransform::drawlocthickmap_sorted(sedm_type* &valid, sedm_type* &locthickmap, uint64_t n_validcenters)
    {
        std::vector<IndexSorting> sorted_centers;
        sorted_centers.reserve(n_validcenters);

        long long int nslice = shape[0]*shape[1];
        long long int nstack = shape[2]*nslice;

        for(uint64_t idx = 0; idx < nstack; idx++)
        {
            if(valid[idx] > 0)
                sorted_centers.push_back({idx,(uint64_t) valid[idx]});
        }
        std::vector<uint64_t> sorted_idx;
        sorted_idx.reserve(n_validcenters);
        std::sort(sorted_centers.begin(), sorted_centers.end());
        for(uint64_t idx = 0; idx < sorted_centers.size(); idx++)
        {
            sorted_idx.push_back(sorted_centers[idx].idx);
        }

        #pragma omp parallel
        {
            #pragma omp for
            for (uint64_t i = 0; i < sorted_idx.size(); i++)
            {
                uint64_t idx = sorted_idx[i];
                if(valid[idx] > 0)
                {
                    int x,y,z;
                    bool change;

                    if (valid[idx] != 0)
                    {
                        aux::idx2xyz(idx,shape,x,y,z);
                        change = evaluateposition(valid,locthickmap,x,y,z,idx,valid[idx]);
                        if (change == false)
                            valid[idx] = 0;
                    }
                }
            }
        }

        return locthickmap;
    }

    sedm_type* LocalThicknessTransform::drawlocthickmap(sedm_type* &valid, sedm_type* &locthickmap)
    {
        uint64_t n_points_evaluated = 0;
        uint64_t removed = 0;

        long long int nslice = shape[0]*shape[1];
        long long int nstack = shape[2]*nslice;

        #pragma omp parallel
        {
            #pragma omp for
            for (uint64_t idx = 0; idx < nstack; idx++)
            {
                int x,y,z;
                bool change;

                if (valid[idx] != 0)
                {
                    aux::idx2xyz(idx,shape,x,y,z);
                    change = evaluateposition(valid,locthickmap,x,y,z,idx,valid[idx]);
                    if (change == false)
                        valid[idx] = 0;
                }
            }
        }

        return locthickmap;
    }
    uint64_t LocalThicknessTransform::positionstoevaluate(sedm_type* &valid)
    {
        long long int nslice = shape[0]*shape[1];
        long long int nstack = shape[2]*nslice;

        uint64_t counter = 0;
        for(int i = 0; i < nstack; i++)
        {
            if (valid[i] != 0)
                counter++;
        }
        return counter;
    }
    bool LocalThicknessTransform::evaluateposition(sedm_type* &valid, sedm_type* &locthickmap, const int &x, const int &y, const int &z, const uint64_t &origin,
                const uint64_t &squared_radius)
    {
        bool change = false;

        double parent_radius = sqrt(squared_radius);
        double child_radius;
        uint64_t squared_radius_reduced = (parent_radius-.5)*(parent_radius-.5); //reduced to consider surface half-way
        int radius = floor(parent_radius);

        uint64_t idx = origin-radius-shape[0]*radius-shape[0]*shape[1]*radius;
        uint64_t d1,d2,squared_distance,distance,radsum;

        for(int r = -radius; r <= radius; r++)
        {
            if (z+r < 0) {idx+= n_slice; continue;}
            if (z+r >= shape[2]) {idx+= n_slice; continue;}
            d1 = r*r;

            for (int q = -radius; q <= radius; q++)
            {
                if (y+q < 0) {idx+= n_row; continue;}
                if (y+q >= shape[1]) {idx+= n_row; continue;}
                d2 = d1+q*q;

                for(int p = -radius; p <= radius; p++)
                {
                    if (x+p < 0) {idx++; continue;}
                    if (x+p >= shape[0]) {idx++; continue;}

                    squared_distance = d2+p*p;

                    if (valid[idx] != 0)
                    {
                        if (r != 0 || q != 0 || p != 0)
                        {
                            //maybe some voxels can be removed
                            distance = sqrt(squared_distance);
                            child_radius = sqrt(valid[idx]);
                            if (parent_radius >= child_radius+distance)
                            {
                                valid[idx] = 0;
                            }
                            if (child_radius >= parent_radius+distance)
                            {
                                valid[origin] = 0;
                                return false;
                            }
                        }

                    }
                    if (squared_distance <= squared_radius_reduced)
                    {
                        if (squared_radius >= locthickmap[idx])
                        {
                            locthickmap[idx] = squared_radius;
                            change = true;
                        }
                    }
                    idx++;
                }
                idx += n_row-(2*radius+1);
            }
            idx += n_slice-(2*radius+1)*n_row;
        }
        return change;
    }
}

