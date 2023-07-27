# LocalThicknessTransform
3D local thickness in C++ for labeled 8bit tif-image sequences

The program should build on any Linux distribution by executing the script *build_localthicknesstransform.sh*.
<br>
The newly create binary needs to be executed with two mandatory arguments (*-i* and *-o* followed by a **full absolute path** to the data / desired output directory. There is now safety check for faulty inputs.)
<br>
The output is a 32bit tif-image sequence providing the 3D local thickness as diameter of a maximal inscribed sphere.

| argument | type | explanation |
|----------|------|-------------|
| -i       | string | full path to input directory which needs to contain an 8bit tif-image sequence. |
| -o       | string | full path to desired output directory. |
| -color   | uint8  | grayscale value of the label to be evaluated. Default is 255. |
| -r_min   | int    | Default is 1. Allows discarding incscribed spheres with small radii. |
| -vxl     | float  | Default 1.0. Allows setting the voxel size for scaling the output to physical length scales. |

@Florian: It might be beneficial to add a CPU limit.
