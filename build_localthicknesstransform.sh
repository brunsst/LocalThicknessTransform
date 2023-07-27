#!/bin/bash

mkdir -p "Build"

echo "building LocalThicknessTransform for CPU"
echo "--------------------------------------------------"
echo "compiling auxiliary.cpp"
g++ -w -fexceptions -O3 -std=c++11 -fopenmp  -c $PWD"/Geometry/auxiliary.cpp" -o $PWD"/Build/auxiliary.o"
echo "compiling hdcommunication.cpp"
g++ -w -fexceptions -O3 -std=c++11 -fopenmp  -c $PWD"/Geometry/hdcommunication.cpp" -o $PWD"/Build/hdcommunication.o"
echo "compiling distancemapping.cpp"
g++ -w -fexceptions -O3 -std=c++11 -fopenmp  -c $PWD"/Geometry/distancemapping.cpp" -o $PWD"/Build/distancemapping.o"
echo "compiling localthicknesstransform.cpp"
g++ -w -fexceptions -O3 -std=c++11 -fopenmp  -c $PWD"/Geometry/localthicknesstransform.cpp" -o $PWD"/Build/localthicknesstransform.o"
echo "compiling main_cpu.cpp"
g++ -w -fexceptions -O3 -std=c++11 -fopenmp  -c $PWD"/main.cpp" -o $PWD"/Build/main.o"
echo "linking LocalTHicknessTransform"
g++  -o $PWD/LocalThicknessTransform $PWD/Build/main.o  $PWD/Build/auxiliary.o $PWD/Build/hdcommunication.o $PWD/Build/distancemapping.o  $PWD/Build/localthicknesstransform.o  -ltiff -lgomp
echo "--------------------------------------------------"

################################################################################################################################################################
################################################################################################################################################################
