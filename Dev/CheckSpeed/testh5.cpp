#include <iostream>
#include <ctime>
#include <fftw3.h>
#include <H5Cpp.h>
#include "Eigen/Dense"

int main (int argc, char *argv[]) 
{
    int i, j;
    int NX,NY;

    NX = 3;
    NY = 4;

    double data[NX][NY]; // buffer for data to write
    for (j = 0; j < NX; j++)
    {
        for (i = 0; i < NY; i++)
        {
            data[j][i] = double(i + j);
        }
            
    }

    // Eigen::MatrixXd data(NX,NY);
    // M.random()

    H5::H5File F("test.h5",H5F_ACC_TRUNC);
    hsize_t dimsf[2]; // dataset dimensions
    dimsf[0] = NX;
    dimsf[1] = NY;
    H5::IntType datatype(H5::PredType::NATIVE_DOUBLE);
    H5::DataSpace dataspace(2, dimsf);
    H5::DataSet D = F.createDataSet("M",datatype,dataspace);
    D.write(data,H5::PredType::NATIVE_DOUBLE);
}