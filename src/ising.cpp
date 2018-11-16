/*
   2d Ising model - importance sampling
Date: 15 Sep 2016
Author: Sachin Krishnan T V (sachin@physics.iitm.ac.in)
*/

#include <iostream>
#include <fstream>
#include <random>
#include <cstring>
#include <cstdio>
#include <png.h>
#include "../pcg_random/pcg_random.hpp"

#define N 80  // Size of the system
#define MCS 10000  // Number of independent configurations
//#define INT 1.5 // Magnitude of interaction potential
#define EXT 0.0 // Magnitude of external field

//#define KBT 1.0

double JOverkBT;

// This captures the configuration. false -> down, true -> up.
bool lattice[N][N];

using namespace std;

void outputToPng(int t)
{
    int x,y;
    int width = N, height = N;
    png_byte bit_depth = 8;

    png_structp png_ptr;
    png_infop info_ptr;
    png_bytep *row_pointers;

    char filename[256];

    sprintf(filename, "dump_%03d.png", t);
    FILE *dumpf = fopen(filename, "wb");
    if(!dumpf)
        cout<<"[outputToPng] File could not be opened for writing."<<endl;

    png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if(!png_ptr)
        cout<<"[outputToPng] png_create_write_struct failed."<<endl;

    info_ptr = png_create_info_struct(png_ptr);
    if(!info_ptr)
        cout<<"[outputToPng] png_create_info_struct failed."<<endl;

    if(setjmp(png_jmpbuf(png_ptr)))
        cout<<"[outputToPng] Error during init_io."<<endl;

    png_init_io(png_ptr, dumpf);

    if(setjmp(png_jmpbuf(png_ptr)))
        cout<<"[outputToPng] Error during writing header."<<endl;

    png_set_IHDR(png_ptr, info_ptr, width, height,
            bit_depth, PNG_COLOR_TYPE_RGBA, PNG_INTERLACE_NONE,
            PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);

    png_write_info(png_ptr, info_ptr);

    if(setjmp(png_jmpbuf(png_ptr)))
        cout<<"[outputToPng] Error while writing bytes."<<endl;

    row_pointers = (png_bytep*) malloc(sizeof(png_bytep) * height);
    for(y = 0;y < height;y++)
    {
        row_pointers[y] = (png_byte*) malloc(png_get_rowbytes(png_ptr, info_ptr));
        for(x = 0; x < width;x++)
        {
            if(lattice[x][y])
            {
                row_pointers[y][x*4] = 0;
                row_pointers[y][x*4+1] = 0;
                row_pointers[y][x*4+2] = 0;
                row_pointers[y][x*4+3] = 255;
            }
            else
            {
                row_pointers[y][x*4] = 255;
                row_pointers[y][x*4+1] = 255;
                row_pointers[y][x*4+2] = 255;
                row_pointers[y][x*4+3] = 255;
            }
        }
    }

    png_write_image(png_ptr, row_pointers);

    if(setjmp(png_jmpbuf(png_ptr)))
        cout<<"[outputToPng] Error during end of write."<<endl;

    png_write_end(png_ptr, NULL);

    for(y = 0;y < height;y++)
        free(row_pointers[y]);
    free(row_pointers);

    fclose(dumpf);
}

int main(int argc, char *argv[])
{
    int i, j, k, s;
    int il, ir;
    int jl, jr;
    int rand_site;

    ofstream logf;

    pcg_extras::seed_seq_from<random_device> seed_source;
    pcg32 rng(seed_source);
    uniform_int_distribution<int> udist(0, 1);
    uniform_int_distribution<int> site(0, N*N-1);
    uniform_real_distribution<double> metro(0.0,1.0);

    double initial_energy, final_energy, dE;
    double energy;
    double mag_per_step;

    if(argc>1)
    {
        JOverkBT = atof(argv[1]);
    }
    else
    {
        cout<<"Usage : "<<argv[0]<<" <kBT>"<<endl;
        return 0;
    }
    logf.open("log.dat", ios::out);

    // Initial state
    // Taking all down/up spin as initial state is not a good idea, because that is the ground state of the system
    for(i = 0; i < N;i++)
        for(j = 0;j < N;j++)
            lattice[i][j] = static_cast<bool>(udist(rng));

    // Main MC Loop
    for(k = 0; k < MCS;k++)
    {
        for(s = 0;s < N*N;s++)
        {
            rand_site = site(rng);
            i = rand_site / N;
            j = rand_site % N;

            il = (i==0)?N-1:i-1;
            ir = (i==N-1)?0:i+1;

            jl = (j==0)?N-1:j-1;
            jr = (j==N-1)?0:j+1;

            initial_energy = ((lattice[i][j] ^ lattice[il][j])?JOverkBT:-JOverkBT);
            initial_energy += ((lattice[i][j] ^ lattice[ir][j])?JOverkBT:-JOverkBT);
            initial_energy += ((lattice[i][j] ^ lattice[i][jl])?JOverkBT:-JOverkBT);
            initial_energy += ((lattice[i][j] ^ lattice[i][jr])?JOverkBT:-JOverkBT);
            initial_energy += -EXT * ((lattice[i][j])?1:-1);

            lattice[i][j] = !lattice[i][j];

            final_energy = ((lattice[i][j] ^ lattice[il][j])?JOverkBT:-JOverkBT);
            final_energy += ((lattice[i][j] ^ lattice[ir][j])?JOverkBT:-JOverkBT);
            final_energy += ((lattice[i][j] ^ lattice[i][jl])?JOverkBT:-JOverkBT);
            final_energy += ((lattice[i][j] ^ lattice[i][jr])?JOverkBT:-JOverkBT);
            final_energy += -EXT * ((lattice[i][j])?1:-1);

            dE = final_energy - initial_energy;

            if(metro(rng) > exp(-dE))
                lattice[i][j] = !lattice[i][j];

        }

        energy = 0.0;
        mag_per_step = 0.0;
        for(i = 0; i < N;i++)
        {
            il = (i==0)?N-1:i-1;
            ir = (i==N-1)?0:i+1;
            for(j = 0; j < N;j++)
            {
                jl = (j==0)?N-1:j-1;
                jr = (j==N-1)?0:j+1;

                energy += ((lattice[i][j] ^ lattice[il][j])?JOverkBT:-JOverkBT);
                energy += ((lattice[i][j] ^ lattice[ir][j])?JOverkBT:-JOverkBT);
                energy += ((lattice[i][j] ^ lattice[i][jl])?JOverkBT:-JOverkBT);
                energy += ((lattice[i][j] ^ lattice[i][jr])?JOverkBT:-JOverkBT);

                energy += -EXT * ((lattice[i][j])?1:-1);
                mag_per_step += ((lattice[i][j])?1:-1);
            }
        }
        mag_per_step /= (static_cast<double>(N*N));
        logf<<k<<"\t"<<energy<<"\t"<<mag_per_step<<endl;

        // Output configuration
        if(k%10==0)
            outputToPng(k/10);
    }

    logf.close();
    return 0;
}
