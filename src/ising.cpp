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

int SIZE;  // Size of the system
int MCS; // Number of independent configurations
double INT; //Interaction strength
double EXT; // Magnitude of external field

// This captures the configuration. false -> down, true -> up.
bool *lattice;

using namespace std;

void outputToPng(int t)
{
    int x,y;
    int width = SIZE, height = SIZE;
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
            if(lattice[x*SIZE+y])
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

    png_destroy_write_struct(&png_ptr,&info_ptr);

    fclose(dumpf);
}

int main(int argc, char *argv[])
{
    if(argc>=5)
    {
        SIZE = atoi(argv[1]);
        MCS = atoi(argv[2]);
        INT = atof(argv[3]);
        EXT = atof(argv[4]);
    }
    else
    {
        cout<<"Usage : "<<argv[0]<<" <Size> <MCSweeps> <InteractionStrength> <ExternalFieldStrength>"<<endl;
        return 0;
    }
    int i, j, k, s;
    int il, ir;
    int jl, jr;
    int rand_site;

    ofstream logf;

    pcg_extras::seed_seq_from<random_device> seed_source;
    pcg32 rng(seed_source);
    uniform_int_distribution<int> udist(0, 1);
    uniform_int_distribution<int> site(0, SIZE*SIZE-1);
    uniform_real_distribution<double> metro(0.0,1.0);

    double initial_energy, final_energy, dE;
    double energy;
    double mag_per_step;

    lattice = (bool *)malloc(SIZE*SIZE*sizeof(bool));
    if(lattice == nullptr)
    {
        cout<<"Could not allocate memory for lattice. Qutting.\n";
        return 1;
    }
    logf.open("log.dat", ios::out);
    logf<<"#timestep\tenergy\tmagnetization_per_site\n";
    // Initial state
    // Taking all down/up spin as initial state is not a good idea, because that is the ground state of the system
    for(i = 0; i < SIZE*SIZE;i++)
        lattice[i] = static_cast<bool>(udist(rng));

    bool lsite, lleft, lright, lup, ldown;

    // Main MC Loop
    for(k = 0; k <= MCS;k++)
    {
        cout<<"\rTimestep : "<<k<<" / "<<MCS<<"\t";
        fflush(stdout);
        for(s = 0;s < SIZE*SIZE;s++)
        {
            rand_site = site(rng);
            i = rand_site / SIZE;
            j = rand_site % SIZE;

            il = (i==0)?SIZE-1:i-1;
            ir = (i==SIZE-1)?0:i+1;

            jl = (j==0)?SIZE-1:j-1;
            jr = (j==SIZE-1)?0:j+1;
            lsite = lattice[rand_site];
            lleft = lattice[i*SIZE+jl];
            lright = lattice[i*SIZE+jr];
            lup = lattice[il*SIZE+j];
            ldown = lattice[ir*SIZE+j];

            initial_energy = ((lsite ^ lup)?INT:-INT);
            initial_energy += ((lsite ^ ldown)?INT:-INT);
            initial_energy += ((lsite ^ lleft)?INT:-INT);
            initial_energy += ((lsite ^ lright)?INT:-INT);
            initial_energy += -EXT * ((lsite)?1:-1);

            lsite = !lsite;

            final_energy = ((lsite ^ lup)?INT:-INT);
            final_energy += ((lsite ^ ldown)?INT:-INT);
            final_energy += ((lsite ^ lleft)?INT:-INT);
            final_energy += ((lsite ^ lright)?INT:-INT);
            final_energy += -EXT * ((lsite)?1:-1);

            dE = final_energy - initial_energy;

            if(metro(rng) < exp(-dE)) // MC trial accepted
                lattice[rand_site] = lsite;
        }

        energy = 0.0;
        mag_per_step = 0.0;
        for(s = 0; s < SIZE*SIZE;s++)
        {
            i = s / SIZE;
            j = s % SIZE;

            il = (i==0)?SIZE-1:i-1;
            ir = (i==SIZE-1)?0:i+1;

            jl = (j==0)?SIZE-1:j-1;
            jr = (j==SIZE-1)?0:j+1;
            lsite = lattice[s];
            lleft = lattice[i*SIZE+jl];
            lright = lattice[i*SIZE+jr];
            lup = lattice[il*SIZE+j];
            ldown = lattice[ir*SIZE+j];

            energy += ((lsite ^ lup)?INT:-INT);
            energy += ((lsite ^ ldown)?INT:-INT);
            energy += ((lsite ^ lleft)?INT:-INT);
            energy += ((lsite ^ lright)?INT:-INT);

            energy += -EXT * ((lsite)?1:-1);
            mag_per_step += ((lsite)?1:-1);
        }
        mag_per_step /= (static_cast<double>(SIZE*SIZE));
        logf<<k<<"\t"<<energy<<"\t"<<mag_per_step<<endl;

        // Output configuration
        if(k%10==0)
            outputToPng(k/10);
    }

    logf.close();
    free(lattice);
    cout<<endl;
    return 0;
}
