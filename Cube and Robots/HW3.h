#pragma once
#ifndef PHYSICSSIM
#define PHYSICSSIM

#include <iostream>
#include <cmath>
#include <vector>
#include <GL/glut.h>
#include <stdarg.h>
#include <numeric>
#include <cstdlib>
#include <ctime>
#include <cstdio>
#include <fstream>

#include <algorithm>




#define LEN 8192  //  Maximum length of text string
#define M_PI 3.1415926
#define Cos(th) cos(M_PI/180*(th))
#define Sin(th) sin(M_PI/180*(th))
#define CUBENUM 30
#define MASSNUM 27
#define SPRINGNUM 351
#define POPNUM 15
#define CHROMENUM 351

std::vector<double> mutation(std::vector<double> child);
void cubeDistance(int popsize);
void best_gene_init();
struct Mass {
    double m;       //mass
    double p[3];    //position
    double v[3];    //velocity
    double a[3];    //acceleration
    double f[3];    //force
};

struct Spring {
    double k;       //constriant
    double L0;      //rest length
    double oriL;
    int m1;
    int m2;
};

struct Cube {
    std::vector<Mass> cube_mass;
    std::vector<Spring> cube_spring;
};

struct Gene {
    double b;
    double c;
};




#endif

//
//#include<GL/glut.h>
////#include <GL/glew.h>
////#include <GL/glu.h>
//#include<cmath>
//#include <vector>
//#include <stdarg.h>
//#include <numeric>
//#include <cstdlib>
//#include <ctime>
//#include <cstdio>
//#include <fstream>
//#include <iostream>
//#include <iterator>
//#include <list>
//
//using namespace std;
//
//int th = 0;            //  Azimuth of view angle
//int ph = 0;           //  Elevation of view angle
//int axes = 1;         //  Display axes
//int light = 1;
//double asp = 1;     // aspect ratio
//int fov = 45;         //  Field of view (for perspective)
//double dim = 0.7;  // size of the workd
//int time_s = 20;
//int emission = 100;  // Emission intensity (%)
//int ambient = 100;  // Ambient intensity (%)
//int diffuse = 100;  // Diffuse intensity (%)
//int specular = 100;  // Specular intensity (%)
//int shininess = 128;  // Shininess (power of two)
//double shiny = 1;    // Shininess (value)
//double white[] = { 1,1,1,1 };
//double black[] = { 0,0,0,1 };
//
//unsigned int grassTexture;
//unsigned int slimeTexture;
//unsigned int skyBoxTexture[10]; // Texture for Sky Box
//
//// Physics Simluator Variables
//
//double mass = 1;
//double mass_2 = 2;
//double length = 0.1;
//double gravity = 9.81;
//double T = 0;
//double springPotentialEnergy = 0;
//double gravityPotentialEnergy = 0;
//double totalPotentialEnergy = 0;
//double kineticEnergy = 0;
//double totalEnergy = 0;
//double groundEnergy = 0;
//double springenergy = 0;
//double gravityenergy = 0;
//double allpotential = 0;
//double kineticenergy = 0;
//double allEnergy = 0;
//double groundenergy = 0;
//double timeStep = 0.002;
//double Nground = 70000;
//double k = 100000;
//double k1 = 20000;
//double k2 = 5000;
//double k3 = 200000;
////double R1 = 0.0 / 255;
////double G1 = 100.0 / 255;
////double B1 = 0.0 / 255;
//
//double R1 = 0.0 / 255;
//double G1 = 0.0 / 255;
//double B1 = 0.0 / 255;
//double dampening = 1;
//int evaluation = 0;
//
//int size = 4;
//int generation = 1000;
//double frictionCoefficient = 0.6;//0.8;
//int oneforce = 0;
//
//double R = 46.0/255;
//double G = 139.0/255;
//double B = 87.0/255;
////int num_sp = 72;
//#define LEN 8192  //  Maximum length of text string
//# define M_PI          3.141592653589793238462643383279502884
//
//#define Cos(th) cos(M_PI/180*(th))
//#define Sin(th) sin(M_PI/180*(th))
//
//
//
//
//
//void Fatal(const char* format , ...)
//{
//   va_list args;
//   va_start(args,format);
//   vfprintf(stderr,format,args);
//   va_end(args);
//   exit(1);
//}
//void ErrCheck(const char* where)
//{
//   int err = glGetError();
//   if (err) fprintf(stderr,"ERROR: %s [%s]\n",gluErrorString(err),where);
//}
//static void Reverse(void* x,const int n)
//{
//   int k;
//   char* ch = (char*)x;
//   for (k=0;k<n/2;k++)
//   {
//      char tmp = ch[k];
//      ch[k] = ch[n-1-k];
//      ch[n-1-k] = tmp;
//   }
//}
//
///*
// *  Load texture from BMP file
// */
//unsigned int LoadTexBMP(const char* file)
//{
//   unsigned int   texture;    // Texture name
//   FILE*          f;          // File pointer
//   unsigned short magic;      // Image magic
//   unsigned int   dx,dy,size; // Image dimensions
//   unsigned short nbp,bpp;    // Planes and bits per pixel
//   unsigned char* image;      // Image data
//   unsigned int   off;        // Image offset
//   unsigned int   k;          // Counter
//   int            max;        // Maximum texture dimensions
//
//   //  Open file
//   f = fopen(file,"rb");
//   if (!f) Fatal("Cannot open file %s\n",file);
//   //  Check image magic
//   if (fread(&magic,2,1,f)!=1) Fatal("Cannot read magic from %s\n",file);
//   if (magic!=0x4D42 && magic!=0x424D) Fatal("Image magic not BMP in %s\n",file);
//   //  Read header
//   if (fseek(f,8,SEEK_CUR) || fread(&off,4,1,f)!=1 ||
//       fseek(f,4,SEEK_CUR) || fread(&dx,4,1,f)!=1 || fread(&dy,4,1,f)!=1 ||
//       fread(&nbp,2,1,f)!=1 || fread(&bpp,2,1,f)!=1 || fread(&k,4,1,f)!=1)
//     Fatal("Cannot read header from %s\n",file);
//   //  Reverse bytes on big endian hardware (detected by backwards magic)
//   if (magic==0x424D)
//   {
//      Reverse(&off,4);
//      Reverse(&dx,4);
//      Reverse(&dy,4);
//      Reverse(&nbp,2);
//      Reverse(&bpp,2);
//      Reverse(&k,4);
//   }
//   //  Check image parameters
//   glGetIntegerv(GL_MAX_TEXTURE_SIZE,&max);
//   if (dx<1 || dx>max) Fatal("%s image width %d out of range 1-%d\n",file,dx,max);
//   if (dy<1 || dy>max) Fatal("%s image height %d out of range 1-%d\n",file,dy,max);
//   if (nbp!=1)  Fatal("%s bit planes is not 1: %d\n",file,nbp);
//   if (bpp!=24) Fatal("%s bits per pixel is not 24: %d\n",file,bpp);
//   if (k!=0)    Fatal("%s compressed files not supported\n",file);
//#ifndef GL_VERSION_2_0
//   //  OpenGL 2.0 lifts the restriction that texture size must be a power of two
//   for (k=1;k<dx;k*=2);
//   if (k!=dx) Fatal("%s image width not a power of two: %d\n",file,dx);
//   for (k=1;k<dy;k*=2);
//   if (k!=dy) Fatal("%s image height not a power of two: %d\n",file,dy);
//#endif
//
//   //  Allocate image memory
//   size = 3*dx*dy;
//   image = (unsigned char*) malloc(size);
//   if (!image) Fatal("Cannot allocate %d bytes of memory for image %s\n",size,file);
//   //  Seek to and read image
//   if (fseek(f,off,SEEK_SET) || fread(image,size,1,f)!=1) Fatal("Error reading data from image %s\n",file);
//   fclose(f);
//   //  Reverse colors (BGR -> RGB)
//   for (k=0;k<size;k+=3)
//   {
//      unsigned char temp = image[k];
//      image[k]   = image[k+2];
//      image[k+2] = temp;
//   }
//
//   //  Sanity check
//   ErrCheck("LoadTexBMP");
//   //  Generate 2D texture
//   glGenTextures(1,&texture);
//   glBindTexture(GL_TEXTURE_2D,texture);
//   //  Copy image
//   glTexImage2D(GL_TEXTURE_2D,0,3,dx,dy,0,GL_RGB,GL_UNSIGNED_BYTE,image);
//   if (glGetError()) Fatal("Error in glTexImage2D %s %dx%d\n",file,dx,dy);
//   //  Scale linearly when image size doesn't match
//   glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
//   glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
//
//   //  Free image memory
//   free(image);
//   //  Return texture name
//   return texture;
//}
//
