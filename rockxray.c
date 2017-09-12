/*************************************************************
 *
 *  rockxray.c - Create an image of an irregular triangle mesh
 *
 *  Mark J. Stock, mstock@umich.edu
 *
 * rocktools - Tools for creating and manipulating triangular meshes
 * Copyright (C) 2004-15  Mark J. Stock
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 *
 *********************************************************** */


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/* define a C preprocessor variable so that when structs.h is included,
 * it will contain extra information used only by this program */
#define MODULE_ROCKXRAY
#include "structs.h"

node_ptr node_head = NULL;
norm_ptr norm_head = NULL;
text_ptr text_head = NULL;

static VEC six_views[] = {
  {0.0000000e+00,-5.2573100e-01,8.5065100e-01},
  {-8.5065100e-01,0.0000000e+00,5.2573100e-01},
  {-5.2573100e-01,-8.5065100e-01,0.0000000e+00},
  {8.5065100e-01,0.0000000e+00,5.2573100e-01},
  {5.2573100e-01,-8.5065100e-01,0.0000000e+00},
  {0.0000000e+00,-5.2573100e-01,-8.5065100e-01} };

static VEC nineteen_views[] = {
  {5.2573095e-01,8.5065091e-01,0.0000000e+00},
  {-5.2573095e-01,8.5065091e-01,0.0000000e+00},
  {0.0000000e+00,-5.2573095e-01,8.5065091e-01},
  {-4.8514077e-01,-3.0380062e-01,8.1996562e-01},
  {-2.8955865e-01,-8.2736144e-01,4.8127834e-01},
  {-8.5065091e-01,0.0000000e+00,5.2573095e-01},
  {-8.0715690e-01,-5.0482313e-01,3.0602508e-01},
  {8.5065091e-01,0.0000000e+00,5.2573095e-01},
  {7.9454876e-01,-5.0682845e-01,3.3439077e-01},
  {-4.6526244e-01,3.3434290e-01,8.1960093e-01},
  {-7.9904235e-01,5.3003109e-01,2.8389852e-01},
  {0.0000000e+00,5.2573095e-01,8.5065091e-01},
  {-3.2024080e-01,8.2031953e-01,4.7383721e-01},
  {2.9830821e-01,-8.0409137e-01,5.1424632e-01},
  {5.3825546e-01,-2.6550337e-01,7.9986813e-01},
  {7.8571118e-01,5.6242547e-01,2.5755686e-01},
  {4.9016516e-01,3.4632091e-01,7.9987496e-01},
  {3.0426116e-01,8.3225059e-01,4.6344806e-01},
  {1.1739148e-02,1.5749720e-02,9.9980705e-01} };

static VEC seventysix_views[] = {
  {1.1739148e-02,1.5749720e-02,9.9980705e-01},
  {2.8862360e-01,-1.3184102e-01,9.4832187e-01},
  {2.6342357e-01,1.8860793e-01,9.4606293e-01},
  {5.3825546e-01,-2.6550337e-01,7.9986813e-01},
  {5.3994539e-01,4.3141774e-02,8.4059370e-01},
  {4.9016516e-01,3.4632091e-01,7.9987496e-01},
  {2.8290323e-01,-4.1572213e-01,8.6437311e-01},
  {6.5021116e-03,-2.6519152e-01,9.6417383e-01},
  {0.0000000e+00,-5.2573095e-01,8.5065091e-01},
  {2.5099641e-01,4.5522802e-01,8.5426474e-01},
  {0.0000000e+00,5.2573095e-01,8.5065091e-01},
  {5.6367603e-03,2.7914556e-01,9.6023225e-01},
  {8.5065091e-01,0.0000000e+00,5.2573095e-01},
  {6.9996735e-01,1.8655229e-01,6.8937940e-01},
  {7.2118494e-01,-1.3201849e-01,6.8004662e-01},
  {-8.0715690e-01,-5.0482313e-01,3.0602508e-01},
  {-9.5208302e-01,-2.6750756e-01,1.4824856e-01},
  {-6.8771854e-01,-7.0930616e-01,1.5468672e-01},
  {-8.5065091e-01,0.0000000e+00,5.2573095e-01},
  {-9.6577859e-01,9.0170044e-04,2.5936634e-01},
  {-8.6214488e-01,-2.5677154e-01,4.3677750e-01},
  {5.2573095e-01,8.5065091e-01,0.0000000e+00},
  {5.2573095e-01,-8.5065091e-01,0.0000000e+00},
  {7.8571118e-01,5.6242547e-01,2.5755686e-01},
  {9.4687193e-01,2.9953633e-01,1.1709623e-01},
  {6.7017257e-01,7.3133035e-01,1.2658845e-01},
  {9.6750783e-01,6.6754047e-03,2.5275290e-01},
  {8.6167278e-01,2.9730842e-01,4.1125142e-01},
  {-2.3818418e-01,1.8230039e-01,9.5395748e-01},
  {-2.4833068e-01,-1.5100953e-01,9.5683227e-01},
  {-4.6526244e-01,3.3434290e-01,8.1960093e-01},
  {-5.0193439e-01,1.3063874e-02,8.6480703e-01},
  {-4.8514077e-01,-3.0380062e-01,8.1996562e-01},
  {-2.4142485e-01,4.4576144e-01,8.6198073e-01},
  {-2.4898458e-01,-4.3140358e-01,8.6712031e-01},
  {-6.9542478e-01,-1.5615940e-01,7.0142613e-01},
  {-6.8885233e-01,1.7534146e-01,7.0337603e-01},
  {1.5330990e-01,-9.5085614e-01,2.6901427e-01},
  {-1.5431035e-01,-9.5589386e-01,2.4991048e-01},
  {2.9830821e-01,-8.0409137e-01,5.1424632e-01},
  {2.0151368e-03,-8.5370638e-01,5.2075077e-01},
  {-2.8955865e-01,-8.2736144e-01,4.8127834e-01},
  {4.2836046e-01,-8.6210984e-01,2.7069157e-01},
  {-4.2237904e-01,-8.7258784e-01,2.4532918e-01},
  {-1.5129978e-01,-7.0387758e-01,6.9402069e-01},
  {1.5503975e-01,-6.8899479e-01,7.0798931e-01},
  {-3.2024080e-01,8.2031953e-01,4.7383721e-01},
  {-8.6180943e-03,8.6951572e-01,4.9383008e-01},
  {-1.7694117e-01,9.5277779e-01,2.4679203e-01},
  {3.0426116e-01,8.3225059e-01,4.6344806e-01},
  {1.5036793e-01,9.5847994e-01,2.4229255e-01},
  {1.5951474e-01,7.1121081e-01,6.8464169e-01},
  {-1.6753602e-01,7.0409521e-01,6.9005915e-01},
  {-4.3644630e-01,8.6486491e-01,2.4803896e-01},
  {4.2826329e-01,8.7300171e-01,2.3336358e-01},
  {4.1703103e-01,6.2074539e-01,6.6389780e-01},
  {5.7359323e-01,7.2711282e-01,3.7722374e-01},
  {6.7566115e-01,4.8044570e-01,5.5915466e-01},
  {4.4210096e-01,-5.6574887e-01,6.9604235e-01},
  {6.9754085e-01,-4.0526050e-01,5.9093206e-01},
  {5.7480746e-01,-6.8660734e-01,4.4515924e-01},
  {7.9454876e-01,-5.0682845e-01,3.3439077e-01},
  {8.5774459e-01,-2.5300626e-01,4.4750648e-01},
  {6.8403117e-01,-7.0806653e-01,1.7533726e-01},
  {-9.4972581e-01,2.8147391e-01,1.3708873e-01},
  {-7.9904235e-01,5.3003109e-01,2.8389852e-01},
  {-6.8301575e-01,7.1557332e-01,1.4643874e-01},
  {-8.6044540e-01,2.8157860e-01,4.2467305e-01},
  {-5.8381854e-01,7.0771419e-01,3.9786497e-01},
  {-4.1333531e-01,6.0638680e-01,6.7930035e-01},
  {-6.6981992e-01,4.5738784e-01,5.8492533e-01},
  {8.3202046e-01,-5.5439683e-01,1.9649561e-02},
  {9.5172571e-01,-2.6305908e-01,1.5817110e-01},
  {-5.7968354e-01,-7.0168269e-01,4.1425643e-01},
  {-6.7940909e-01,-4.2677217e-01,5.9688257e-01},
  {-4.1250361e-01,-5.9563265e-01,6.8924780e-01} };

extern int write_xray(tri_pointer,VEC,double*,double*,double*,int,double,int,double,int,double,double,int,int,int,int,char*,char*,int);
int Usage(char[255],int);

int main(int argc,char **argv) {

   int i,do_volume,do_fade,max_size,force_square,quality,write_hibit;
   int force_num_threads = -1;
   int num_layers;				// how many layers to render to?
   int do6 = FALSE;
   int do19 = FALSE;
   int do76 = FALSE;
   char infile[255];				/* name of input file */
   char progname[255];				/* name of binary executable */
   char output_format[4];			/* file format extension for output */
   char* out_prefix = NULL;			/* prefix to use for output files */
   double thickness;				/* thickness of mesh, world coords */
   double border;				/* thickness of image border, fraction */
   double peak_crop;				/* muliplier on image value to crop peaks */
   double gamma;				/* output image gamma value */
   double xb[3],yb[3],zb[3];			/* image bounds, in world units, [t/f,min,max] */
   VEC viewp;
   tri_pointer tri_head = NULL;

   viewp.x = 0.0;
   viewp.y = -1.0;
   viewp.z = 0.0;
   max_size = 512;				// default is sane, now
   force_square = FALSE;
   num_layers = 1;
   quality = 0;					// 0 is default, 1 higher, 2 very high
   write_hibit = FALSE;
   do_volume = FALSE;
   do_fade = FALSE;
   border = 0.1;
   peak_crop = 0.8;
   gamma = 1.0;
   thickness = -1.0;
   xb[0] = -1;					// negative means "do not use bounds"
   xb[1] = 0;
   xb[2] = 0;
   yb[0] = -1;					// negative means "do not use bounds"
   yb[1] = 0;
   yb[2] = 0;
   zb[0] = -1;					// negative means "do not use bounds"
   zb[1] = 0;
   zb[2] = 0;

   /* Parse command-line args */
   (void) strcpy(progname,argv[0]);
   if (argc < 2) {
      fprintf(stderr,"\nFirst argument must be input file.\n");
      fflush(stderr);
      (void) Usage(progname,0);
   }
   if (strncmp(argv[1], "-", 1) == 0) {
      fprintf(stderr,"\nFirst argument must be input file.\n");
      fflush(stderr);
      (void) Usage(progname,0);
   }
   (void) strcpy(infile,argv[1]);
   for (i=2; i<argc; i++) {
      if (strncmp(argv[i], "-xb", 3) == 0) {
         xb[0] = +1.0;
         xb[1] = atof(argv[++i]);
         xb[2] = atof(argv[++i]);
      } else if (strncmp(argv[i], "-yb", 3) == 0) {
         yb[0] = +1.0;
         yb[1] = atof(argv[++i]);
         yb[2] = atof(argv[++i]);
      } else if (strncmp(argv[i], "-zb", 3) == 0) {
         zb[0] = +1.0;
         zb[1] = atof(argv[++i]);
         zb[2] = atof(argv[++i]);
      } else if (strncmp(argv[i], "-s", 2) == 0) {
         do_volume = FALSE;
      } else if (strncmp(argv[i], "-fade", 3) == 0) {
         do_fade = TRUE;
      } else if (strncmp(argv[i], "-f", 2) == 0) {
         force_square = TRUE;
      } else if (strncmp(argv[i], "-v", 2) == 0) {
         do_volume = TRUE;
      } else if (strncmp(argv[i], "-qqq", 4) == 0) {
         quality = 3;
      } else if (strncmp(argv[i], "-qq", 3) == 0) {
         quality = 2;
      } else if (strncmp(argv[i], "-q", 2) == 0) {
         quality = 1;
      } else if (strncmp(argv[i], "-b", 2) == 0) {
         border = atof(argv[++i]);
      } else if (strncmp(argv[i], "-pc", 3) == 0) {
         peak_crop = atof(argv[++i]);
      } else if (strncmp(argv[i], "-g", 2) == 0) {
         gamma = atof(argv[++i]);
      } else if (strncmp(argv[i], "-t", 2) == 0) {
         thickness = atof(argv[++i]);
      } else if (strncmp(argv[i], "-r", 2) == 0) {
         max_size = atoi(argv[++i]);
      } else if (strncmp(argv[i], "-8", 2) == 0) {
         write_hibit = FALSE;
      } else if (strncmp(argv[i], "-16", 3) == 0) {
         write_hibit = TRUE;
      } else if (strncmp(argv[i], "-do6", 4) == 0) {
         do6 = TRUE;
      } else if (strncmp(argv[i], "-do19", 5) == 0) {
         do19 = TRUE;
      } else if (strncmp(argv[i], "-do76", 5) == 0) {
         do76 = TRUE;
      } else if (strncmp(argv[i], "-d", 2) == 0) {
         viewp.x = atof(argv[++i]);
         viewp.y = atof(argv[++i]);
         viewp.z = atof(argv[++i]);
      } else if (strncmp(argv[i], "-o", 2) == 0) {
         strncpy(output_format,argv[i]+2,4);
      } else if (strncmp(argv[i], "-layers", 2) == 0) {
         num_layers = atoi(argv[++i]);
      } else if (strncmp(argv[i], "-prefix", 3) == 0) {
         out_prefix = (char*)malloc(255*sizeof(char));
         strcpy(out_prefix,argv[++i]);
      } else if (strncmp(argv[i], "-n", 2) == 0) {
         force_num_threads = atoi(argv[++i]);
      } else {
         fprintf(stderr,"\nUnrecognized argument (%s)\n",argv[i]);
         fflush(stderr);
         (void) Usage(progname,0);
      }
   }
   if (max_size > MAX_IMAGE) {
      fprintf(stderr,"Maximum image size is %d, reducing.\n",MAX_IMAGE);
      max_size = MAX_IMAGE;
   }

   /* If asking for multiple images, and no prefix is present, assign a default here */
   if (!out_prefix && (num_layers > 1 || do6 || do19 || do76)) {
      out_prefix = (char*)malloc(255*sizeof(char));
      strcpy(out_prefix,"out");
   }

   /* Read the input file */
   tri_head = read_input(infile,FALSE,NULL);

   if (do6) {
      for (int i=0; i<6; i++) {
         // append an index to the output prefix
         char new_prefix[255];
         sprintf(new_prefix, "%s_v%02d", out_prefix, i);
         fprintf(stderr,"\nRendering number %d of %d to %s\n", i, 6, new_prefix);
         // set the viewpoint
         viewp = six_views[i];
         // render the image
         (void) write_xray(tri_head,viewp,xb,yb,zb,max_size,thickness,force_square,
                           border,quality,peak_crop,gamma,write_hibit,do_volume,do_fade,
                           num_layers,new_prefix,output_format,force_num_threads);
      }

   } else if (do19) {
      for (int i=0; i<19; i++) {
         // append an index to the output prefix
         char new_prefix[255];
         sprintf(new_prefix, "%s_v%02d", out_prefix, i);
         fprintf(stderr,"\nRendering number %d of %d to %s\n", i, 19, new_prefix);
         // set the viewpoint
         viewp = nineteen_views[i];
         // render the image
         (void) write_xray(tri_head,viewp,xb,yb,zb,max_size,thickness,force_square,
                           border,quality,peak_crop,gamma,write_hibit,do_volume,do_fade,
                           num_layers,new_prefix,output_format,force_num_threads);
      }

   } else if (do76) {
      for (int i=0; i<76; i++) {
         // append an index to the output prefix
         char new_prefix[255];
         sprintf(new_prefix, "%s_v%02d", out_prefix, i);
         fprintf(stderr,"\nRendering number %d of %d to %s\n", i, 76, new_prefix);
         // set the viewpoint
         viewp = seventysix_views[i];
         // render the image
         (void) write_xray(tri_head,viewp,xb,yb,zb,max_size,thickness,force_square,
                           border,quality,peak_crop,gamma,write_hibit,do_volume,do_fade,
                           num_layers,new_prefix,output_format,force_num_threads);
      }

   } else {
      /* Just write one image to stdout */
      (void) write_xray(tri_head,viewp,xb,yb,zb,max_size,thickness,force_square,
                        border,quality,peak_crop,gamma,write_hibit,do_volume,do_fade,
                        num_layers,out_prefix,output_format,force_num_threads);
   }

   fprintf(stderr,"Done.\n");
   exit(0);
}


/*
 * This function writes basic usage information to stderr,
 * and then quits. Too bad.
 */
int Usage(char progname[255],int status) {

   /* Usage for rockxray */
   static char **cpp, *help_message[] =
   {
       "where [-options] are one or more of the following:                         ",
       "                                                                           ",
       "   -d [x y z]  specify the view direction, object will always be centered; ",
       "               view up vector will always be in the +z direction           ",
       "                                                                           ",
       "   -xb [lo hi] specify x-image bounds in world units, may not work well    ",
       "                                                                           ",
       "   -yb [lo hi] specify y-image bounds in world units, may not work well    ",
       "                                                                           ",
       "   -zb [lo hi] specify z-image bounds in world units, may not work well    ",
       "                                                                           ",
       "   -layers n   write multiple z-layers (does NOT write to stdout)          ",
       "                                                                           ",
       "   -prefix p   prefix every z-layer with this ('out' writes 'out_00.png')  ",
       "                                                                           ",
       "   -f          force output image to be a square                           ",
       "                                                                           ",
       "   -b frac     size of border around geometry, as fraction of image size,  ",
       "               default=0.1                                                 ",
       "                                                                           ",
       "   -q          create higher quality image at the expense of compute time  ",
       "                                                                           ",
       "   -qq         creates even higher quality image                           ",
       "                                                                           ",
       "   -qqq        creates highest quality image---takes very long!            ",
       "                                                                           ",
       "   -t num      thickness of the mesh in world coordinates,                 ",
       "               default = (1 layer), regardless of resolution               ",
       "                                                                           ",
       "   -s          image only the shell of the mesh (default behavior)         ",
       "                                                                           ",
       "   -v          image the volume of the mesh                                ",
       "                                                                           ",
       "   -fade       fade intensity from front to back                           ",
       "                                                                           ",
       "   -8          write an 8-bit grey image (default)                         ",
       "                                                                           ",
       "   -16         write a 16-bit grey image                                   ",
       "                                                                           ",
       "   -pc frac    peak cropping: peak value in image will be this number times",
       "               the actual peak value computed, default=0.8, 1.0 means no   ",
       "               data is lost, but a smaller number may produce a better     ",
       "               image                                                       ",
       "                                                                           ",
       "   -g gamma    gamma value (exponent) for output image                     ",
       "                                                                           ",
       "   -r [res]    pixel resolution of the long edge of the image, default=512 ",
       "                                                                           ",
       "   -n num      force this number of threads (default is number of cores)   ",
       "                                                                           ",
       "   -okey       specify output format, key= pgm, png, default = png         ",
       "                                                                           ",
       "   -help       returns this help information                               ",
       " ",
       "The input file can be of .obj, .raw, or .tin format, and the program requires",
       "   the input file to use its valid 3-character filename extension.",
       " ",
       "Options may be abbreviated to an unambiguous length (duh).",
       "Output is to stdout",
       NULL
   };

   fprintf(stderr, "\nusage:\n  %s infile [-options]\n\n", progname);
   for (cpp = help_message; *cpp; cpp++)
      fprintf(stderr, "%s\n", *cpp);
      fflush(stderr);
   exit(status);
   return(0);
}
