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

extern int write_xray(tri_pointer,VEC,double*,double*,int,double,int,double,int,double,double,int,int,char*,int);
int Usage(char[80],int);

int main(int argc,char **argv) {

   int i,do_volume,max_size,force_square,quality,write_hibit;
   int force_num_threads = -1;
   char infile[80];				/* name of input file */
   char progname[80];				/* name of binary executable */
   char output_format[4];			/* file format extension for output */
   double thickness;				/* thickness of mesh, world coords */
   double border;				/* thickness of image border, fraction */
   double peak_crop;				/* muliplier on image value to crop peaks */
   double gamma;				/* output image gamma value */
   double xb[3],yb[3];				/* image bounds, in world units, [t/f,min,max] */
   VEC viewp;
   tri_pointer tri_head = NULL;

   viewp.x = 0.0;
   viewp.y = -1.0;
   viewp.z = 0.0;
   max_size = 512;				// default is sane, now
   force_square = FALSE;
   quality = 0;					// 0 is default, 1 higher, 2 very high
   write_hibit = FALSE;
   do_volume = FALSE;
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

   /* Parse command-line args */
   (void) strcpy(progname,argv[0]);
   if (argc < 2) (void) Usage(progname,0);
   if (strncmp(argv[1], "-", 1) == 0)
      (void) Usage(progname,0);
   (void) strcpy(infile,argv[1]);
   for (i=2; i<argc; i++) {
      if (strncmp(argv[i], "-d", 2) == 0) {
         viewp.x = atof(argv[++i]);
         viewp.y = atof(argv[++i]);
         viewp.z = atof(argv[++i]);
      } else if (strncmp(argv[i], "-xb", 3) == 0) {
         xb[0] = +1.0;
         xb[1] = atof(argv[++i]);
         xb[2] = atof(argv[++i]);
      } else if (strncmp(argv[i], "-yb", 3) == 0) {
         yb[0] = +1.0;
         yb[1] = atof(argv[++i]);
         yb[2] = atof(argv[++i]);
      } else if (strncmp(argv[i], "-s", 2) == 0) {
         do_volume = FALSE;
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
      } else if (strncmp(argv[i], "-o", 2) == 0) {
         strncpy(output_format,argv[i]+2,4);
      } else if (strncmp(argv[i], "-n", 2) == 0) {
         force_num_threads = atoi(argv[++i]);
      } else
         (void) Usage(progname,0);
   }
   if (max_size > MAX_IMAGE) {
      fprintf(stderr,"Maximum image size is %d, reducing.\n",MAX_IMAGE);
      max_size = MAX_IMAGE;
   }

   /* Read the input file */
   tri_head = read_input(infile,FALSE,NULL);

   /* Write the image to stdout */
   (void) write_xray(tri_head,viewp,xb,yb,max_size,thickness,force_square,
                     border,quality,peak_crop,gamma,write_hibit,do_volume,
                     output_format,force_num_threads);

   fprintf(stderr,"Done.\n");
   exit(0);
}


/*
 * This function writes basic usage information to stderr,
 * and then quits. Too bad.
 */
int Usage(char progname[80],int status) {

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
       "   -f          force output image to be a square                           ",
       "                                                                           ",
       "   -b frac     size of border around geometry, as fraction of image size,  ",
       "               default=0.1                                                 ",
       "                                                                           ",
       "   -q          create higher quality image at the expense of compute time  ",
       "                                                                           ",
       "   -t num      thickness of the mesh in world coordinates,                 ",
       "               default = (1 layer), regardless of resolution               ",
       "                                                                           ",
       "   -s          image only the shell of the mesh (default behavior)         ",
       "                                                                           ",
       "   -v          image the volume of the mesh                                ",
       "                                                                           ",
       "   -8          write an 8-bit grey image                                   ",
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
       "   -okey       specify output format, key= pgm, png, default = png         ",
       "                                                                           ",
       "   -n num      force this number of threads                                ",
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

   fprintf(stderr, "usage:\n  %s infile [-options]\n\n", progname);
   for (cpp = help_message; *cpp; cpp++)
      fprintf(stderr, "%s\n", *cpp);
      fflush(stderr);
   exit(status);
   return(0);
}
