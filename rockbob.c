/*************************************************************
 *
 *  rockbob.c - Create a brick of bytes or floats from a mesh
 *
 *  Mark J. Stock, mstock@umich.edu
 *
 * rocktools - Tools for creating and manipulating triangular meshes
 * Copyright (C) 2004-14  Mark J. Stock
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
#define MODULE_ROCKBOB
#include "structs.h"

node_ptr node_head = NULL;
norm_ptr norm_head = NULL;
text_ptr text_head = NULL;

extern int write_bob(tri_pointer,double*,double*,double*,double,double,int,double,double,char*);
int Usage(char[MAX_FN_LEN],int);

int main(int argc,char **argv) {

   int i;
   char infile[MAX_FN_LEN];				/* name of input file */
   char progname[MAX_FN_LEN];				/* name of binary executable */
   char output_format[4];			/* file format extension for output */
   int diffuseSteps = 0;			/* number of steps to diffuse */
   double dx;					/* voxel size */
   double thickness;				/* thickness of mesh, world coords */
   double repose;				/* angle of repose (45-90), negative turns off */
   double erode;				/* number of cells to erode the volume, neg is dilate */
   double xb[3],yb[3],zb[3];			/* bounds, in world units, [t/f,min,max] */
   tri_pointer tri_head = NULL;

   // negative means "not being used"
   dx = -1.0;
   thickness = -1.0;
   repose = -45.0;
   erode = 0.0;
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
   if (argc < 2) (void) Usage(progname,0);
   if (strncmp(argv[1], "-", 1) == 0)
      (void) Usage(progname,0);
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
      } else if (strncmp(argv[i], "-dx", 3) == 0) {
         dx = atof(argv[++i]);
         if (dx <= 0.0) {
            (void) Usage(progname,0);
            fprintf(stderr,"Error: -d must be positive\n");
         }
      } else if (strncmp(argv[i], "-t", 2) == 0) {
         thickness = atof(argv[++i]);
         if (thickness <= 0.0) {
            (void) Usage(progname,0);
            fprintf(stderr,"Error: -t must be positive\n");
         }
      } else if (strncmp(argv[i], "-diffuse", 3) == 0) {
         diffuseSteps = atoi(argv[++i]);
      } else if (strncmp(argv[i], "-repose", 3) == 0) {
         repose = atof(argv[++i]);
      } else if (strncmp(argv[i], "-erode", 3) == 0) {
         erode = atof(argv[++i]);
      } else if (strncmp(argv[i], "-o", 2) == 0) {
         strncpy(output_format,argv[i]+2,3);
      } else
         (void) Usage(progname,0);
   }

   /* Read the input file */
   tri_head = read_input(infile,FALSE,NULL);

   /* Write the image to stdout */
   (void) write_bob(tri_head,xb,yb,zb,dx,thickness,diffuseSteps,repose,erode,output_format);

   fprintf(stderr,"Done.\n");
   exit(0);
}


/*
 * This function writes basic usage information to stderr,
 * and then quits. Too bad.
 */
int Usage(char progname[MAX_FN_LEN],int status) {

   /* Usage for rockbob */
   static char **cpp, *help_message[] =
   {
       "where [-options] are one or more of the following:                         ",
       "                                                                           ",
       "   -xb [lo hi] specify x-image bounds in world units                       ",
       "                                                                           ",
       "   -yb [lo hi] specify y-image bounds in world units                       ",
       "                                                                           ",
       "   -zb [lo hi] specify z-image bounds in world units                       ",
       "                                                                           ",
       "   -t num      thickness of the mesh in world coordinates,                 ",
       "               default is two voxels                                       ",
       "                                                                           ",
       "   -dx num     size of voxel in world coordinates,                         ",
       "               default is 1/100th of largest dimension                     ",
       "                                                                           ",
       "   -diffuse num                                                            ",
       "               perform num diffusion/smoothing iterations (default=0)      ",
       "                                                                           ",
       "   -repose angle                                                           ",
       "               create volume supports under object (45-90) (default=none)  ",
       "                                                                           ",
       "   -erode distance                                                         ",
       "               number of cells to erode (shrink) the volume, negative      ",
       "               will dilate (grow) instead (default=0)                      ",
       "                                                                           ",
       "   -okey       specify output format, key= bob, bof, default = bob         ",
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
   for (cpp = help_message; *cpp; cpp++) fprintf(stderr, "%s\n", *cpp);
   fflush(stderr);
   exit(status);
   return(0);
}
