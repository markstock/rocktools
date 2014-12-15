/*************************************************************
 *
 *  rockpng.c - Read a PNG heightfield, write an OBJ mesh
 *
 *  Mark J. Stock, mstock@umich.edu
 *
 *
 * rocktools - Tools for creating and manipulating triangular meshes
 * Copyright (C) 1999,2003-2004,2006-2007,2010-14  Mark J. Stock
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
#include <ctype.h>
#include "structs.h"

node_ptr node_head = NULL;
norm_ptr norm_head = NULL;
text_ptr text_head = NULL;

int num_tri = 0;

void rescale_nodes (double);
int Usage(char[160],int);

extern float** read_png (char*, float, float, int*, int*);
extern tri_pointer generate_heightmesh (tri_pointer, float**, int, int, int, int, double, int, int, double, double, int);

int main(int argc,char **argv) {

   int i,j;
   char infile[160];		/* name of input file */
   char progname[160];		/* name of binary executable */
   char output_format[4];	/* file format extension for output */
   tri_pointer tri_head = NULL;
   int nx = -1;
   int ny = -1;
   int do_rescale = FALSE;
   int do_bottom = FALSE;
   int do_trans = FALSE;
   int do_legs = FALSE;
   int do_walls = FALSE;
   int do_texture_coords = FALSE;
   float **hf = NULL;
   double depth = 0.01;
   double thick = 0.01;
   double inset = 0.2;
   double scale = 1.;
   double hmin = 0.;
   double hmax = 1.;
   double initmin,initmax = 1.;


   /* Parse command-line args */
   (void) strcpy(progname,argv[0]);
   if (argc < 2) (void) Usage(progname,0);
   if (strncmp(argv[1], "-help", 2) == 0)
      (void) Usage(progname,0);
   (void) strcpy(infile,argv[1]);
   for (i=2; i<argc; i++) {
      if (strncmp(argv[i], "-o", 2) == 0) {
         // specify output format (Shapeways takes obj now)
         strncpy(output_format,argv[i]+2,4);
      } else if (strncmp(argv[i], "-finalscale", 2) == 0) {
         // after all geom is created, scale up to millimeters
         if (argc > i) {
            if (!isalpha(argv[i+1][1])) {
               do_rescale = TRUE;
               scale = atof(argv[++i]);
            }
         }
      } else if (strncmp(argv[i], "-bottom", 2) == 0) {
         do_bottom = TRUE;
      } else if (strncmp(argv[i], "-transparency", 3) == 0) {
         do_bottom = TRUE;
         do_trans = TRUE;
      } else if (strncmp(argv[i], "-legs", 2) == 0) {
         do_legs = TRUE;
         // legs and walls are exclusive
         do_walls = FALSE;
         // always do a bottom if we need legs
         do_bottom = TRUE;
      } else if (strncmp(argv[i], "-walls", 2) == 0) {
         do_walls = TRUE;
         // legs and walls are exclusive
         do_legs = FALSE;
         // always do a bottom if we need walls
         do_bottom = TRUE;
      } else if (strncmp(argv[i], "-elevation", 2) == 0) {
         if (argc > i) {
            if (!isalpha(argv[i+1][1])) {
               // first number after elevation is minimum height above ground
               hmin = atof(argv[++i]);
            }
         }
         if (argc > i) {
            if (!isalpha(argv[i+1][1])) {
               hmax = atof(argv[++i]);
            } else {
               // if only one number after elevation, it's max height
               hmax = hmin;
               hmin = 0.0;
            }
         }
      } else if (strncmp(argv[i], "-depth", 2) == 0) {
         if (argc > i) {
            if (!isalpha(argv[i+1][1])) {
               depth = atof(argv[++i]);
            }
         }
      } else if (strncmp(argv[i], "-thickness", 3) == 0) {
         if (argc > i) {
            if (!isalpha(argv[i+1][1])) {
               thick = atof(argv[++i]);
            }
         }
      } else if (strncmp(argv[i], "-inset", 2) == 0) {
         if (argc > i) {
            if (!isalpha(argv[i+1][1])) {
               inset = atof(argv[++i]);
            }
         }
      } else if (strncmp(argv[i], "-texture", 3) == 0) {
         do_texture_coords = TRUE;
      } else {
         (void) Usage(progname,0);
      }
   }


   // Read the input PNG, initially scale between 0 and 1
   hf = read_png (infile, 0., 1., &nx, &ny);
   fprintf(stderr,"Read %d x %d png file\n",nx,ny);
   fflush(stderr);

   // if we're doing two-sided, scale the intensities to a log scale
   //   so that light traveling through the piece is correct
   // later

   // if we're doing a bottom, and the layer thickness would touch the ground, fix it
   if (do_bottom && !do_trans) {
      if (hmin - depth < 0.0) {
         // move the hmax up the same amount
         hmax = 2*depth + (hmax-hmin);
         // and reset the hmin to leave a gap the size of depth
         hmin = 2*depth;
         fprintf(stderr,"Reset -elev to %g %g to account for depth\n",hmin,hmax);
         fflush(stderr);
      }
   }

   // if we're doing a transparency, scale the heightfield to logarithms
   // this thickness will be doubled on the underside
   if (do_trans) {
      for (i=0; i<nx; i++) {
         for (j=0; j<ny; j++) {
            //hf[i][j] = (double)log(0.05+hf[i][j]);
            // add a tad to make sure we don't get zero
            hf[i][j] += 0.03;
            // ungamma-correct
            //hf[i][j] = expf(0.45*logf(hf[i][j]));
            // scale value to thickness
            hf[i][j] = logf(1./hf[i][j]);
            // min thickness will be applied later
         }
      }
   }

   // convert PNG to proper height scale
   // find data max and min
   initmax = 0.0;
   initmin = 1.0;
   for (i=0; i<nx; i++) {
      for (j=0; j<ny; j++) {
         if (hf[i][j] > initmax) initmax = (double)hf[i][j];
         if (hf[i][j] < initmin) initmin = (double)hf[i][j];
      }
   }
   fprintf(stderr,"Found min/max of %g / %g\n",initmin,initmax);
   fprintf(stderr,"Scaling height to %g : %g\n",hmin,hmax);
   fflush(stderr);
   // now scale
   for (i=0; i<nx; i++) {
      for (j=0; j<ny; j++) {
         hf[i][j] = (float)(hmin + ((double)hf[i][j]-initmin)*(hmax-hmin)/(initmax-initmin));
      }
   }

   // generate the mesh
   fprintf(stderr,"Generating mesh\n");
   fflush(stderr);
   tri_head = generate_heightmesh (tri_head,hf,nx,ny,do_bottom,do_trans,depth,
                                   do_legs,do_walls,thick,inset,do_texture_coords);

   // apply uniform scaling transformation
   if (do_rescale) {
      fprintf(stderr,"Scaled mesh by %g\n",scale);
      fflush(stderr);
      rescale_nodes (scale);
   }

   // Write triangles to stdout
   (void) write_output (tri_head,output_format,TRUE,argc,argv);

   fprintf(stderr,"Done.\n");
   exit(0);
}


/*
 * Scale all nodes
 */
void rescale_nodes (double scale) {

   node_ptr this_node;

   // then, perturb each according to the normal
   this_node = node_head;
   while (this_node) {
      this_node->loc.x *= scale;
      this_node->loc.y *= scale;
      this_node->loc.z *= scale;
      this_node = this_node->next_node;
   }

   return;
}


/*
 * This function writes basic usage information to stderr,
 * and then quits. Too bad.
 */
int Usage(char progname[80],int status) {

   /* Usage for rockpng */
   static char **cpp, *help_message[] =
   {
       "where [-options] are one or more of the following:                         ",
       "                                                                           ",
       "   -okey       specify output format, key= raw, rad, pov, obj, tin, rib    ",
       "               default = raw; surface normal vectors are not supported     ",
       "                                                                           ",
       "   -elevation minval maxval                                                ",
       "               if maximum horizontal dimension is 1.0, these are the       ",
       "               minimum and maximum extents of the heightfield geometry     ",
       "                                                                           ",
       "   -finalscale val                                                         ",
       "               scale geometry by the given factors, default=1              ",
       "                                                                           ",
       "   -bottom                                                                 ",
       "               generate a lower-facing mesh beneath the top mesh           ",
       "                                                                           ",
       "   -transparency                                                           ",
       "               generate a two-sided mesh, with relief on both sides        ",
       "                                                                           ",
       "   -depth val                                                              ",
       "               depth (in non-dimensional units) of top-to-bottom layer     ",
       "               default is 0.01                                             ",
       "                                                                           ",
       "   -legs                                                                   ",
       "               include geometry for legs                                   ",
       "                                                                           ",
       "   -walls                                                                  ",
       "               include geometry for curtain walls around base of model     ",
       "                                                                           ",
       "   -thick val                                                              ",
       "               thickness (in non-dimensional units) of walls or legs       ",
       "               default is 0.01                                             ",
       "                                                                           ",
       "   -inset val                                                              ",
       "               wall or leg inset, (in non-dimensional units)               ",
       "               default is 0.2                                              ",
       "                                                                           ",
       "   -texture                                                                ",
       "               generate and write texture coordinates for the top surface  ",
       "                                                                           ",
       "   -help       (in place of infile) returns this help information          ",
       " ",
       "The input file can be of .obj, .raw, .msh, or .tin format, and the program",
       "   requires the input file to use its valid 3-character filename extension ",
       " ",
       "Options may be abbreviated to an unambiguous length",
       "Output is to stdout, so redirect it to a file using '>'",
       " ",
       "rockpng creates a triangle mesh from a PNG heightfield                     ",
       "                                                                           ",
       "examples:                                                                  ",
       "   rockpng in.tin -oobj > out.obj                                      ",
       "      Read in the file in.tin and write a Wavefront .obj format version    ",
       "                                                                           ",
       "   rockpng in.raw -s 2 -orad > out.rad                                 ",
       "      Read in the input file and scale all nodes by 2.0 in x, y, and z;    ",
       "      finally write the result in Radiance format                          ",
       "                                                                           ",
       NULL
   };

   fprintf(stderr, "usage:\n  %s infile [-options]\n\n", progname);
   for (cpp = help_message; *cpp; cpp++)
      fprintf(stderr, "%s\n", *cpp);
      fflush(stderr);
   exit(status);
   return(0);
}
