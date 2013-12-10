/*************************************************************
 *
 *  rockmarker.c - Write a marker object in place of each
 *	triangle of an irregular triangle mesh
 *
 *  Mark J. Stock, mstock@umich.edu
 *
 *
 * rocktools - Tools for creating and manipulating triangular meshes
 * Copyright (C) 1999,2003-2006  Mark J. Stock
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
#include "structs.h"


node_ptr node_head = NULL;

int num_tri = 0;

int write_markers(tri_pointer, MARKER, char[4]);
int Usage(char[80],int);


int main(int argc,char **argv) {

   int i;
   int num_wrote;
   char infile[80];				/* name of input file */
   char progname[80];				/* name of binary executable */
   char output_format[4];			/* file format extension for output */
   tri_pointer tri_head = NULL;
   MARKER marker;

   // set defaults
   marker.marker_type = sphere;
   marker.marker_size = 1.;
   marker.randomize_size = FALSE;
   marker.min_size = 0.6;
   marker.max_size = 1.6;
   marker.size_range = marker.max_size - marker.min_size;
   marker.density_type = by_elem;
   marker.density = 1.;
   marker.randomize_radius = FALSE;
   marker.min_radius = 0.;
   marker.max_radius = 1.;
   marker.radius_range = marker.max_radius - marker.min_radius;
   marker.randomize_height = FALSE;
   marker.min_height = 0.0;
   marker.max_height = 1.0;
   marker.height_range = marker.max_height - marker.min_height;
   marker.randomize_rotation = FALSE;
   marker.randomize_normal = FALSE;
   marker.normal_pert = 1.0;

   /* Parse command-line args */
   (void) strcpy(progname,argv[0]);
   if (argc < 2) (void) Usage(progname,0);
   if (strncmp(argv[1], "-help", 2) == 0)
      (void) Usage(progname,0);
   (void) strcpy(infile,argv[1]);
   for (i=2; i<argc; i++) {
      if (strncmp(argv[i], "-o", 2) == 0) {
         strncpy(output_format,argv[i]+2,4);
      } else if (strncmp(argv[i], "-de", 5) == 0) {
         marker.density_type = by_elem;
         if (i < argc-1) {
            if (strncmp(argv[i+1], "-", 1) != 0) {
                marker.density = atof(argv[++i]);
            }
         }
      } else if (strncmp(argv[i], "-da", 5) == 0) {
         marker.density_type = by_area;
         if (i < argc-1) {
            if (strncmp(argv[i+1], "-", 1) != 0) {
               marker.density = atof(argv[++i]);
            }
         }
      } else if (strncmp(argv[i], "-rs", 3) == 0) {
         marker.randomize_size = TRUE;
         if (i < argc-2) {
            if (strncmp(argv[i+1], "-", 1) != 0) {
               marker.min_size = atof(argv[++i]);
            }
            if (strncmp(argv[i+1], "-", 1) != 0) {
               marker.max_size = atof(argv[++i]);
            }
         }
      } else if (strncmp(argv[i], "-rp", 3) == 0) {
         marker.randomize_radius = TRUE;
         marker.randomize_height = TRUE;
         if (i < argc-2) {
            if (strncmp(argv[i+1], "-", 1) != 0) {
               marker.max_radius = atof(argv[++i]);
            }
            if (strncmp(argv[i+1], "-", 1) != 0) {
               marker.max_height = atof(argv[++i]);
            }
         }
      } else if (strncmp(argv[i], "-rr", 3) == 0) {
         marker.randomize_rotation = TRUE;
      } else if (strncmp(argv[i], "-rn", 3) == 0) {
         marker.randomize_normal = TRUE;
         if (i < argc-1) {
            if (strncmp(argv[i+1], "-", 1) != 0) {
               marker.normal_pert = atof(argv[++i]);
            }
         }
      } else if (strncmp(argv[i], "-r", 2) == 0) {
         marker.marker_type = rectangle;
         if (i < argc-1) {
            if (strncmp(argv[i+1], "-", 1) != 0) {
               marker.marker_size = atof(argv[++i]);
            }
         }
      } else if (strncmp(argv[i], "-s", 2) == 0) {
         marker.marker_type = sphere;
         if (i < argc-1) {
            if (strncmp(argv[i+1], "-", 1) != 0) {
               marker.marker_size = atof(argv[++i]);
            }
         }
      } else if (strncmp(argv[i], "-c", 2) == 0) {
         marker.marker_type = dualcone;
         if (i < argc-1) {
            if (strncmp(argv[i+1], "-", 1) != 0) {
               marker.marker_size = atof(argv[++i]);
            }
         }
      } else {
         (void) Usage(progname,0);
      }
   }

   // recalcuate ranges
   marker.size_range = marker.max_size - marker.min_size;
   marker.radius_range = marker.max_radius - marker.min_radius;
   marker.height_range = marker.max_height - marker.min_height;

   // Read the input file
   tri_head = read_input(infile,FALSE,NULL);

   // Write markers to stdout
   num_wrote = write_markers(tri_head,marker,output_format);
   fprintf(stderr,"Wrote %d markers.\n",num_wrote);

   fprintf(stderr,"Done.\n");
   exit(0);
}


/*
 * This function writes basic usage information to stderr,
 * and then quits. Too bad.
 */
int Usage(char progname[80],int status) {

   /* Usage for rockmarker */
   static char **cpp, *help_message[] =
   {
       "where [-options] are one or more of the following:                         ",
       "                                                                           ",
       "   -s [rad]    use spheres as markers, rad is non-dimensional radius with  ",
       "               1 being the size of the element, default = 1.0              ",
       "                                                                           ",
       "   -r [size]   use rectangles as markers, size is non-dimensional size     ",
       "               with 1 being the size of the element, default = 1.0         ",
       "                                                                           ",
       "   -c [size]   use dual-cones as markers, size is non-dimensional size with",
       "               1 being the size of the element, default = 1.0              ",
       "                                                                           ",
       "   -de [dens]  element-wise density, dens is mean number of markers per    ",
       "               element, default is 1.0                                     ",
       "                                                                           ",
       "   -da [dens]  area-wise density, dens is mean number of markers per unit  ",
       "               area, default is 1.0                                        ",
       "                                                                           ",
       "   -rs [small large]  randomize size between given bounds,                 ",
       "               default = 0.6 to 1.6                                        ",
       "                                                                           ",
       "   -rp [planar normal]  randomize positions relative to element center,    ",
       "               radius in plane of element, distance normal to plane,       ",
       "               non-dimensional units, default = 1.0 and 0.0                ",
       "                                                                           ",
       "   -rr         randomize rotation (in-plane)                               ",
       "                                                                           ",
       "   -rn [pert]  randomize normal direction, pert is perturbation scale,     ",
       "               where 0 is no perturbation, and inf is total randomness,    ",
       "               default is 1.0                                              ",
       "                                                                           ",
       "   -okey       specify output format, key= obj, rad                        ",
       "               default = rad; example \"-oobj\" or \"-orad\"               ",
       "                                                                           ",
       "   -help       (in place of infile) returns this help information          ",
       " ",
       "The input file can be of .obj, .raw, .msh, or .tin format, and the program",
       "   requires the input file to use its valid 3-character filename extension.",
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
