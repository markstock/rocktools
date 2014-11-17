/*************************************************************
 *
 *  rockerode.c - Apply fluvial erosion to a tri mesh
 *
 *  Mark J. Stock, mstock@umich.edu
 *
 *
 * rocktools - Tools for creating and manipulating triangular meshes
 * Copyright (C) 1999,2002-03,06,14  Mark J. Stock
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
norm_ptr norm_head = NULL;

int num_tri = 0;

int Usage(char[80],int);
extern int find_flow(tri_pointer);
extern int fill_basins(tri_pointer);
extern int erode_surface(tri_pointer, double);

int main(int argc,char **argv) {

   int i;
   int step = 0;
   int total_steps = 1;			/* total number of erosion steps */
   int do_erosion= 0;			/* perturb nodes to smooth shape? */
   double erosion_factor = 0.0;		/* amount of smoothing to take place */
   //double v_thresh = 1.0;		/* threshhold for common normal, convex edge */
   //double c_thresh = 1.0;		/* threshhold for common normal, concave edge */
   char infile[80];			/* name of input file, TIN only for starters */
   char progname[80];			/* name of binary executable */
   char output_format[4];		/* file format extension for output */
   tri_pointer tri_head = NULL;


   /* Parse command-line args */
   (void) strcpy(progname,argv[0]);
   if (argc < 2) (void) Usage(progname,0);
   if (strncmp(argv[1], "-help", 2) == 0)
      (void) Usage(progname,0);
   (void) strcpy(infile,argv[1]);
   for (i=2; i<argc; i++) {
      if (strncmp(argv[i], "-e", 2) == 0) {
         do_erosion = 1;
         erosion_factor = 1.0;
         if (i < argc-1)
            if (strncmp(argv[i+1], "-", 1) != 0)
               erosion_factor = atof(argv[++i]);
      } else if (strncmp(argv[i], "-s", 2) == 0) {
         total_steps = atoi(argv[++i]);
      } else if (strncmp(argv[i], "-o", 2) == 0) {
         strncpy(output_format,argv[i]+2,4);
      } else
         (void) Usage(progname,0);
   }


   /* Read the input file */
   tri_head = read_input(infile,FALSE,NULL);


   /* Create the linked list of node connectivity (just uses node_head) */
   (void) set_node_connectivity();


   // first, fill basins to prevent flow (rivers) from disappearing
   fprintf(stderr,"Filling basins");
   for (step=0; step<100; step++) {
      fprintf(stderr,".");
      (void) find_flow(tri_head);
      (void) fill_basins(tri_head);
   }
   fprintf(stderr,"\n");


   // Then, run one or several erosion steps
   if (do_erosion)
   for (step=0; step<total_steps; step++) {

      fprintf(stderr,"Eroding %d.\n",step);

      // Calculate the flow
      (void) find_flow(tri_head);

      // Fill the basins (so flow goes somewhere)
      (void) fill_basins(tri_head);

      // Apply erosion
      (void) find_flow(tri_head);
      (void) erode_surface(tri_head,erosion_factor);
   }


   /* Write triangles to stdout */
   (void) write_output(tri_head,output_format,TRUE,argc,argv);


   /* Somehow, write flow data somewhere */
   /* (void) write_flow_data(); */

   fprintf(stderr,"Done.\n");
   exit(0);
}


/*
 * This function writes basic usage information to stderr,
 * and then quits. Too bad.
 */
int Usage(char progname[80],int status) {

   /* Usage for rockerode */
   static char **cpp, *help_message[] =
   {
       "where [-options] are one or more of the following:                         ",
       "                                                                           ",
       "   -e [val]    do an erosion calculation, optionally change the erosion    ",
       "               rate, which defaults to 1.0                                 ",
       "                                                                           ",
       "   -s num      compute num erosion steps                                   ",
       "                                                                           ",
       "   -okey       specify output format, key= raw, rad, pov, obj, tin, rib    ",
       "               default = raw                                               ",
       "                                                                           ",
       "   -help       (in place of infile) returns this help information          ",
       " ",
       "The input file can be of .obj, .raw, or .tin format, and the program       ",
       "   requires a valid 3-character extension for the input file.",
       " ",
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
