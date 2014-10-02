/*************************************************************
 *
 *  rockcreate.c - Create a general rock shape from a convex
 *	hull calculation
 *
 *  Mark J. Stock, mstock@umich.edu
 *
 *
 * rocktools - Tools for creating and manipulating triangular meshes
 * Copyright (C) 1999,2006-7,9  Mark J. Stock
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

int Usage(char[80],int);
extern int read_files_for_nodes (int,char**);
extern int create_cubic_nodes (int);
extern int create_gaussian_nodes (int);
extern int create_random_walk_nodes (int);
extern int sphericalize_nodes (int);
extern int pack_into_cube ();
extern tri_pointer create_convex_hull ();

//extern int convert_simplex_to_tris(simplex*);
//extern int write_output (tri_pointer, char[4], int, char**);


int main(int argc,char **argv) {

   int i;
   int rand_seed = 1;			// seed for the random number generator
   int num_input = 0;			// number of nodes from files
   int num_cube = 0;			// number of random cube nodes
   int num_gauss = 0;			// number of Gaussian random nodes
   int num_walk = 0;			// number of random walk nodes
   int tot_nodes;
   double roundness = 0.0;		// relative roundness, 1.0 = average
   char progname[80];			// name of binary executable
   char output_format[4];		// file format extension for output
   int num_its;
   tri_pointer tri_head = NULL;

   /* Parse command-line args */
   (void) strcpy(progname,argv[0]);
   if (argc < 2) {
      (void) Usage(progname,0);
   } else {
      if (strncmp(argv[1], "-help", 2) == 0)
         (void) Usage(progname,0);
      for (i=1; i<argc; i++) {
         if (strncmp(argv[i], "-s", 2) == 0) {
            rand_seed = atoi(argv[++i]);
            srand((unsigned int) rand_seed);
         } else if (strncmp(argv[i], "-n", 2) == 0) {
            if (isdigit(argv[i+1][0])) {
               num_cube = atoi(argv[++i]);
            } else {
               num_cube = 6+(int)(18.0*((rand()+1.0)/RAND_MAX));
            }
         } else if (strncmp(argv[i], "-g", 2) == 0) {
            if (isdigit(argv[i+1][0])) {
               num_gauss = atoi(argv[++i]);
            } else {
               num_gauss = 6+(int)(18.0*((rand()+1.0)/RAND_MAX));
            }
         } else if (strncmp(argv[i], "-w", 2) == 0) {
            if (isdigit(argv[i+1][0])) {
               num_walk = atoi(argv[++i]);
            } else {
               num_walk = 6+(int)(18.0*((rand()+1.0)/RAND_MAX));
            }
         } else if (strncmp(argv[i], "-r", 2) == 0) {
            roundness = atof(argv[++i]);
         } else if (strncmp(argv[i], "-o", 2) == 0) {
            strncpy(output_format,argv[i]+2,4);
         } else if (strncmp(argv[i], "-", 1) == 0) {
            (void) Usage(progname,0);
         }
      }
   }


   // Read in all nodes from files on the command-line
   num_input = read_files_for_nodes (argc,argv);
   if (num_input > 0) {
      fprintf(stderr,"After reading inputs, have %d nodes\n",num_input);
      fflush(stderr);
   }

   // Is this enough, or too many nodes?
   tot_nodes = num_input+num_cube+num_gauss+num_walk;
   if (tot_nodes > 100000) {
      fprintf(stderr,"This might take a while---rockcreate prefers fewer than 100000 nodes.\n");
   } else if (tot_nodes < 4) {
      fprintf(stderr,"Rockcreate cannot support under 4 nodes, increasing count to 4.\n");
      num_cube += 4-tot_nodes;
   }


   // Create the procedural nodes
   if (num_cube > 0) (void) create_cubic_nodes (num_cube);
   if (num_gauss > 0) (void) create_gaussian_nodes (num_gauss);
   if (num_walk > 0) (void) create_random_walk_nodes (num_walk);


   // move the nodes into a more-evenly-distributed arrangement
   //num_its = (int)(roundness*1000.0/sqrt((float)num_sites));
   num_its = (int)(roundness*50.0);
   if (num_its > 0) {
      fprintf(stderr,"Using %d point skewing iterations on %d sites...\n",num_its,tot_nodes);
      fflush(stderr);
      (void) sphericalize_nodes (num_its);
   }


   // call QuickHull routine to generate trimesh
   fprintf(stderr,"Creating a rock from %d initial sites...\n",tot_nodes); fflush(stderr);
   tri_head = create_convex_hull ();


   // pack into a unit cube
   //(void) pack_into_cube ();


   // Write triangles to stdout
   (void) write_output(tri_head,output_format,TRUE,argc,argv);


   fprintf(stderr,"Done.\n");
   exit(0);
}


/*
 * This function writes basic usage information to stderr,
 * and then quits. Too bad.
 */
int Usage(char progname[80],int status) {

   /* Usage for rockcreate */
   static char **cpp, *help_message[] =
   {
       "where [-options] are one or more of the following:                         ",
       "                                                                           ",
       "   filename    read nodes from an input file, must be either .obj, .tin,   ",
       "               or .raw, and must use the proper extension                  ",
       "                                                                           ",
       "   -n [val]    create val random nodes within a unit cube, if no number    ",
       "               is given, program will choose between 6 and 24              ",
       "                                                                           ",
       "   -g [val]    create val random nodes using Gaussian distribution, if no  ",
       "               number is given, program will choose between 6 and 24       ",
       "                                                                           ",
       "   -w [val]    create val nodes, each one step in a random walk starting   ",
       "               at the origin, if no number is given, the program will      ",
       "               choose between 6 and 24                                     ",
       "                                                                           ",
       "               filename, -u, -g, and -w options are not exclusive; nodes   ",
       "               may be added using any or all of those methods              ",
       "                                                                           ",
       "               if no node-creation options are listed, -u will be used     ",
       "                                                                           ",
       "   -r val      change the relative roundness of the resulting shape,       ",
       "               val <1.0 results in more small faces and greater variety,   ",
       "               while >1.0 makes the shape closer to a sphere, default=0.0  ",
       "                                                                           ",
       "   -s val      seed the random number generator with an unsigned integer,  ",
       "               default=1                                                   ",
       "                                                                           ",
       "   -okey       specify output format, key= raw, rad, pov, obj, tin, rib    ",
       "               default = raw; surface normal vectors are not written       ",
       "                                                                           ",
       "   -help       print usage information                                     ",
       " ",
       "Options may be abbreviated to an unambiguous length",
       "Output is to stdout, so redirect it to a file using '>'",
       " ",
       "rockcreate creates a closed triangle mesh of a rock-like object by         ",
       "   performing a convex hull calculation over a set of input points         ",
       "                                                                           ",
       "examples:                                                                  ",
       "   rockcreate -oobj > out.obj                                              ",
       "      Just generate a random rock in Wavefront .obj format using 6-24      ",
       "      random points in a cube                                              ",
       "                                                                           ",
       "   rockcreate -g 200 -r 3 -s 1234 -otin > out.tin                          ",
       "      Generate a .tin-format rock file from a convex hull of 200 points    ",
       "      in a random Gaussian distribution, but smooth or round out the       ",
       "      points before performing the convex hull algorithm                   ",
       "                                                                           ",
       "   rockcreate mydots.obj -w 100 more.raw -oraw > out.raw                   ",
       "      Read data points from a .obj file called mydots.obj (will read any   ",
       "      line beginning with 'v ') and from a raw text file called more.raw   ",
       "      (will read any line beginning with three floats) and then add 100    ",
       "      more points along a random walk, then wrap a convex hull around      ",
       "      the entire collection and write the resulting tri mesh in raw format ",
       "                                                                           ",
       NULL
   };

   fprintf(stderr, "usage:\n  %s [-options]\n\n", progname);
   for (cpp = help_message; *cpp; cpp++)
      fprintf(stderr, "%s\n", *cpp);
      fflush(stderr);
   exit(status);
   return(0);
}
