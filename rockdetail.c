/*************************************************************
 *
 *  rockdetail.c - Recursively detail a rock surface composed
 *	of an irregular triangle mesh
 *
 *  Mark J. Stock, mstock@umich.edu
 *
 *
 * rocktools - Tools for creating and manipulating triangular meshes
 * Copyright (C) 1999,2003,4,6  Mark J. Stock
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
//#include <malloc.h>
#include "structs.h"


node_ptr node_head = NULL;

int num_tri = 0;
double normal_shake = 0.1;	/* the non-dimensionalized maximum perturbation of a midpoint
                                 * in the direction normal to the average plane of the two
                                 * triangles sharing the edge */
double normal_exponent = 0.5;	/* exponent for new point perturbation, normal to surface plane */
double normal_bias = 0.0;	/* the amount which to add to the normal perturbation before scaling */

double base_shake = 0.1;	/* the absolute non-dimensionalized maximum perturbation of a midpoint
                                 * unless the _3 or _4 methods are used, then this is the
                                 * non-dimensionalized maximum perturbation in the average
                                 * plane of the two triangles sharing the edge */
double base_exponent = 0.5;	/* exponent for new point perturbation, surface plane */

int clamp_edges = FALSE;	/* clamp edges to maintain straight lines (TRUE) or 
                                 * allow them to become jagged as they are split */
int use_hex_splitting = FALSE;	/* use hexagonal-like split, each tri makes 3, not 4 */
int use_spline = TRUE;		// new nodes are placed using cubic splines
int perturb_older_nodes = TRUE;	// turning this off recreates the old behavior where
				//   older nodes are not additionally perturbed
int use_gaussian_random = FALSE;	// use Gaussian random numbers, range is then
				// standard deviation


int use_thresh = FALSE;		/* if TRUE, impose area threshhold limit on triangle splitting */
double area_thresh = 0.0001;	/* area below which to *not* split triangle */

int use_dist = FALSE;		/* if TRUE, impose distance/area threshhold on tri splitting */
double distance_thresh = 0.01;	/* ratio of width:distance above which to continue splitting */
VEC viewp = {0.0, 0.0, 0.0};	/* viewpoint to use as origin for distances */

int force_sphere = FALSE;	// force the all points to be the same radius from the center
 				// of the rock---all rocks become spheres
double sphere_rad = -1.;	// this is the radius; if negative, compute the radius


int Usage(char [80],int );
extern tri_pointer split_tri(int,tri_pointer);
extern tri_pointer split_tri_hex(int,tri_pointer);
extern tri_pointer split_tri_5(int,tri_pointer);


//#define _GNU_SOURCE
//#include <fenv.h>
//#include <signal.h>
//#include <stdlib.h>

//void fpehandler(int sig_num)
//{
//        signal(SIGFPE, fpehandler);
//        printf("SIGFPE: floating point exception occured, exiting.\n");
//        abort();
//}


int main(int argc,char **argv) {

   int i;
   int depth;					/* currently-processing recursion depth */
   int end_depth = 1;				/* number of recursion levels to calculate */
   int rand_seed = 1;				/* seed for the random number generator */
   //int num_wrote = 0;				/* number of triangles written out */
   double area;
   char infile[80];				/* name of input file */
   char progname[80];				/* name of binary executable */
   char output_format[4];			/* file format extension for output */
   tri_pointer tri_head = NULL;
   tri_pointer curr = NULL;

        //int feenableexcept();
        //feenableexcept(FE_ALL_EXCEPT);
        //signal(SIGFPE, fpehandler);

   /* Parse command-line args */
   (void) strcpy(progname,argv[0]);
   if (argc < 2) (void) Usage(progname,0);
   if (strncmp(argv[1], "-help", 2) == 0)
      (void) Usage(progname,0);
   (void) strcpy(infile,argv[1]);
   for (i=2; i<argc; i++) {
      if (strncmp(argv[i], "-dt", 3) == 0) {
         use_dist = TRUE;
         distance_thresh = atof(argv[++i]);
         viewp.x = atof(argv[++i]);
         viewp.y = atof(argv[++i]);
         viewp.z = atof(argv[++i]);
      } else if (strncmp(argv[i], "-d", 2) == 0) {
         end_depth = atoi(argv[++i]);
         if (end_depth < 0) {
            fprintf(stderr,"Recursion depth can not be negative, resetting to 0\n");
            end_depth = 0;
         } else if (end_depth > 10) {
            fprintf(stderr,"Recursion depth should not be more than 8, resetting to 10\n");
            end_depth = 10;
         }
      } else if (strncmp(argv[i], "-be", 3) == 0) {
         base_exponent = atof(argv[++i]);
      } else if (strncmp(argv[i], "-b", 2) == 0) {
         base_shake = atof(argv[++i]);
      } else if (strncmp(argv[i], "-ne", 3) == 0) {
         normal_exponent = atof(argv[++i]);
      } else if (strncmp(argv[i], "-nb", 3) == 0) {
         normal_bias = atof(argv[++i]);
      } else if (strncmp(argv[i], "-n", 2) == 0) {
         normal_shake = atof(argv[++i]);
      } else if (strncmp(argv[i], "-mid", 4) == 0) {
         use_spline = FALSE;
      } else if (strncmp(argv[i], "-spl", 4) == 0) {
         use_spline = TRUE;
      } else if (strncmp(argv[i], "-sph", 4) == 0) {
         force_sphere = TRUE;
         if (argv[i+1][0]!='-') sphere_rad = atof(argv[++i]);
      } else if (strncmp(argv[i], "-se", 3) == 0) {
         rand_seed = atoi(argv[++i]);
      } else if (strncmp(argv[i], "-o", 2) == 0) {
         strncpy(output_format,argv[i]+2,4);
      } else if (strncmp(argv[i], "-3", 2) == 0) {
         use_hex_splitting = TRUE;
      } else if (strncmp(argv[i], "-4", 2) == 0) {
         use_hex_splitting = FALSE;
      } else if (strncmp(argv[i], "-ce", 2) == 0) {
         clamp_edges = TRUE;
      } else if (strncmp(argv[i], "-gr", 2) == 0) {
         use_gaussian_random = TRUE;
      } else if (strncmp(argv[i], "-at", 3) == 0) {
         area_thresh = atof(argv[++i]);
         use_thresh = TRUE;
      } else
         (void) Usage(progname,0);
   }
   srand((unsigned int) rand_seed);


   // Read the input file
   tri_head = read_input(infile,FALSE,NULL);

   // we are using the splittable flag to denote splittability, TRUE=can split it
   // initially flag all triangles that are not allowed to split at all
   curr = tri_head;
   while (curr) {
      curr->splittable = TRUE;
      if (use_thresh) {
         area = find_area(curr);
         if (area < area_thresh) curr->splittable = FALSE;
      }
      if (use_dist) {
         area = find_area(curr);
         if (sqrt(area)/find_tri_dist(curr,viewp) < distance_thresh)
            curr->splittable = FALSE;
      }
      curr = curr->next_tri;
   }


   // Recursively detail the mesh
   for (depth=0; depth<end_depth; depth++) {

      // take the old triangle list, and split each one
      if (use_hex_splitting) {
         tri_head = split_tri_hex (depth,tri_head);
      } else {
         if (perturb_older_nodes) {
            // this method places new nodes, then perturbs ALL nodes
            tri_head = split_tri_5 (depth,tri_head);
         } else {
            // this is the old method, it doesn't move the older nodes
            tri_head = split_tri (depth,tri_head);
         }
      }
   }


   // we don't know normals anymore, so don't print them
   //if (!force_sphere) {
   //  curr = tri_head;
   //  while (curr) {
   //    curr->use_norm = FALSE;
   //    curr = curr->next_tri;
   //  }
   //}

   /* Write triangles to stdout */
   (void) write_output(tri_head,output_format,argc,argv);

   fprintf(stderr,"Done.\n");
   exit(0);
}


/*
 * This function writes basic usage information to stderr,
 * and then quits. Too bad.
 */
int Usage(char progname[80],int status) {

   /* Usage for rockdetail */
   static char **cpp, *help_message[] =
   {
       "where [-options] are one or more of the following:                         ",
       "   -d val      specify levels of recursion to perform, default=1           ",
       "               each level creates 4 times as many triangles as the last    ",
       "               note: the usable values are currently 0 to 10               ",
       "                                                                           ",
       "   -4          use standard 1-triangle-becomes-4-triangles scheme, default ",
       "                                                                           ",
       "   -3          use hexagonal-like subdivision, each tri makes 3, not 4     ",
       "               subtriangles; first proposed by B. Mandelbrot               ",
       "                                                                           ",
       "   -spl        new node placement uses cubic spline interpolation, default ",
       "                                                                           ",
       "   -mid        new node placement uses linear interpolation (midpoint)     ",
       "                                                                           ",
       "   -b val      set the base perturbation amount for the first recursion,   ",
       "               default = 0.1 this is the maximum amount the placed midpoint ",
       "               will be in the average surface plane, setting this above 0.5",
       "               will create unevenly-sized triangles                        ",
       "                                                                           ",
       "   -be val     set the exponent on the perturbation, <1 more uniform, >1 is",
       "               rougher this exponent only affects the planar perturbation, ",
       "               default = 0.5                                               ",
       "                                                                           ",
       "   -n val      set the normal perturbation amount for the first recursion, ",
       "               default = 0.1, this is the maximum amount the placed midpoint ",
       "               will be normal to the average surface plane                 ",
       "                                                                           ",
       "   -ne val     set the exponent on the perturbation, <1 is smoother, >1    ",
       "               makes spikes this exponent only affects the normal          ",
       "               perturbation, default = 0.5                                 ",
       "                                                                           ",
       "   -nb val     set the normal bias to val, this is a value added to the    ",
       "               normal perturbation before scaling, a positive value        ",
       "               makes the resulting shape more bubbly, default = 0.0,       ",
       "               this option is not necessary if spline interpolation is used",
       "                                                                           ",
       "   -ce         clamp the surface's open edges to maintain them as straight ",
       "               lines, dafault is to allow open edges to become jagged      ",
       "                                                                           ",
       "   -at val     use triangle area threshhold, val is minimum area of any    ",
       "               triangle to be split, default = 0.0001                      ",
       "                                                                           ",
       "   -dt val x y z   use distance threshhold to control triangle splitting,  ",
       "               val is ratio of width:distance above which a triangle       ",
       "               will be split, x, y, and z are the coordinates of the       ",
       "               point from which to measure the distance, default = 0.01    ",
       "                                                                           ",
       "   -sph        force the shape to become a sphere at every level           ",
       "                                                                           ",
       "   -seed vel   seed the random number generator with an unsigned integer,  ",
       "               defaul t=1                                                  ",
       "                                                                           ",
       "   -gr         use Gaussian random numbers; all perturbations become       ",
       "               standard deviations; normal bias is scaled to std dev.      ",
       "                                                                           ",
       "   -okey       specify output format, key= raw, rad, pov, obj, tin, rib    ",
       "               default = raw; surface normal vectors are not supported     ",
       "                                                                           ",
       "   -help       (in place of infile) returns this help information          ",
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
