/*************************************************************
 *
 *  rocksmooth.c - Calculate surface normals of an arbitrary
 *	triangle mesh for smoothing purposes
 *
 *  Mark J. Stock, mstock@umich.edu
 *
 *
 * rocktools - Tools for creating and manipulating triangular meshes
 * Copyright (C) 1999,2004-06,14  Mark J. Stock
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

/* some essential global variables */
node_ptr node_head = NULL;
norm_ptr norm_head = NULL;

int num_tri = 0;

int Usage(char[80],int);
extern int three_d_laplace(tri_pointer,int);
extern int three_d_surface_tension(tri_pointer,double);
extern int compute_normals_2(tri_pointer,int);
extern int compute_normals_3(tri_pointer,int,int,double);
extern int grow_surface_along_normal(tri_pointer,double);


int main(int argc,char **argv) {

   int i;
   int do_laplace = FALSE;		/* perturb nodes to smooth shape? */
   int laplace_factor = 1;		/* amount of smoothing to take place */

   int do_tension = FALSE;		/* perturb nodes to smooth shape? */
   double tension_factor = 1.0;		/* amount of smoothing to take place */

   int do_normals = FALSE;		/* compute and write out normal vectors? */
   //double v_thresh = 1.0;		/* threshhold for common normal, convex edge */
   //double c_thresh = 1.0;		/* threshhold for common normal, concave edge */

   int do_invert = FALSE;		/* invert normal upon writing */

   int allow_sharp_edges = FALSE;	/* identify sharp edges with unique nodes
 					 * and normals */
   double edge_thresh = 45.;		/* threshhold in degrees for sharp edges */

   int grow_boundary = FALSE;		/* grow all boundaries? */
   double grow_distance = 0.01;		/* grow them by this much (can be negative) */

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
      if (strncmp(argv[i], "-s", 2) == 0) {
         do_laplace = TRUE;
         if (i < argc-1)
            if (strncmp(argv[i+1], "-", 1) != 0)
               laplace_factor = atoi(argv[++i]);
      } else if (strncmp(argv[i], "-a", 2) == 0) {
         allow_sharp_edges = TRUE;
         if (i < argc-1)
            if (strncmp(argv[i+1], "-", 1) != 0)
               edge_thresh = atof(argv[++i]);
         do_normals = TRUE;
      } else if (strncmp(argv[i], "-t", 2) == 0) {
         do_tension = TRUE;
         if (i < argc-1)
            if (strncmp(argv[i+1], "-", 1) != 0)
               tension_factor = atof(argv[++i]);
      } else if (strncmp(argv[i], "-n", 2) == 0) {
         do_normals = TRUE;
      } else if (strncmp(argv[i], "-i", 2) == 0) {
         do_invert = TRUE;
      } else if (strncmp(argv[i], "-grow", 2) == 0) {
         grow_boundary = TRUE;
         do_normals = TRUE;
         grow_distance = atof(argv[++i]);
      } else if (strncmp(argv[i], "-o", 2) == 0) {
         strncpy(output_format,argv[i]+2,4);
      } else
         (void) Usage(progname,0);
   }

   /* Run a check on the output format, not all support normal vectors */
   if (do_normals && strncmp(output_format, "rad", 3) == 0) {
      fprintf(stderr,"Cannot write smoothed triangles to .rad format, will write to .obj\n");
      strncpy(output_format,"obj\0",4);
   //} else if (do_normals && strncmp(output_format, "raw", 3) == 0) {
   //   fprintf(stderr,"Cannot write smoothed triangles to .raw format, will write to .tin\n");
   //   strncpy(output_format,"tin\0",4);
   }

   // Read the input file
   tri_head = read_input(infile,do_invert,NULL);

   // first, set node connectivity
   if (do_laplace) (void) set_node_connectivity();

   // Optionally smooth the surface by moving node locations - no normals needed
   if (do_laplace) (void) three_d_laplace(tri_head,laplace_factor);
   if (do_tension) (void) three_d_surface_tension(tri_head,tension_factor);

   // Define sharp edges by splitting nodes along the edges, thus any
   //   normal-finding routine that uses all connected tris will not see
   //   the offending tris, and the normals will be smoother
   //if (allow_sharp_edges) (void) define_sharp_edges(tri_head,edge_thresh);

   // grow surface
   if (grow_boundary) {
      // first, find *continuous* normals on the nodes;
      // but enforce that all normals at same node are same!
      // this prevents creating holes during the growth procedure
      (void) compute_normals_2(tri_head,3);
      // then, grow the surface along the normals
      (void) grow_surface_along_normal(tri_head,grow_distance);
   }

   // finally, find the normals on a per-element basis (not per-node!)
   // paying attention to the sharp edges, if necessary
   if (do_normals) (void) compute_normals_3(tri_head,3,allow_sharp_edges,edge_thresh);

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

  /* Usage for rocksmooth */
  static char **cpp, *help_message[] =
  {
    "where [-options] are one or more of the following:                         ",
    "                                                                           ",
    "   -s [val]    smooth the triangle mesh by flattening the more protruding  ",
    "               nodes, always done before calculating normals (-n)          ",
    "               Optional value indicates number of smoothing passes,        ",
    "               integer values only, default = 1                            ",
    "                                                                           ",
    "   -t [val]    smooth surface using surface tension algorithm, optional    ",
    "               argument is coefficient of surface tension, default = 1.0   ",
    "                                                                           ",
    "   -a [val]    allow edges between triangles with normal vectors varying   ",
    "               by more than val degrees to be considered sharp, and the    ",
    "               normal vector will be made discontinuous across that edge,  ",
    "               default = 45; use of this option implies -n                 ",
    "                                                                           ",
    "   -grow val   grow the surface val units in the normal direction, negative",
    "               values are allowed; use of this option implies -n; no       ",
    "               default exists, a value must be given on command-line       ",
    "                                                                           ",
    "   -i          flip surface normals                                        ",
    "                                                                           ",
    "   -n          calculate and output surface normals for all nodes,         ",
    "               output formats supported are: tin, pov, obj, rib            ",
    "                                                                           ",
    "   -okey       specify output format, key= raw, rad, pov, obj, tin, rib    ",
    "               default = raw, example: -oobj                               ",
    "                                                                           ",
    "   -help       (in place of infile) returns this help information          ",
    " ",
    "The input file can be of .raw, .tin, .obj format, and the program requires ",
    "   a valid lowercase 3-character filename extension.",
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
