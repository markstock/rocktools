/*************************************************************
 *
 *  rockinfo.c - Dump out some key info for each input file
 *
 *  Mark J. Stock, mstock@umich.edu
 *
 *
 * rocktools - Tools for creating and manipulating triangular meshes
 * Copyright (C) 1999,2004,2006,2008 Mark J. Stock
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
#include <ctype.h>

#define INFO
#include "structs.h"


// This needs to be here to use routines in utils.c
node_ptr node_head = NULL;

// this subroutine appears in this file
int Usage(char[80],int);

// from inout.c
extern int find_mesh_stats(char*, VEC*, VEC*, int*, int*);
extern int get_tri(FILE*, int, tri_pointer, char[512]);


int main(int argc,char **argv) {

   int num_nodes = 0;
   int num_tris = 0;
   VEC bmax,bmin;		// mesh bounds
   char infile[80];		/* name of input file */
   char progname[80];		/* name of binary executable */

   bmin.x = 9.9e+9;
   bmax.x = -9.9e+9;
   bmin.y = 9.9e+9;
   bmax.y = -9.9e+9;
   bmin.z = 9.9e+9;
   bmax.z = -9.9e+9;

   /* Parse command-line args */
   (void) strcpy(progname,argv[0]);
   if (argc > 2) (void) Usage(progname,0);
   if (strncmp(argv[1], "-help", 2) == 0)
      (void) Usage(progname,0);
   (void) strcpy(infile,argv[1]);


   // external subroutine does all the work
   find_mesh_stats(infile,&bmin,&bmax,&num_tris,&num_nodes);

   fprintf(stdout,"nodes %d tris %d ",num_nodes,num_tris);
   fprintf(stdout,"x %g %g y %g %g z %g %g\n",bmin.x,bmax.x,bmin.y,bmax.y,bmin.z,bmax.z);

   exit(0);
}


/*
 * This function writes basic usage information to stderr,
 * and then quits. Too bad.
 */
int Usage(char progname[80],int status) {

   /* Usage for rockinfo */
   static char **cpp, *help_message[] =
   {
       " ",
       "The input file can be raw, tin, rad, or obj format, and the program requires",
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
