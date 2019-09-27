/*************************************************************
 *
 *  rockinfo.c - Dump out some key info for each input file
 *
 *  Mark J. Stock, mstock@umich.edu
 *
 *
 * rocktools - Tools for creating and manipulating triangular meshes
 * Copyright (C) 1999,2004,2006,2008,2015 Mark J. Stock
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

node_ptr node_head = NULL;
norm_ptr norm_head = NULL;
text_ptr text_head = NULL;

// this subroutine appears in this file
int Usage(char[255],int);

// from inout.c
//extern int find_mesh_stats(char*, VEC*, VEC*, int, VEC*, int*, int*);


int main(int argc,char **argv) {

   int doCM = FALSE;
   int num_nodes = 0;
   int num_tris = 0;
   VEC bmax,bmin;		// mesh bounds
   VEC cm;			// center of mass
   float volume;		// enclosed volume
   char infile[255];		// name of input file
   char progname[255];	// name of binary executable

   /* Parse command-line args */
   (void) strcpy(progname,argv[0]);
   if (argc < 2) (void) Usage(progname,0);
   if (strncmp(argv[1], "-help", 2) == 0)
      (void) Usage(progname,0);
   for (int i=2; i<argc; i++) {
      if (strncmp(argv[i], "-cm", 2) == 0) {
         doCM = TRUE;
      } else {
         (void) Usage(progname,0);
      }
   }
   (void) strcpy(infile,argv[1]);

   // external subroutine does all the work
   find_mesh_stats(infile,&bmin,&bmax,doCM,&cm,&volume,&num_tris,&num_nodes);

   fprintf(stdout,"nodes %d tris %d ",num_nodes,num_tris);
   fprintf(stdout,"x %g %g y %g %g z %g %g",bmin.x,bmax.x,bmin.y,bmax.y,bmin.z,bmax.z);
   if (doCM) fprintf(stdout," cm %g %g %g vol %g",cm.x,cm.y,cm.z,volume);
   fprintf(stdout,"\n");

   exit(0);
}


/*
 * This function writes basic usage information to stderr,
 * and then quits. Too bad.
 */
int Usage(char progname[255],int status) {

   /* Usage for rockinfo */
   static char **cpp, *help_message[] =
   {
       "   -cm         also compute and print center of mass                       ",
       "                                                                           ",
       "   -help       (in place of infile) returns this help information          ",
       " ",
       "The input file can be raw, tin, rad, or obj format, and the program requires",
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
