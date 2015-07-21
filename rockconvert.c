/*************************************************************
 *
 *  rockconvert.c - Perform a file format conversion
 *	of an irregular triangle mesh
 *
 *  Mark J. Stock, mstock@umich.edu
 *
 *
 * rocktools - Tools for creating and manipulating triangular meshes
 * Copyright (C) 1999,2003-2004,2006-2007,14  Mark J. Stock
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

void rescale_nodes (VEC);
void translate_nodes (VEC);
int Usage(char[255],int);


int main(int argc,char **argv) {

   int i;
   char infile[255];		/* name of input file */
   char progname[255];		/* name of binary executable */
   char output_format[4];	/* file format extension for output */
   tri_pointer tri_head = NULL;
   int keep_normals = TRUE;
   int do_rescale = FALSE;
   int do_trans = FALSE;
   int do_trans_first = FALSE;
   VEC scale = {1., 1., 1.};
   VEC trans = {0., 0., 0.};


   /* Parse command-line args */
   (void) strcpy(progname,argv[0]);
   if (argc < 2) (void) Usage(progname,0);
   if (strncmp(argv[1], "-help", 2) == 0)
      (void) Usage(progname,0);
   (void) strcpy(infile,argv[1]);
   for (i=2; i<argc; i++) {
      if (strncmp(argv[i], "-o", 2) == 0) {
         strncpy(output_format,argv[i]+2,4);
      } else if (strncmp(argv[i], "-i", 2) == 0) {
         keep_normals = FALSE;
      } else if (strncmp(argv[i], "-s", 2) == 0) {
         do_rescale = TRUE;
         if (argv[i+1][0] != '-') {
            scale.x = atof(argv[++i]);
         }
         if (argv[i+1][0] != '-') {
            scale.y = atof(argv[++i]);
         } else {
            scale.y = scale.x;
         }
         if (argv[i+1][0] != '-') {
            scale.z = atof(argv[++i]);
         } else {
            scale.z = scale.y;
         }
      } else if (strncmp(argv[i], "-t", 2) == 0) {
         do_trans = TRUE;
         if (!do_rescale) do_trans_first = TRUE;
         if (argc > i) {
            if (!isalpha(argv[i+1][1])) {
               trans.x = atof(argv[++i]);
               trans.y = trans.x;
               trans.z = trans.y;
            }
         }
         if (argc > i) {
            if (!isalpha(argv[i+1][1])) {
               trans.y = atof(argv[++i]);
               trans.z = trans.y;
            }
         }
         if (argc > i) {
            if (!isalpha(argv[i+1][1])) {
               trans.z = atof(argv[++i]);
            }
         }
      } else {
         (void) Usage(progname,0);
      }
   }


   // Read the input file
   tri_head = read_input (infile,FALSE,NULL);

   // apply any transformations
   if (do_trans_first) translate_nodes (trans);
   if (do_rescale) rescale_nodes (scale);
   if (do_trans && !do_trans_first) translate_nodes (trans);

   // Write triangles to stdout
   (void) write_output (tri_head,output_format,keep_normals,argc,argv);

   fprintf(stderr,"Done.\n");
   exit(0);
}


/*
 * Scale all nodes
 */
void rescale_nodes (VEC scale) {

   node_ptr this_node;

   // then, perturb each according to the normal
   this_node = node_head;
   while (this_node) {
      this_node->loc.x *= scale.x;
      this_node->loc.y *= scale.y;
      this_node->loc.z *= scale.z;
      this_node = this_node->next_node;
   }

   return;
}


/*
 * Translate all nodes
 */
void translate_nodes (VEC trans) {

   node_ptr this_node;

   // then, perturb each according to the normal
   this_node = node_head;
   while (this_node) {
      this_node->loc.x += trans.x;
      this_node->loc.y += trans.y;
      this_node->loc.z += trans.z;
      this_node = this_node->next_node;
   }

   return;
}


/*
 * This function writes basic usage information to stderr,
 * and then quits. Too bad.
 */
int Usage(char progname[255],int status) {

   /* Usage for rockconvert */
   static char **cpp, *help_message[] =
   {
       "where [-options] are one or more of the following:                         ",
       "                                                                           ",
       "   -okey       specify output format, key= raw,rad,pov,obj,tin,rib,seg     ",
       "               default = raw; surface normal vectors are not supported     ",
       "                                                                           ",
       "   -s [x [y [z]]]                                                          ",
       "               scale geometry by the given x,y,z factors, default=1,1,1    ",
       "                                                                           ",
       "   -t [x [y [z]]]                                                          ",
       "               translate geometry by the given x,y,z, default=0,0,0        ",
       "                                                                           ",
       "   -ignore     ignore normals in output (strips normals)                   ",
       "                                                                           ",
       "   -help       (in place of infile) returns this help information          ",
       " ",
       "The input file can be of .obj, .raw, .msh, or .tin format, and the program",
       "   requires the input file to use its valid 3-character filename extension ",
       " ",
       "Options may be abbreviated to an unambiguous length",
       "Output is to stdout, so redirect it to a file using '>'",
       " ",
       "rockconvert converts tri meshes from one file format to another            ",
       "                                                                           ",
       "examples:                                                                  ",
       "   rockconvert in.tin -oobj > out.obj                                      ",
       "      Read in the file in.tin and write a Wavefront .obj format version    ",
       "                                                                           ",
       "   rockconvert in.raw -s 2 -orad > out.rad                                 ",
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
