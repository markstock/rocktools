/*************************************************************
 *
 *  rockbalance.c - Compute the convex hull of the input and
 *                  determine on which faces it will balance
 *
 *  Mark J. Stock, mstock@umich.edu
 *
 *
 * rocktools - Tools for creating and manipulating triangular meshes
 * Copyright (C) 1999,2006-7,9,14-15  Mark J. Stock
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

int Usage(char[MAX_FN_LEN],int);
extern tri_pointer create_convex_hull ();

//extern int write_output (tri_pointer, char[4], int, char**);

int main(int argc,char **argv) {

   int writemoststable = TRUE;		// write out most- or least-stable balanced conformation?
   char progname[MAX_FN_LEN];			// name of binary executable
   char infile[MAX_FN_LEN];			// name of the file to act upon
   char output_format[4];		// file format extension for output
   tri_pointer tri_head = NULL;

   /* Parse command-line args */
   (void) strcpy(progname,argv[0]);

   if (strncmp(argv[1], "-help", 2) == 0)
      (void) Usage(progname,0);
   (void) strcpy(infile,argv[1]);

   if (argc < 2) {
      (void) Usage(progname,0);
   } else {
      for (int i=2; i<argc; i++) {
         if (strncmp(argv[i], "-help", 2) == 0) {
            (void) Usage(progname,0);
         } else if (strncmp(argv[i], "-most", 2) == 0) {
            writemoststable = TRUE;
         } else if (strncmp(argv[i], "-least", 2) == 0) {
            writemoststable = FALSE;
         } else if (strncmp(argv[i], "-o", 2) == 0) {
            strncpy(output_format,argv[i]+2,3);
         } else if (strncmp(argv[i], "-", 1) == 0) {
            (void) Usage(progname,0);
         }
      }
   }


   // Read in the geometry from the command-line
   tri_head = read_input (infile,FALSE,NULL);

   // Compute the center of mass
   VEC cm = find_cm (tri_head);
   fprintf(stderr,"Center of mass is at %g %g %g\n", cm.x, cm.y, cm.z); fflush(stderr);

   // Translate object such that cm is origin
   node_ptr this_node = node_head;
   int num_nodes = 0;
   while (this_node) {
      this_node->loc.x -= cm.x;
      this_node->loc.y -= cm.y;
      this_node->loc.z -= cm.z;
      num_nodes++;
      this_node = this_node->next_node;
   }

   // Save the original geometry
   tri_pointer original = tri_head;

   // remove all original triangles
   //while (tri_head) {
   //   tri_head = delete_tri(tri_head);
   //}
   tri_head = NULL;

   // perform QuickHull algorithm to generate new trimesh
   fprintf(stderr,"Creating convex hull");
   if (num_nodes > 100000) fprintf(stderr," (this may take a while)");
   fprintf(stderr,"...\n");
   fflush(stderr);
   tri_head = create_convex_hull ();

   // Check each face to see if it will balance the object
   // and identify the most stable one
   tri_pointer this = tri_head;
   double maxStability = -1.0;
   tri_pointer stableTri = NULL;
   double minStability = 9.9e+9;
   tri_pointer leastStableTri = NULL;
   while (this) {
      // how far is CM from plane of this triangle?
      VEC thisnorm = find_tri_normal(this);
      //fprintf(stderr,"Elem with norm %g %g %g\n", thisnorm.x, thisnorm.y, thisnorm.z);
      //fprintf(stderr,"Elem\n");

      double height = dot(this->node[0]->loc, thisnorm);
      //fprintf(stderr,"  cm height is %g\n", height);

      // is CM within prism of this triangle?
      VEC localCM = vscale(height, thisnorm);
      double edge1 = dot( thisnorm, cross( norm(from(this->node[0]->loc, this->node[1]->loc)),
                                           from(this->node[0]->loc, localCM)));
      double edge2 = dot( thisnorm, cross( norm(from(this->node[1]->loc, this->node[2]->loc)),
                                           from(this->node[1]->loc, localCM)));
      double edge3 = dot( thisnorm, cross( norm(from(this->node[2]->loc, this->node[0]->loc)),
                                           from(this->node[2]->loc, localCM)));
      //fprintf(stderr,"  edges are %g %g %g\n", edge1, edge2, edge3);
      if (edge1 > 0.0 && edge2 > 0.0 && edge3 > 0.0) {
         //fprintf(stderr,"Elem with norm %g %g %g\n", thisnorm.x, thisnorm.y, thisnorm.z);
         //fprintf(stderr,"  cm height is %g\n", height);
         //fprintf(stderr,"  WILL BALANCE\n");
         double minedge = fmin(edge1, fmin(edge2, edge3));
         //fprintf(stderr,"  stability distance %g\n", minedge);
         double stability = minedge / height;
         if (stability > maxStability) {
            maxStability = stability;
            stableTri = this;
            //fprintf(stderr,"  MOST STABLE %g\n", maxStability);
         }
         if (stability < minStability) {
            minStability = stability;
            leastStableTri = this;
            //fprintf(stderr,"  LEAST STABLE %g\n", minStability);
         }
      }

      this = this->next_tri;
   }

   if (leastStableTri == NULL || stableTri == NULL) {
      fprintf(stderr,"No faces of convex hull will balance the object! Quitting without writing.\n");
      fflush(stderr);
      exit(1);
   }

   // Rotate first to one balanceable conformation (+z is up)
   VEC newx, newy, newz;
   if (writemoststable) {
      newx = norm(from(stableTri->node[0]->loc, stableTri->node[1]->loc));
      newz = vscale(-1.0, find_tri_normal(stableTri));
      newy = cross(newz, newx);
   } else {
      newx = norm(from(leastStableTri->node[0]->loc, leastStableTri->node[1]->loc));
      newz = vscale(-1.0, find_tri_normal(leastStableTri));
      newy = cross(newz, newx);
   }
   //fprintf(stderr,"New basis x is %g %g %g\n", newx.x, newx.y, newx.z);
   //fprintf(stderr,"New basis y is %g %g %g\n", newy.x, newy.y, newy.z);
   //fprintf(stderr,"New basis z is %g %g %g\n", newz.x, newz.y, newz.z);

   // Translate such that z=0 is floor
   VEC trans;
   trans.x = 0.0;
   trans.y = 0.0;
   if (writemoststable) {
      trans.z = -dot(stableTri->node[0]->loc, newz);
   } else {
      trans.z = -dot(leastStableTri->node[0]->loc, newz);
   }
   //fprintf(stderr,"Translation is %g %g %g\n", trans.x, trans.y, trans.z);

   this_node = node_head;
   while (this_node) {
      VEC newpos;
      newpos.x = dot(this_node->loc, newx);
      newpos.y = dot(this_node->loc, newy);
      newpos.z = dot(this_node->loc, newz);
      this_node->loc.x = newpos.x + trans.x;
      this_node->loc.y = newpos.y + trans.y;
      this_node->loc.z = newpos.z + trans.z;
      this_node = this_node->next_node;
   }

   // Write one balanceable conformation to the output
   //(void) write_output(tri_head,output_format,TRUE,argc,argv);
   (void) write_output(original,output_format,TRUE,argc,argv);

   // ...aaaaand done.
   fprintf(stderr,"Done.\n");
   exit(0);
}


/*
 * This function writes basic usage information to stderr,
 * and then quits. Too bad.
 */
int Usage(char progname[MAX_FN_LEN],int status) {

   /* Usage for rockbalance */
   static char **cpp, *help_message[] =
   {
       "where [-options] are one or more of the following:                         ",
       "                                                                           ",
       "   -most       write out the most stable conformation, z>0, cm_x=cm_y=0    ",
       "                                                                           ",
       "   -least      write out the least stable conformation, z>0, cm_x=cm_y=0   ",
       "                                                                           ",
       "   -okey       specify output format, key= raw, rad, pov, obj, tin, rib    ",
       "               default = raw; surface normal vectors are not written       ",
       "                                                                           ",
       "   -help       print usage information                                     ",
       " ",
       "The input file can be of .obj, .raw, .msh, or .tin format, and the program",
       "   requires the input file to use its valid 3-character filename extension ",
       " ",
       "Options may be abbreviated to an unambiguous length",
       "Output is to stdout, so redirect it to a file using '>'",
       " ",
       "rockcreate creates a closed triangle mesh of a rock-like object by         ",
       "   performing a convex hull calculation over a set of input points         ",
       "                                                                           ",
       NULL
   };

   fprintf(stderr, "usage:\n  %s infile [-options]\n\n", progname);
   for (cpp = help_message; *cpp; cpp++) fprintf(stderr, "%s\n", *cpp);
   fflush(stderr);
   exit(status);
   return(0);
}
