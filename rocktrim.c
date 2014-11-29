/*************************************************************
 *
 *  rocktrim.c - Selectively cut triangles from a given
 *		irregular triangle mesh
 *
 *  Mark J. Stock, mstock@umich.edu
 *
 *
 * rocktools - Tools for creating and manipulating triangular meshes
 * Copyright (C) 1999,2004,2006-8,14  Mark J. Stock
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

#define TRIM
#include "structs.h"


/* This needs to be here to use routines in utils.c */
node_ptr node_head = NULL;
norm_ptr norm_head = NULL;

double find_area(tri_pointer);
int Usage(char[80],int);

int main(int argc,char **argv) {

   int i;
   int input_format = 0;	/* integer flag for input file type */
   int output_format = 0;	/* integer flag for input file type */
   int num_read = 0;
   int num_wrote = 0;		/* number of triangles written out */
   //int have_normals = 0;
   //int inclusive = TRUE;	/* use the inclusive scheme for trimming */
   int scheme = 3;		/* 1 means use the inclusive scheme for trimming */
   				/* 2 means use exclusive */
   				/* 3 means use smooth trimming */
   int keep_tri;		/* keep the current tri or not? */
   int trim_tri;		/* trim the current tri or not? */
   int node_is_good[3];
   double thresh = 0.;		/* the cutoff threshhold for whatever axis */
   double frac = 0.;		/* more dealing with the cutoff */
   double x_min = 0.0;		/* cut all tris completely below x_min */
   double x_max = 1.0;		/* cut all tris completely above x_max */
   double y_min = 0.0;		/* cut all tris completely below y_min */
   double y_max = 1.0;		/* cut all tris completely above y_max */
   double z_min = 0.0;		/* cut all tris completely below z_min */
   double z_max = 1.0;		/* cut all tris completely above z_max */
   double a_min = 0.0;		/* cut all tris with area below a_min */
   double tri_area;		/* area of the tested triangle */
   int use_x_min = 0;
   int use_x_max = 0;
   int use_y_min = 0;
   int use_y_max = 0;
   int use_z_min = 0;
   int use_z_max = 0;
   int use_a_min = 0;
   int vn1,vn2,vn3;
   char infile[80];		/* name of input file */
   char extension[4];		/* filename extension if infile */
   char output_string[4] = "raw"; /* format extension for the output */
   char progname[80];		/* name of binary executable */
   tri_pointer the_tri,ttri1,ttri2;
   node_ptr the_nodes[3],tnode1,tnode2;
   FILE *ifp;


   /* Parse command-line args */
   (void) strcpy(progname,argv[0]);
   if (argc < 2) (void) Usage(progname,0);
   if (strncmp(argv[1], "-help", 2) == 0)
      (void) Usage(progname,0);
   (void) strcpy(infile,argv[1]);
   for (i=2; i<argc; i++) {
      if (strncmp(argv[i], "-x", 2) == 0) {
         x_min = atof(argv[++i]);
         use_x_min = 1;
      } else if (strncmp(argv[i], "+x", 2) == 0) {
         x_max = atof(argv[++i]);
         use_x_max = 1;
      } else if (strncmp(argv[i], "-y", 2) == 0) {
         y_min = atof(argv[++i]);
         use_y_min = 1;
      } else if (strncmp(argv[i], "+y", 2) == 0) {
         y_max = atof(argv[++i]);
         use_y_max = 1;
      } else if (strncmp(argv[i], "-z", 2) == 0) {
         z_min = atof(argv[++i]);
         use_z_min = 1;
      } else if (strncmp(argv[i], "+z", 2) == 0) {
         z_max = atof(argv[++i]);
         use_z_max = 1;
      } else if (strncmp(argv[i], "-i", 2) == 0) {
         scheme = 1;
      } else if (strncmp(argv[i], "-e", 2) == 0) {
         scheme = 2;
      } else if (strncmp(argv[i], "-s", 2) == 0) {
         scheme = 3;
      } else if (strncmp(argv[i], "-o", 2) == 0) {
         strncpy(output_string,argv[i]+2,4);
      } else if (strncmp(argv[i], "-ma", 3) == 0) {
         a_min = atof(argv[++i]);
         use_a_min = 1;
      } else
         (void) Usage(progname,0);
   }

   /* Read the input file */
   // tri_head = read_input(infile,FALSE,NULL);


   /* Determine the input file format from the .XXX extension, and read it */
   strncpy(extension,infile+strlen(infile)-3,4);
   if (strncmp(extension, "raw", 3) == 0)
      input_format = 1;
   else if (strncmp(extension, "tin", 1) == 0)
      input_format = 2;
   else if (strncmp(extension, "rad", 1) == 0)
      input_format = 3;
   else {
      fprintf(stderr,"Input filename extension is assumed to be (%s)\n",extension);
      fprintf(stderr,"This is either not a supported input file format for rocktrim, or you\n");
      fprintf(stderr,"   need to use a proper extension (last 3 characters in filename).\n");
      fprintf(stderr,"Supported input file formats are: .raw, .tin, and .rad\n");
      exit(0);
   }


   /* Determine and set the output format key from the string */
   if (strncmp(output_string, "raw", 3) == 0)
      output_format = 1;
   else if (strncmp(output_string, "tin", 3) == 0)
      output_format = 2;
   else if (strncmp(output_string, "rad", 3) == 0)
      output_format = 3;
   else {
      fprintf(stderr,"Output filename extension is assumed to be (%s)\n",output_string);
      fprintf(stderr,"This is either not a supported output file format for rocktrim, or you\n");
      fprintf(stderr,"   need to use a proper extension (standard 3-character extension for format).\n");
      fprintf(stderr,"Supported output file formats are: .raw, .tin, and .rad\n");
      exit(0);
   }


   /* Set up memory space for the working triangle, and the three nodes */
   the_tri = (TRI *)malloc(sizeof(TRI));
   for (i=0; i<3; i++) {
      the_nodes[i] = (NODE *)malloc(sizeof(NODE));
      the_tri->node[i] = the_nodes[i];
      the_tri->norm[i] = NULL;
   }
   /* these are for the nodes and triangles in case a triangle needs a corner trimmed */
   tnode1 = (NODE *)malloc(sizeof(NODE));
   tnode2 = (NODE *)malloc(sizeof(NODE));
   ttri1 = (TRI *)malloc(sizeof(TRI));
   ttri2 = (TRI *)malloc(sizeof(TRI));
   for (i=0; i<3; i++) {
      ttri1->norm[i] = NULL;
      ttri2->norm[i] = NULL;
   }


   /* Read the input file */

   /* open the file for reading */
   ifp = fopen(infile,"r");
   if (ifp==NULL) {
      fprintf(stderr,"Could not open input file %s\n",infile);
      exit(0);
   }
   fprintf(stderr,"Opening file %s, trimming",infile);
   fflush(stderr);


   /* as long as there are triangles available, operate */
   while (get_tri(ifp,input_format,the_tri) == 1) {

      num_read++;
      keep_tri = FALSE;
      trim_tri = FALSE;

      if (scheme == 1) {

         /* use the inclusive scheme: */
         /* if any of the three nodes are inside of the boundary, keep it */
         /* so, the point of the 'if's is to prove that the tri should be kept */

         for (i=0; i<3; i++) node_is_good[i] = TRUE;

         if (use_x_min)
            for (i=0; i<3; i++)
               if (the_tri->node[i]->loc.x < x_min)
                  node_is_good[i] = FALSE;

         if (use_x_max)
            for (i=0; i<3; i++)
               if (the_tri->node[i]->loc.x > x_max)
                  node_is_good[i] = FALSE;

         if (use_y_min)
            for (i=0; i<3; i++)
               if (the_tri->node[i]->loc.y < y_min)
                  node_is_good[i] = FALSE;

         if (use_y_max)
            for (i=0; i<3; i++)
               if (the_tri->node[i]->loc.y > y_max)
                  node_is_good[i] = FALSE;

         if (use_z_min)
            for (i=0; i<3; i++)
               if (the_tri->node[i]->loc.z < z_min)
                  node_is_good[i] = FALSE;

         if (use_z_max)
            for (i=0; i<3; i++)
               if (the_tri->node[i]->loc.z > z_max)
                  node_is_good[i] = FALSE;

         if (node_is_good[0] || node_is_good[1] || node_is_good[2])
            keep_tri = TRUE;
         else
            keep_tri = FALSE;

      } else if (scheme == 2) {

         /* if any of the three nodes are outside of the boundary, trim it */
         keep_tri = TRUE;

         if (keep_tri && use_x_min)
            for (i=0; i<3; i++)
               if (the_tri->node[i]->loc.x < x_min)
                  keep_tri = FALSE;

         if (keep_tri && use_x_max)
            for (i=0; i<3; i++)
               if (the_tri->node[i]->loc.x > x_max)
                  keep_tri = FALSE;

         if (keep_tri && use_y_min)
            for (i=0; i<3; i++)
               if (the_tri->node[i]->loc.y < y_min)
                  keep_tri = FALSE;

         if (keep_tri && use_y_max)
            for (i=0; i<3; i++)
               if (the_tri->node[i]->loc.y > y_max)
                  keep_tri = FALSE;

         if (keep_tri && use_z_min)
            for (i=0; i<3; i++)
               if (the_tri->node[i]->loc.z < z_min)
                  keep_tri = FALSE;

         if (keep_tri && use_z_max)
            for (i=0; i<3; i++)
               if (the_tri->node[i]->loc.z > z_max)
                  keep_tri = FALSE;

      } else if (scheme == 3) {

         /* use the smooth trimming scheme */
         /* begin by assuming we lose all triangles. Any triangles
          * that have all 3 corners inside are kept, any that have
          * all 3 outside are dropped, and the rest straddle the edge,
          * and need to be trimmed */

         keep_tri = FALSE;
         trim_tri = FALSE;

         for (i=0; i<3; i++) node_is_good[i] = TRUE;

         if (use_x_min)
            for (i=0; i<3; i++)
               if (the_tri->node[i]->loc.x < x_min)
                  node_is_good[i] = FALSE;

         if (use_x_max)
            for (i=0; i<3; i++)
               if (the_tri->node[i]->loc.x > x_max)
                  node_is_good[i] = FALSE;

         if (use_y_min)
            for (i=0; i<3; i++)
               if (the_tri->node[i]->loc.y < y_min)
                  node_is_good[i] = FALSE;

         if (use_y_max)
            for (i=0; i<3; i++)
               if (the_tri->node[i]->loc.y > y_max)
                  node_is_good[i] = FALSE;

         if (use_z_min)
            for (i=0; i<3; i++)
               if (the_tri->node[i]->loc.z < z_min)
                  node_is_good[i] = FALSE;

         if (use_z_max)
            for (i=0; i<3; i++)
               if (the_tri->node[i]->loc.z > z_max)
                  node_is_good[i] = FALSE;

         if (node_is_good[0] && node_is_good[1] && node_is_good[2]) {
            keep_tri = TRUE;
            trim_tri = FALSE;
         } else if (node_is_good[0] || node_is_good[1] || node_is_good[2]) {
            keep_tri = FALSE;
            trim_tri = TRUE;
         } else {
            keep_tri = FALSE;
            trim_tri = FALSE;
         }
      }

      /* check the minimum area criterium */
      if (use_a_min && (keep_tri || trim_tri)) {
         tri_area = find_area(the_tri);
         if (tri_area < a_min) {
            keep_tri = FALSE;
            /* fprintf(stderr,"\narea is %lf, under %lf",tri_area,a_min); */
         }
      }

      /* if the triangle survived the battery of tests, print it */
      if (keep_tri) {
         /* Keep the triangle wholly */
         write_tri(stdout,output_format,the_tri);
         num_wrote++;

      } else if (trim_tri) {
         /* Keep only part of the triangle */
         /* trim_tri(stdout,output_format,the_tri,node_is_good,num_wrote); */
         /* v3 is the index of the odd node */
         vn3 = -1;
         if (node_is_good[0] == node_is_good[1]) vn3 = 2;
         if (node_is_good[1] == node_is_good[2]) vn3 = 0;
         if (node_is_good[2] == node_is_good[0]) vn3 = 1;
         vn2 = mod(vn3+2,3);
         vn1 = mod(vn2+2,3);
         /* what are the two intersection points? */
         if (use_x_min) thresh = x_min;
         if (use_x_max) thresh = x_max;
         if (use_y_min) thresh = y_min;
         if (use_y_max) thresh = y_max;
         if (use_z_min) thresh = z_min;
         if (use_z_max) thresh = z_max;

         /* between node vn1 and vn3 */
         if (use_x_min || use_x_max) {
            frac = (thresh - the_tri->node[vn1]->loc.x) /
                   (the_tri->node[vn3]->loc.x - the_tri->node[vn1]->loc.x);
         } else if (use_y_min || use_y_max) {
            frac = (thresh - the_tri->node[vn1]->loc.y) /
                   (the_tri->node[vn3]->loc.y - the_tri->node[vn1]->loc.y);
         } else if (use_z_min || use_z_max) {
            frac = (thresh - the_tri->node[vn1]->loc.z) /
                   (the_tri->node[vn3]->loc.z - the_tri->node[vn1]->loc.z);
         }
         tnode1->loc.x = the_tri->node[vn1]->loc.x + frac*
                   (the_tri->node[vn3]->loc.x - the_tri->node[vn1]->loc.x);
         tnode1->loc.y = the_tri->node[vn1]->loc.y + frac*
                   (the_tri->node[vn3]->loc.y - the_tri->node[vn1]->loc.y);
         tnode1->loc.z = the_tri->node[vn1]->loc.z + frac*
                   (the_tri->node[vn3]->loc.z - the_tri->node[vn1]->loc.z);

         /* between node vn2 and vn3 */
         if (use_x_min || use_x_max) {
            frac = (thresh - the_tri->node[vn2]->loc.x) /
                   (the_tri->node[vn3]->loc.x - the_tri->node[vn2]->loc.x);
         } else if (use_y_min || use_y_max) {
            frac = (thresh - the_tri->node[vn2]->loc.y) /
                   (the_tri->node[vn3]->loc.y - the_tri->node[vn2]->loc.y);
         } else if (use_z_min || use_z_max) {
            frac = (thresh - the_tri->node[vn2]->loc.z) /
                   (the_tri->node[vn3]->loc.z - the_tri->node[vn2]->loc.z);
         }
         tnode2->loc.x = the_tri->node[vn2]->loc.x + frac*
                   (the_tri->node[vn3]->loc.x - the_tri->node[vn2]->loc.x);
         tnode2->loc.y = the_tri->node[vn2]->loc.y + frac*
                   (the_tri->node[vn3]->loc.y - the_tri->node[vn2]->loc.y);
         tnode2->loc.z = the_tri->node[vn2]->loc.z + frac*
                   (the_tri->node[vn3]->loc.z - the_tri->node[vn2]->loc.z);

         if (node_is_good[vn3]) {
            /* keep the single triangle tip */
            ttri1->node[0] = the_tri->node[vn3];
            ttri1->node[1] = tnode1;
            ttri1->node[2] = tnode2;

            /* check the minimum area criterium */
            if (use_a_min) {
               tri_area = find_area(ttri1);
               if (tri_area > a_min) {
                  write_tri(stdout,output_format,ttri1);
                  num_wrote++;
               }
            } else {
               write_tri(stdout,output_format,ttri1);
               num_wrote++;
            }

         } else {
            /* keep the 4-sided triangle base */
            ttri1->node[0] = the_tri->node[vn1];
            ttri1->node[1] = tnode2;
            ttri1->node[2] = tnode1;
            ttri2->node[0] = the_tri->node[vn2];
            ttri2->node[1] = tnode2;
            ttri2->node[2] = the_tri->node[vn1];

            /* check the minimum area criterium */
            if (use_a_min) {
               tri_area = find_area(ttri1);
               if (tri_area > a_min) {
                  write_tri(stdout,output_format,ttri1);
                  num_wrote++;
               }
               tri_area = find_area(ttri2);
               if (tri_area > a_min) {
                  write_tri(stdout,output_format,ttri2);
                  num_wrote++;
               }
            } else {
               write_tri(stdout,output_format,ttri1);
               write_tri(stdout,output_format,ttri2);
               num_wrote+=2;
            }
         }

      } else {
         /* nothing, the triangle is trimmed off */
      }

      if (num_read/DOTPER == (num_read+DPMO)/DOTPER) fprintf(stderr,".");
   }
   fprintf(stderr,"\n");

   fclose(ifp);

   fprintf(stderr,"Read %d triangles, wrote %d\n",num_read,num_wrote);
   /* fprintf(stderr,"Done.\n"); */
   exit(0);
}


/*
 * trim_tri will take 1 triangle and trim 1 or 2 corners off of it,
 * potentially making 2 or 1 triangle to write
 */
/*
int trim_tri(FILE *out,int output_format,tri_pointer the_tri,
   node_ptr the_nodes[3],int node_good[3],
   int num_wrote) {

   if (1==1) then {
      write_tri(out,output_format,the_tri);
   } else {
      write_tri(out,output_format,the_tri);
      write_tri(out,output_format,the_other_tri);
   }
}
*/


/*
 * This function writes basic usage information to stderr,
 * and then quits. Too bad.
 */
int Usage(char progname[80],int status) {

   /* Usage for rocktrim */
   static char **cpp, *help_message[] =
   {
       "where [-options] are one or more of the following:                         ",
       "                                                                           ",
       "   -i          inclusive cutting                                           ",
       "   -e          exclusive cutting                                           ",
       "   -s          exclusive cutting with smooth trimming                      ",
       "                                                                           ",
       "   -x val      trim off every triangle that exists completely below        ",
       "               this value of x, think of it as the minimum x value         ",
       "                                                                           ",
       "   +x val      trim off every triangle that exists completely above        ",
       "               this value of x                                             ",
       "                                                                           ",
       "               -y, +y, -z, and +z all work in the same way as -x and +x    ",
       "               All criteria must be met by a triangle in order to keep it  ",
       "                                                                           ",
       "   -ma val     trim all triangles with an area below the area given        ",
       "                                                                           ",
       "   -okey       specify output format, key= raw, rad, or tin                ",
       "               default = raw;                                              ",
       "                                                                           ",
       "   -help       (in place of infile) returns this help information          ",
       " ",
       "The input file can be of .raw, .tin, or .rad format, and the program requires",
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
