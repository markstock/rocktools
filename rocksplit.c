/*************************************************************
 *
 *  rocksplit.c - Smoothly split a trimesh into two files
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

#define SPLIT
#include "structs.h"

node_ptr node_head = NULL;
norm_ptr norm_head = NULL;
text_ptr text_head = NULL;

void find_arc_intersection (int, double, VEC, VEC, VEC, VEC, VEC*, VEC*);
int Usage(char[MAX_FN_LEN],int);


int main(int argc,char **argv) {

   int i;
   int input_format = 0;	// integer flag for input file type
   int output_format = 0;	// integer flag for output file type
   int num_read = 0;
   int num_wrote_1 = 0;		// number of triangles written to file 1
   int num_wrote_2 = 0;		// number of triangles written to file 2
   int num_tris = 0;
   int num_nodes = 0;
   //int have_normals = 0;
   int tri_is_low;		// keep the current tri or not?
   int trim_tri;		// trim the current tri or not?
   int node_is_low[3];
   double thresh = 0.;		// the cutoff threshhold for whatever axis
   double frac = 0.;		// more dealing with the cutoff
   double x_split = 0.0;
   double y_split = 0.0;
   double z_split = 0.0;
   VEC bmin,bmax,cm;
   enum use_dir_type {
      pick_longest,
      x,
      y,
      z } use_dir = pick_longest;
   int vn1,vn2,vn3;
   int use_given_root = FALSE;
   char infile[MAX_FN_LEN];		/* name of input file */
   char extension[4];		/* filename extension if infile */
   //char axis_char1 = 'x';
   //char axis_char2 = 'X';
   char output_string[4] = "raw"; /* format extension for the output */
   char output_root[MAX_FN_LEN-5];	/* filename root for the output */
   char output_1[MAX_FN_LEN];		/* filename for the output */
   char output_2[MAX_FN_LEN];		/* filename for the output */
   char progname[MAX_FN_LEN];		/* name of binary executable */
   tri_pointer the_tri,ttri1,ttri2;
   node_ptr the_nodes[3],tnode1,tnode2;
   VEC tnorm1,tnorm2;
   FILE *ifp;
   FILE *ofp1;
   FILE *ofp2;


   /* Parse command-line args */
   (void) strcpy(progname,argv[0]);
   if (argc < 2) (void) Usage(progname,0);
   if (strncmp(argv[1], "-help", 2) == 0)
      (void) Usage(progname,0);
   (void) strcpy(infile,argv[1]);
   for (i=2; i<argc; i++) {
      if (strncmp(argv[i], "-x", 2) == 0) {
         x_split = atof(argv[++i]);
         use_dir = x;
      } else if (strncmp(argv[i], "-y", 2) == 0) {
         y_split = atof(argv[++i]);
         use_dir = x;
      } else if (strncmp(argv[i], "-z", 2) == 0) {
         z_split = atof(argv[++i]);
         use_dir = x;
      } else if (strncmp(argv[i], "-l", 2) == 0) {
         use_dir = pick_longest;
      } else if (strncmp(argv[i], "-o", 2) == 0) {
         strncpy(output_string,argv[i]+2,3);
      } else if (strncmp(argv[i], "-root", 2) == 0) {
         strcpy(output_root,argv[++i]);
         use_given_root = TRUE;
      } else
         (void) Usage(progname,0);
   }

   // if an output filename root has not been given, determine it
   if (!use_given_root) {
      strncpy(output_root,infile,strlen(output_root));
      output_root[strlen(infile)-4] = '\0';
      //printf("output root is (%s)\n",output_root);
   }

   /* Determine the input file format from the .XXX extension, and read it */
   strncpy(extension,infile+strlen(infile)-3,3);
   if (strncmp(extension, "raw", 3) == 0)
      input_format = 1;
   else if (strncmp(extension, "tin", 1) == 0)
      input_format = 2;
   else if (strncmp(extension, "rad", 1) == 0)
      input_format = 3;
   else {
      fprintf(stderr,"Input filename extension is assumed to be (%s)\n",extension);
      fprintf(stderr,"This is either not a supported input file format for rocksplit, or you\n");
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
      fprintf(stderr,"This is either not a supported output file format for rocksplit, or you\n");
      fprintf(stderr,"   need to use a proper extension (standard 3-character extension for format).\n");
      fprintf(stderr,"Supported output file formats are: .raw, .tin, and .rad\n");
      exit(0);
   }

   // determine the longest side ------------------------------------------
   if (use_dir == pick_longest) {
      // call a routine that scans the file and returns the edge lengths
      float tempf;
      find_mesh_stats(infile,&bmin,&bmax,FALSE,&cm,&tempf,&num_tris,&num_nodes);

      // then, compare them to find the splitting direction and value
      if (bmax.x-bmin.x > bmax.y-bmin.y) {
         if (bmax.z-bmin.z > bmax.x-bmin.x) {
            use_dir = z;
            z_split = 0.5*(bmax.z+bmin.z);
         } else {
            use_dir = x;
            x_split = 0.5*(bmax.x+bmin.x);
         }
      } else {
         if (bmax.z-bmin.z > bmax.y-bmin.y) {
            use_dir = z;
            z_split = 0.5*(bmax.z+bmin.z);
         } else {
            use_dir = y;
            y_split = 0.5*(bmax.y+bmin.y);
         }
      }

      if (num_tris == 0) exit(0);
   }

   // determine the output file names
   if (use_dir == x) {
      //axis_char1 = 'l';
      //axis_char2 = 'X';
      fprintf(stderr,"Splitting at x = %g\n",x_split);
   } else if (use_dir == y) {
      //axis_char1 = 'l';
      //axis_char2 = 'Y';
      fprintf(stderr,"Splitting at y = %g\n",y_split);
   } else {
      //axis_char1 = 'l';
      //axis_char2 = 'Z';
      fprintf(stderr,"Splitting at z = %g\n",z_split);
   }
   if (output_format == 1) {
      strcpy(extension,"raw");
   } else if (output_format == 1) {
      strcpy(extension,"tin");
   } else {
      strcpy(extension,"rad");
   }

   sprintf(output_1,"%s%c.%s",output_root,'l',extension);
   sprintf(output_2,"%s%c.%s",output_root,'r',extension);

   //fprintf(stdout,"outfile 1 is (%s)\n",output_1);
   //fprintf(stdout,"outfile 2 is (%s)\n",output_2);
   //exit(0);

   // finally, open the files and begin the splitting -----------------------

   ofp1 = fopen(output_1,"w");
   if (ofp1==NULL) {
      fprintf(stderr,"Could not open output file %s\n",output_1);
      exit(0);
   }
   fprintf(stdout,"Opening file %s for writing.\n",output_1);
   fflush(stdout);

   ofp2 = fopen(output_2,"w");
   if (ofp2==NULL) {
      fprintf(stderr,"Could not open output file %s\n",output_2);
      exit(0);
   }
   fprintf(stdout,"Opening file %s for writing.\n",output_2);
   fflush(stdout);

   //fprintf(ofp1,"file 1\n");
   //fprintf(ofp2,"file 2\n");
   //fclose(ofp1);
   //fclose(ofp2);
   //exit(0);

   /* Set up memory space for the working triangle, and the three nodes */
   the_tri = alloc_new_tri();
   for (i=0; i<3; i++) {
      the_nodes[i] = (NODE *)malloc(sizeof(NODE));
      the_tri->node[i] = the_nodes[i];
   }
   /* these are for the nodes and triangles in case a triangle needs a corner trimmed */
   tnode1 = (NODE *)malloc(sizeof(NODE));
   tnode2 = (NODE *)malloc(sizeof(NODE));
   ttri1 = alloc_new_tri();
   for (i=0; i<3; i++) ttri1->norm[i] = (NORM *)malloc(sizeof(NORM));
   ttri2 = alloc_new_tri();
   for (i=0; i<3; i++) ttri2->norm[i] = (NORM *)malloc(sizeof(NORM));
   tnorm1.x = 0.;
   tnorm1.y = 0.;
   tnorm1.z = 0.;
   tnorm2.x = 0.;
   tnorm2.y = 0.;
   tnorm2.z = 0.;

   // Read the input file ---------------------------------------------------

   /* open the file for reading */
   ifp = fopen(infile,"r");
   if (ifp==NULL) {
      fprintf(stderr,"Could not open input file %s\n",infile);
      exit(0);
   }
   fprintf(stderr,"Opening file %s, splitting",infile);
   fflush(stderr);


   /* as long as there are triangles available, operate */
   while (get_tri(ifp,input_format,the_tri) == 1) {

      num_read++;
      tri_is_low = FALSE;
      trim_tri = FALSE;

      for (i=0; i<3; i++) node_is_low[i] = TRUE;

      if (use_dir == x)
         for (i=0; i<3; i++)
            if (the_tri->node[i]->loc.x > x_split)
               node_is_low[i] = FALSE;

      if (use_dir == y)
         for (i=0; i<3; i++)
            if (the_tri->node[i]->loc.y > y_split)
               node_is_low[i] = FALSE;

      if (use_dir == z)
         for (i=0; i<3; i++)
            if (the_tri->node[i]->loc.z > z_split)
               node_is_low[i] = FALSE;

      if (node_is_low[0] && node_is_low[1] && node_is_low[2]) {
         tri_is_low = TRUE;
         trim_tri = FALSE;
      } else if (!node_is_low[0] && !node_is_low[1] && !node_is_low[2]) {
         tri_is_low = FALSE;
         trim_tri = FALSE;
      } else {
         tri_is_low = FALSE;
         trim_tri = TRUE;
      }

      // put the triangle where it belongs
      if (trim_tri) {

         // parts of the triangle go to each side

         // v3 is the index of the odd node
         vn3 = -1;
         if (node_is_low[0] == node_is_low[1]) vn3 = 2;
         if (node_is_low[1] == node_is_low[2]) vn3 = 0;
         if (node_is_low[2] == node_is_low[0]) vn3 = 1;
         vn2 = mod(vn3+2,3);
         vn1 = mod(vn2+2,3);

         // what are the two intersection points?
         if (use_dir == x) { thresh = x_split; i = 0; }
         if (use_dir == y) { thresh = y_split; i = 1; }
         if (use_dir == z) { thresh = z_split; i = 2; }

         // between node vn1 and vn3
         if (the_tri->norm[vn1] && the_tri->norm[vn3]) {
            find_arc_intersection(i,thresh,
                                  the_tri->node[vn1]->loc, the_tri->norm[vn1]->norm,
                                  the_tri->node[vn3]->loc, the_tri->norm[vn3]->norm,
                                  &tnode1->loc,            &tnorm1);
         } else {
            if (use_dir == x) {
               frac = (thresh - the_tri->node[vn1]->loc.x) /
                      (the_tri->node[vn3]->loc.x - the_tri->node[vn1]->loc.x);
            } else if (use_dir == y) {
               frac = (thresh - the_tri->node[vn1]->loc.y) /
                      (the_tri->node[vn3]->loc.y - the_tri->node[vn1]->loc.y);
            } else if (use_dir == z) {
               frac = (thresh - the_tri->node[vn1]->loc.z) /
                      (the_tri->node[vn3]->loc.z - the_tri->node[vn1]->loc.z);
            }
            tnode1->loc.x = the_tri->node[vn1]->loc.x + frac*
                      (the_tri->node[vn3]->loc.x - the_tri->node[vn1]->loc.x);
            tnode1->loc.y = the_tri->node[vn1]->loc.y + frac*
                      (the_tri->node[vn3]->loc.y - the_tri->node[vn1]->loc.y);
            tnode1->loc.z = the_tri->node[vn1]->loc.z + frac*
                      (the_tri->node[vn3]->loc.z - the_tri->node[vn1]->loc.z);
            //tnorm1.x = the_tri->norm[vn1].x + frac*
            //          (the_tri->norm[vn3].x - the_tri->norm[vn1].x);
            //tnorm1.y = the_tri->norm[vn1].y + frac*
            //          (the_tri->norm[vn3].y - the_tri->norm[vn1].y);
            //tnorm1.z = the_tri->norm[vn1].z + frac*
            //          (the_tri->norm[vn3].z - the_tri->norm[vn1].z);
         }

         // between node vn2 and vn3
         if (the_tri->norm[vn2] && the_tri->norm[vn3]) {
            find_arc_intersection(i,thresh,
                                  the_tri->node[vn2]->loc, the_tri->norm[vn2]->norm,
                                  the_tri->node[vn3]->loc, the_tri->norm[vn3]->norm,
                                  &tnode2->loc,            &tnorm2);
         } else {
            if (use_dir == x) {
               frac = (thresh - the_tri->node[vn2]->loc.x) /
                      (the_tri->node[vn3]->loc.x - the_tri->node[vn2]->loc.x);
            } else if (use_dir == y) {
               frac = (thresh - the_tri->node[vn2]->loc.y) /
                      (the_tri->node[vn3]->loc.y - the_tri->node[vn2]->loc.y);
            } else if (use_dir == z) {
               frac = (thresh - the_tri->node[vn2]->loc.z) /
                      (the_tri->node[vn3]->loc.z - the_tri->node[vn2]->loc.z);
            }
            tnode2->loc.x = the_tri->node[vn2]->loc.x + frac*
                      (the_tri->node[vn3]->loc.x - the_tri->node[vn2]->loc.x);
            tnode2->loc.y = the_tri->node[vn2]->loc.y + frac*
                      (the_tri->node[vn3]->loc.y - the_tri->node[vn2]->loc.y);
            tnode2->loc.z = the_tri->node[vn2]->loc.z + frac*
                      (the_tri->node[vn3]->loc.z - the_tri->node[vn2]->loc.z);
            //tnorm2.x = the_tri->norm[vn2].x + frac*
            //          (the_tri->norm[vn3].x - the_tri->norm[vn2].x);
            //tnorm2.y = the_tri->norm[vn2].y + frac*
            //          (the_tri->norm[vn3].y - the_tri->norm[vn2].y);
            //tnorm2.z = the_tri->norm[vn2].z + frac*
            //          (the_tri->norm[vn3].z - the_tri->norm[vn2].z);
         }

         if (node_is_low[vn3]) {
            // the single triangle tip goes to file 1
            ttri1->node[0] = the_tri->node[vn3];
            ttri1->node[1] = tnode1;
            ttri1->node[2] = tnode2;
            if (the_tri->norm[0] && the_tri->norm[1] && the_tri->norm[2]) {
               ttri1->norm[0]->norm.x = the_tri->norm[vn3]->norm.x;
               ttri1->norm[0]->norm.y = the_tri->norm[vn3]->norm.y;
               ttri1->norm[0]->norm.z = the_tri->norm[vn3]->norm.z;
               ttri1->norm[1]->norm.x = tnorm1.x;
               ttri1->norm[1]->norm.y = tnorm1.y;
               ttri1->norm[1]->norm.z = tnorm1.z;
               ttri1->norm[2]->norm.x = tnorm2.x;
               ttri1->norm[2]->norm.y = tnorm2.y;
               ttri1->norm[2]->norm.z = tnorm2.z;
            }

            write_tri(ofp1,output_format,ttri1);
            num_wrote_1++;

            // the 4-sided triangle base goes to file 2
            ttri1->node[0] = the_tri->node[vn1];
            ttri1->node[1] = tnode2;
            ttri1->node[2] = tnode1;
            ttri2->node[0] = the_tri->node[vn2];
            ttri2->node[1] = tnode2;
            ttri2->node[2] = the_tri->node[vn1];
            if (the_tri->norm[0] && the_tri->norm[1] && the_tri->norm[2]) {
               ttri1->norm[0]->norm.x = the_tri->norm[vn1]->norm.x;
               ttri1->norm[0]->norm.y = the_tri->norm[vn1]->norm.y;
               ttri1->norm[0]->norm.z = the_tri->norm[vn1]->norm.z;
               ttri1->norm[1]->norm.x = tnorm2.x;
               ttri1->norm[1]->norm.y = tnorm2.y;
               ttri1->norm[1]->norm.z = tnorm2.z;
               ttri1->norm[2]->norm.x = tnorm1.x;
               ttri1->norm[2]->norm.y = tnorm1.y;
               ttri1->norm[2]->norm.z = tnorm1.z;
               ttri2->norm[0]->norm.x = the_tri->norm[vn2]->norm.x;
               ttri2->norm[0]->norm.y = the_tri->norm[vn2]->norm.y;
               ttri2->norm[0]->norm.z = the_tri->norm[vn2]->norm.z;
               ttri2->norm[1]->norm.x = tnorm2.x;
               ttri2->norm[1]->norm.y = tnorm2.y;
               ttri2->norm[1]->norm.z = tnorm2.z;
               ttri2->norm[2]->norm.x = the_tri->norm[vn1]->norm.x;
               ttri2->norm[2]->norm.y = the_tri->norm[vn1]->norm.y;
               ttri2->norm[2]->norm.z = the_tri->norm[vn1]->norm.z;
            }

            write_tri(ofp2,output_format,ttri1);
            write_tri(ofp2,output_format,ttri2);
            num_wrote_2 += 2;

         } else {
            // the 4-sided triangle base goes to file 1
            ttri1->node[0] = the_tri->node[vn1];
            ttri1->node[1] = tnode2;
            ttri1->node[2] = tnode1;
            ttri2->node[0] = the_tri->node[vn2];
            ttri2->node[1] = tnode2;
            ttri2->node[2] = the_tri->node[vn1];
            if (the_tri->norm[0] && the_tri->norm[1] && the_tri->norm[2]) {
               ttri1->norm[0]->norm.x = the_tri->norm[vn1]->norm.x;
               ttri1->norm[0]->norm.y = the_tri->norm[vn1]->norm.y;
               ttri1->norm[0]->norm.z = the_tri->norm[vn1]->norm.z;
               ttri1->norm[1]->norm.x = tnorm2.x;
               ttri1->norm[1]->norm.y = tnorm2.y;
               ttri1->norm[1]->norm.z = tnorm2.z;
               ttri1->norm[2]->norm.x = tnorm1.x;
               ttri1->norm[2]->norm.y = tnorm1.y;
               ttri1->norm[2]->norm.z = tnorm1.z;
               ttri2->norm[0]->norm.x = the_tri->norm[vn2]->norm.x;
               ttri2->norm[0]->norm.y = the_tri->norm[vn2]->norm.y;
               ttri2->norm[0]->norm.z = the_tri->norm[vn2]->norm.z;
               ttri2->norm[1]->norm.x = tnorm2.x;
               ttri2->norm[1]->norm.y = tnorm2.y;
               ttri2->norm[1]->norm.z = tnorm2.z;
               ttri2->norm[2]->norm.x = the_tri->norm[vn1]->norm.x;
               ttri2->norm[2]->norm.y = the_tri->norm[vn1]->norm.y;
               ttri2->norm[2]->norm.z = the_tri->norm[vn1]->norm.z;
            }

            write_tri(ofp1,output_format,ttri1);
            write_tri(ofp1,output_format,ttri2);
            num_wrote_1 += 2;

            // the single triangle tip goes to file 2
            ttri1->node[0] = the_tri->node[vn3];
            ttri1->node[1] = tnode1;
            ttri1->node[2] = tnode2;
            if (the_tri->norm[0] && the_tri->norm[1] && the_tri->norm[2]) {
               ttri1->norm[0]->norm.x = the_tri->norm[vn3]->norm.x;
               ttri1->norm[0]->norm.y = the_tri->norm[vn3]->norm.y;
               ttri1->norm[0]->norm.z = the_tri->norm[vn3]->norm.z;
               ttri1->norm[1]->norm.x = tnorm1.x;
               ttri1->norm[1]->norm.y = tnorm1.y;
               ttri1->norm[1]->norm.z = tnorm1.z;
               ttri1->norm[2]->norm.x = tnorm2.x;
               ttri1->norm[2]->norm.y = tnorm2.y;
               ttri1->norm[2]->norm.z = tnorm2.z;
            }

            write_tri(ofp2,output_format,ttri1);
            num_wrote_2++;
         }

      } else if (tri_is_low) {

         // write the whole triangle to output_1
         write_tri(ofp1,output_format,the_tri);
         num_wrote_1++;

      } else {

         // write the whole triangle to output_2
         write_tri(ofp2,output_format,the_tri);
         num_wrote_2++;
      }

      if (num_read/DOTPER == (num_read+DPMO)/DOTPER) fprintf(stderr,".");
   }
   fprintf(stderr,"\n");

   fclose(ifp);
   fclose(ofp1);
   fclose(ofp2);

   fprintf(stderr,"Read %d triangles, wrote %d and %d\n",num_read,num_wrote_1,num_wrote_2);
   /* fprintf(stderr,"Done.\n"); */
   exit(0);
}


/*
 * find_arc_intersection uses spline interpolation and a Newton solver
 * to determine the intersection of the ideal, smoothed edge with the
 * cutting plane
 */
void find_arc_intersection (int dim, double thresh, VEC x1, VEC n1,
                            VEC x2, VEC n2, VEC *xp, VEC *np) {

  int i,j;
  int cnt = 0;
  double frac,error, lfrac,lerror, ufrac,uerror;
  //VEC planenorm,planeinter;
  double dl[3],fp[2][3][3],p1[3],p2[3],a[4];
  double xp1[3],xp2[3],norm1[3],norm2[3];

  // normalize normals
  n1 = norm(n1);
  n2 = norm(n2);

  // convert "VEC" to array
  xp1[0] = x1.x;
  xp1[1] = x1.y;
  xp1[2] = x1.z;
  norm1[0] = n1.x;
  norm1[1] = n1.y;
  norm1[2] = n1.z;
  xp2[0] = x2.x;
  xp2[1] = x2.y;
  xp2[2] = x2.z;
  norm2[0] = n2.x;
  norm2[1] = n2.y;
  norm2[2] = n2.z;

  // set up error vector
  //planenorm.x = 0.;
  //planenorm.y = 0.;
  //planenorm.z = 0.;
  //if (dim==0) planenorm.x = 1.;
  //if (dim==1) planenorm.y = 1.;
  //if (dim==2) planenorm.z = 1.;
  //planeinter.x = 0.;
  //planeinter.y = 0.;
  //planeinter.z = 0.;
  //if (dim==0) planeinter.x = thresh;
  //if (dim==1) planeinter.y = thresh;
  //if (dim==2) planeinter.z = thresh;
  //fprintf(stderr,"\nx1 %g %g %g  n1 %g %g %g\n",x1.x,x1.y,x1.z,n1.x,n1.y,n1.z);
  //fprintf(stderr,"x2 %g %g %g  n2 %g %g %g\n",x2.x,x2.y,x2.z,n2.x,n2.y,n2.z);

  // compute the tangential operator for each node (P = I - nn^T)
  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      fp[0][i][j] = 0.;
      fp[1][i][j] = 0.;
    }
    fp[0][i][i] = 1.;
    fp[1][i][i] = 1.;
    for (j=0; j<3; j++) {
      fp[0][i][j] -= norm1[i]*norm1[j];
      fp[1][i][j] -= norm2[i]*norm2[j];
    }
  }
  // find the vector product of each of these with dl
  for (i=0; i<3; i++) {
    p1[i] = 0.;
    p2[i] = 0.;
  }
  for (i=0; i<3; i++) dl[i] = xp2[i] - xp1[i];
  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      p1[i] += dl[j]*fp[0][i][j];
      p2[i] += dl[j]*fp[1][i][j];
    }
  }
  // calculate the spline for only the splitting plane normal axis
  // a[0] = f_0
  a[0] = xp1[dim];
  // a[1] = f'_0
  a[1] = p1[dim];
  // a[2] = 3 (f_1 - f_0) / h^2 - (f'_1 + 2 f'_0) / h  (h=1 for us)
  a[2] = 3. * (xp2[dim] - xp1[dim]) - (p2[dim] + 2. * p1[dim]);
  // a[3] = 2 (f_0 - f_1) / h^3 + (f'_1 + f'_0) / h^2
  a[3] = 2. * (xp1[dim] - xp2[dim]) + (p1[dim] + p2[dim]);


  // first guesses
  frac = 0.5;
  lfrac = 0.0;
  lerror = a[0] + a[1]*lfrac + a[2]*lfrac*lfrac + a[3]*lfrac*lfrac*lfrac - thresh;
  //fprintf(stderr,"guess %g has error %g\n",lfrac,lerror);

  ufrac = 1.0;
  uerror = a[0] + a[1]*ufrac + a[2]*ufrac*ufrac + a[3]*ufrac*ufrac*ufrac - thresh;
  //fprintf(stderr,"guess %g has error %g\n",ufrac,uerror);

  if (lerror*uerror > 0.) {
    fprintf(stderr,"degenerate edge?\n");
    frac = (thresh - xp1[dim]) / (xp2[dim] - xp1[dim]);
    xp->x = xp1[0] + frac*(xp2[0]-xp1[0]);
    xp->y = xp1[1] + frac*(xp2[1]-xp1[1]);
    xp->z = xp1[2] + frac*(xp2[2]-xp1[2]);
    np->x = norm1[0] + frac*(norm2[0]-norm1[0]);
    np->y = norm1[1] + frac*(norm2[1]-norm1[1]);
    np->z = norm1[2] + frac*(norm2[2]-norm1[2]);
    *np = norm(*np);
    return;
  }

  // loop until the error is tiny
  error = 1.;
  while (fabs(error) > 1.e-7 && ++cnt < 100) {

    // take new guess
    frac = lfrac - lerror*(ufrac-lfrac)/(uerror-lerror);

    // find error
    error = a[0] + a[1]*frac + a[2]*frac*frac + a[3]*frac*frac*frac - thresh;
    //fprintf(stderr,"guess %g has error %g (%g %g)\n",frac,error,lerror,uerror);

    // replace upper or lower bound
    if (error*lerror > 0.) {
      // current and lower errors have same sign
      lerror = error;
      lfrac = frac;
    } else {
      uerror = error;
      ufrac = frac;
    }

  }

  // if we finished poorly, take a lame guess and quit
  if (cnt > 100) {
    fprintf(stderr,"degenerate edge?\n");
    frac = (thresh - xp1[dim]) / (xp2[dim] - xp1[dim]);
    xp->x = xp1[0] + frac*(xp2[0]-xp1[0]);
    xp->y = xp1[1] + frac*(xp2[1]-xp1[1]);
    xp->z = xp1[2] + frac*(xp2[2]-xp1[2]);
    np->x = norm1[0] + frac*(norm2[0]-norm1[0]);
    np->y = norm1[1] + frac*(norm2[1]-norm1[1]);
    np->z = norm1[2] + frac*(norm2[2]-norm1[2]);
    *np = norm(*np);
    return;
  }

  // if we finished well, use frac to find the full 3D point
  for (i=0; i<3; i++) {
    a[0] = xp1[i];
    a[1] = p1[i];
    a[2] = 3. * (xp2[i] - xp1[i]) - (p2[i] + 2. * p1[i]);
    a[3] = 2. * (xp1[i] - xp2[i]) + (p1[i] + p2[i]);
    xp1[i] = a[0] + a[1]*frac + a[2]*frac*frac + a[3]*frac*frac*frac;
    norm1[i] = norm1[i] + frac*(norm2[i]-norm1[i]);
  }
  xp->x = xp1[0];
  xp->y = xp1[1];
  xp->z = xp1[2];
  np->x = norm1[0];
  np->y = norm1[1];
  np->z = norm1[2];
  // normalize normal
  *np = norm(*np);

  //fprintf(stderr,"xp %g %g %g  np %g %g %g\n",xp->x,xp->y,xp->z,np->x,np->y,np->z);

  //exit(0);
  return;
}


/*
 * This function writes basic usage information to stderr,
 * and then quits. Too bad.
 */
int Usage(char progname[MAX_FN_LEN],int status) {

   /* Usage for rocksplit */
   static char **cpp, *help_message[] =
   {
       "where [-options] are one or more of the following:                         ",
       "                                                                           ",
       "   -x val      split the triangle mesh into components completely existing ",
       "               above and below the given value for x                       ",
       "                                                                           ",
       "   -y val      same as for -x, but in y-direction                          ",
       "                                                                           ",
       "   -z val      same as for -x, but in z-direction                          ",
       "                                                                           ",
       "   -l          Split along middle of longest side, whether x, y, or z      ",
       "                                                                           ",
       "   -root name  use 'name' instead of input file root as root of new files  ",
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
       " ",
       "Program will write two files, only runtime information goes to stdout.",
       "  Example:",
       "    rocksplit test.raw -x 0.0 -oraw",
       "  will produce two files: testl.raw and testr.raw",
       " ",
       NULL
   };

   fprintf(stderr, "usage:\n  %s infile [-options]\n\n", progname);
   for (cpp = help_message; *cpp; cpp++) fprintf(stderr, "%s\n", *cpp);
   fflush(stderr);
   exit(status);
   return(0);
}
