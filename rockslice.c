/*************************************************************
 *
 *  rockslice.c - Cut a 2D cross-section from a trimesh
 *
 *  Mark J. Stock, mstock@umich.edu
 *
 *
 * rocktools - Tools for creating and manipulating triangular meshes
 * Copyright (C) 1999,2004,2006-8,2013-4  Mark J. Stock
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
int Usage(char[80],int);
int write_seg_slice(FILE*, char*, double, node_ptr, node_ptr, VEC*, VEC*);

// from inout.c
extern int find_mesh_stats(char*, VEC*, VEC*, int, VEC*, int*, int*);

int main(int argc,char **argv) {

   int i;
   int num_read = 0;
   int num_wrote_1 = 0;		// number of triangles written to file 1
   int num_tris = 0;
   int num_nodes = 0;
   int trim_tri;		// trim the current tri or not?
   int node_is_low[3];
   double thresh = 0.;		// the cutoff threshhold for whatever axis
   double frac = 0.;		// more dealing with the cutoff
   double split_val = 0.0;
   double line_width = 0.01;
   VEC bmin,bmax,cm;
   enum use_dir_type {
      pick_shortest,
      x,
      y,
      z } use_dir = pick_shortest;
   int vn1,vn2,vn3;
   char infile[255];		/* name of input file */
   char extension[4];		/* filename extension if infile */
   char output_format[4] = "raw"; /* format extension for the output */
   //char output_root[255];	/* filename root for the output */
   //char output_1[255];		/* filename for the output */
   char progname[255];		/* name of binary executable */
   //char normal_string[512];
   tri_pointer curr,tri_head;
   node_ptr tnode1,tnode2;
   VEC tnorm1,tnorm2;
   FILE *ofp1;


   /* Parse command-line args */
   (void) strcpy(progname,argv[0]);
   if (argc < 2) (void) Usage(progname,0);
   if (strncmp(argv[1], "-help", 2) == 0)
      (void) Usage(progname,0);
   (void) strcpy(infile,argv[1]);
   for (i=2; i<argc; i++) {
      if (strncmp(argv[i], "-x", 2) == 0) {
         split_val = atof(argv[++i]);
         use_dir = x;
      } else if (strncmp(argv[i], "-y", 2) == 0) {
         split_val = atof(argv[++i]);
         use_dir = y;
      } else if (strncmp(argv[i], "-z", 2) == 0) {
         split_val = atof(argv[++i]);
         use_dir = z;
      } else if (strncmp(argv[i], "-s", 2) == 0) {
         use_dir = pick_shortest;
      } else if (strncmp(argv[i], "-r", 2) == 0) {
         line_width = atof(argv[++i]);
      } else if (strncmp(argv[i], "-o", 2) == 0) {
         strncpy(output_format,argv[i]+2,4);
      } else
         (void) Usage(progname,0);
   }

   /* Determine the input file format from the .XXX extension, and read it */
   strncpy(extension,infile+strlen(infile)-3,4);

   /* Determine and set the output format key from the string */
   if (strncmp(output_format, "raw", 3) != 0 &&
       strncmp(output_format, "seg", 3) != 0) {
      fprintf(stderr,"Output filename extension is assumed to be (%s)\n",output_format);
      fprintf(stderr,"This is either not a supported output file format for rockslice, or you\n");
      fprintf(stderr,"   need to use a proper extension (standard 3-character extension for format).\n");
      fprintf(stderr,"Supported output file formats are: .raw and .seg\n");
      exit(0);
   }
   //sprintf(output_1,"%s%c.%s",output_root,'l',output_format);

   // determine the longest side ------------------------------------------
   if (use_dir == pick_shortest) {
      // call a routine that scans the file and returns the edge lengths
      find_mesh_stats(infile,&bmin,&bmax,FALSE,&cm,&num_tris,&num_nodes);

      // then, compare them to find the splitting direction and value
      if (bmax.x-bmin.x < bmax.y-bmin.y) {
         if (bmax.z-bmin.z < bmax.x-bmin.x) {
            use_dir = z;
            split_val = 0.5*(bmax.z+bmin.z);
         } else {
            use_dir = x;
            split_val = 0.5*(bmax.x+bmin.x);
         }
      } else {
         if (bmax.z-bmin.z < bmax.y-bmin.y) {
            use_dir = z;
            split_val = 0.5*(bmax.z+bmin.z);
         } else {
            use_dir = y;
            split_val = 0.5*(bmax.y+bmin.y);
         }
      }

      if (num_tris == 0) exit(0);
   }

   if (use_dir == x) {
      fprintf(stderr,"Generating slice at x = %g\n",split_val);
   } else if (use_dir == y) {
      fprintf(stderr,"Generating slice at y = %g\n",split_val);
   } else {
      fprintf(stderr,"Generating slice at z = %g\n",split_val);
   }


   // finally, open the files and begin the splitting -----------------------

   ofp1 = stdout;
   //ofp1 = fopen(output_1,"w");
   //if (ofp1==NULL) {
   //   fprintf(stderr,"Could not open output file %s\n",output_1);
   //   exit(0);
   //}
   //fprintf(stdout,"Opening file %s for writing.\n",output_1);
   fflush(stdout);

   /* these are for the nodes and normals */
   tnode1 = (NODE *)malloc(sizeof(NODE));
   tnode2 = (NODE *)malloc(sizeof(NODE));
   tnorm1.x = 0.;
   tnorm1.y = 0.;
   tnorm1.z = 0.;
   tnorm2.x = 0.;
   tnorm2.y = 0.;
   tnorm2.z = 0.;

   // Read the input file ---------------------------------------------------

   tri_head = read_input(infile,FALSE,NULL);

   /* as long as there are triangles available, operate */
   curr = tri_head;
   while (curr) {

      num_read++;
      trim_tri = FALSE;

      for (i=0; i<3; i++) node_is_low[i] = TRUE;

      if (use_dir == x) {
         for (i=0; i<3; i++)
            if (curr->node[i]->loc.x > split_val)
               node_is_low[i] = FALSE;

      } else if (use_dir == y) {
         for (i=0; i<3; i++)
            if (curr->node[i]->loc.y > split_val)
               node_is_low[i] = FALSE;

      } else if (use_dir == z) {
         for (i=0; i<3; i++)
            if (curr->node[i]->loc.z > split_val)
               node_is_low[i] = FALSE;
      }

      if (node_is_low[0] && node_is_low[1] && node_is_low[2]) {
         trim_tri = FALSE;
      } else if (!node_is_low[0] && !node_is_low[1] && !node_is_low[2]) {
         trim_tri = FALSE;
      } else {
         trim_tri = TRUE;
      }

      // this triangle intersects the slice plane, dump a segment
      if (trim_tri) {

         // v3 is the index of the odd node
         vn3 = -1;
         if (node_is_low[0] == node_is_low[1]) vn3 = 2;
         if (node_is_low[1] == node_is_low[2]) vn3 = 0;
         if (node_is_low[2] == node_is_low[0]) vn3 = 1;
         vn2 = mod(vn3+2,3);
         vn1 = mod(vn2+2,3);

         // what are the two intersection points?
         if (use_dir == x) { thresh = split_val; i = 0; }
         if (use_dir == y) { thresh = split_val; i = 1; }
         if (use_dir == z) { thresh = split_val; i = 2; }

         if (curr->norm[0] && curr->norm[1] && curr->norm[2]) {
            // use splines
            // between node vn1 and vn3
            find_arc_intersection(i,thresh,
                                  curr->node[vn1]->loc, curr->norm[vn1]->norm,
                                  curr->node[vn3]->loc, curr->norm[vn3]->norm,
                                  &tnode1->loc,         &tnorm1);
            // between node vn2 and vn3
            find_arc_intersection(i,thresh,
                                  curr->node[vn2]->loc, curr->norm[vn2]->norm,
                                  curr->node[vn3]->loc, curr->norm[vn3]->norm,
                                  &tnode2->loc,         &tnorm2);
         } else {
            // use midpoint
            // between node vn1 and vn3
            if (use_dir == x) {
               frac = (thresh - curr->node[vn1]->loc.x) /
                      (curr->node[vn3]->loc.x - curr->node[vn1]->loc.x);
            } else if (use_dir == y) {
               frac = (thresh - curr->node[vn1]->loc.y) /
                      (curr->node[vn3]->loc.y - curr->node[vn1]->loc.y);
            } else {
               frac = (thresh - curr->node[vn1]->loc.z) /
                      (curr->node[vn3]->loc.z - curr->node[vn1]->loc.z);
            }
            tnode1->loc.x = curr->node[vn1]->loc.x + frac*
                      (curr->node[vn3]->loc.x - curr->node[vn1]->loc.x);
            tnode1->loc.y = curr->node[vn1]->loc.y + frac*
                      (curr->node[vn3]->loc.y - curr->node[vn1]->loc.y);
            tnode1->loc.z = curr->node[vn1]->loc.z + frac*
                      (curr->node[vn3]->loc.z - curr->node[vn1]->loc.z);

            // between node vn2 and vn3
            if (use_dir == x) {
               frac = (thresh - curr->node[vn2]->loc.x) /
                      (curr->node[vn3]->loc.x - curr->node[vn2]->loc.x);
            } else if (use_dir == y) {
               frac = (thresh - curr->node[vn2]->loc.y) /
                      (curr->node[vn3]->loc.y - curr->node[vn2]->loc.y);
            } else if (use_dir == z) {
               frac = (thresh - curr->node[vn2]->loc.z) /
                      (curr->node[vn3]->loc.z - curr->node[vn2]->loc.z);
            }
            tnode2->loc.x = curr->node[vn2]->loc.x + frac*
                      (curr->node[vn3]->loc.x - curr->node[vn2]->loc.x);
            tnode2->loc.y = curr->node[vn2]->loc.y + frac*
                      (curr->node[vn3]->loc.y - curr->node[vn2]->loc.y);
            tnode2->loc.z = curr->node[vn2]->loc.z + frac*
                      (curr->node[vn3]->loc.z - curr->node[vn2]->loc.z);
         }

         // swap the coordinates so that the output is x-y
         if (use_dir == x) {
            tnode1->loc.x = tnode1->loc.y;
            tnode1->loc.y = tnode1->loc.z;
            tnode1->loc.z = thresh;
            tnode2->loc.x = tnode2->loc.y;
            tnode2->loc.y = tnode2->loc.z;
            tnode2->loc.z = thresh;
         } else if (use_dir == y) {
            tnode1->loc.y = tnode1->loc.z;
            tnode1->loc.z = thresh;
            tnode2->loc.y = tnode2->loc.z;
            tnode2->loc.z = thresh;
         }

         // just write the segment!
         if (write_seg_slice(ofp1,output_format,line_width,tnode1,tnode2,&tnorm1,&tnorm2)) {
            fprintf(stderr,"\n");
            fprintf(stderr,"Quitting.\n");
            fclose(ofp1);
            exit(0);
         }

         num_wrote_1++;

      } // end if trim_tri

      if (num_read/DOTPER == (num_read+DPMO)/DOTPER) fprintf(stderr,".");

      curr = curr->next_tri;
   }
   fprintf(stderr,"\n");

   fclose(ofp1);

   fprintf(stderr,"Read %d triangles, wrote %d edges\n",num_read,num_wrote_1);
   /* fprintf(stderr,"Done.\n"); */
   exit(0);
}


/*
 * write_seg_slice - write a segment to the output file
 */
int write_seg_slice(FILE* ofp, char* output_format, double width,
              node_ptr n1, node_ptr n2, VEC *norm1, VEC *norm2) {

   static long int num_nodes = 0;
   static long int num_tangents = 0;
   static int first_time = TRUE;

   // write any necessary headers
   if (first_time) {
      if (strncmp(output_format, "seg", 3) == 0) {
         fprintf (ofp,"d 2\n");
         fprintf (ofp,"gr %g\n",width);
      }
      first_time = FALSE;
   }

   // write the segment in the proper format
   if (strncmp(output_format, "raw", 3) == 0) {
      fprintf(ofp,"%g %g %g %g\n", n1->loc.x, n1->loc.y, n2->loc.x, n2->loc.y);
   } else if (strncmp(output_format, "seg", 3) == 0) {
      // print the nodes, tangents, and segment
      fprintf(ofp,"v %g %g\n", n1->loc.x, n1->loc.y);
      //fprintf(ofp,"vn %g %g\n", norm1->x, norm1->y);
      fprintf(ofp,"v %g %g\n", n2->loc.x, n2->loc.y);
      //fprintf(ofp,"vn %g %g\n", norm2->x, norm2->y);
      num_nodes += 2;
      num_tangents += 2;
      fprintf(ofp,"s %ld %ld\n", num_nodes-1, num_nodes);
      //fprintf(ofp,"s %ld//%ld %ld//%ld\n", num_nodes, num_nodes+1, num_tangents, num_tangents+1);
   //} else if (strncmp(output_format, "svg", 3) == 0) {
   } else {
      fprintf(stderr,"Output file format (%s) unsupported!\n",output_format);
      return(-1);
   }

   return(0);
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

  return;
}


/*
 * This function writes basic usage information to stderr,
 * and then quits. Too bad.
 */
int Usage(char progname[80],int status) {

   /* Usage for rockslice */
   static char **cpp, *help_message[] =
   {
       "where [-options] are one or more of the following:                         ",
       "                                                                           ",
       "   -x val      slice the triangle mesh at the given x value                ",
       "                                                                           ",
       "   -y val      same as for -x, but in y-direction                          ",
       "                                                                           ",
       "   -z val      same as for -x, but in z-direction                          ",
       "                                                                           ",
       "   -s          slice along middle of shortest side                         ",
       "                                                                           ",
       "   -r rad      set radius/width of segments on sliced plane                ",
       "                                                                           ",
       "   -okey       specify output format, key= raw or seg                      ",
       "               default = seg (usable by stickkit)                          ",
       "                                                                           ",
       "   -help       (in place of infile) returns this help information          ",
       " ",
       "The input file can be of .raw, .tin, or .obj format, and the program requires",
       "   the input file to use its valid 3-character filename extension.",
       " ",
       "Options may be abbreviated to an unambiguous length (duh).",
       " ",
       "Program will send output to stdout.",
       "  Example:",
       "    rockslice test.raw -x 0.5 -r 0.01 -oseg > test.seg",
       " ",
       NULL
   };

   fprintf(stderr, "usage:\n  %s infile [-options]\n\n", progname);
   for (cpp = help_message; *cpp; cpp++)
      fprintf(stderr, "%s\n", *cpp);
      fflush(stderr);
   exit(status);
   return(0);
}
