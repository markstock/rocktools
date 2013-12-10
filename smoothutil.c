/*************************************************************
 *
 *  smoothutil.c - Utility subroutines for rocksmooth
 *
 *  Mark J. Stock, mstock@umich.edu
 *
 *
 * rocktools - Tools for creating and manipulating triangular meshes
 * Copyright (C) 1999,2004  Mark J. Stock
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
 ********************************************************** */


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
//#include <malloc.h>
#include "structs.h"

int three_d_laplace(tri_pointer,int);
int three_d_surface_tension(tri_pointer,double);
int compute_normals_2(tri_pointer,int);
int grow_surface_along_normal(tri_pointer,double);
int define_sharp_edges(tri_pointer,double);


/*
 * Run a laplacian-like function over the entire surface to smooth the nodes
 */
int three_d_laplace(tri_pointer tri_head,int num_cycles) {

   int i,j,index;
   //int num = 0;
   VEC sum,tri_norm;
   node_ptr curr_node = node_head;
   //node_ptr test_node;
   tri_pointer test_tri;


   /* Loop this routine a number of times */
   if (num_cycles > 0) fprintf(stderr,"Smoothing surface");
   for (i=0; i<num_cycles; i++) {

      fprintf(stderr,".");
      fflush(stderr);

      /* save all current node locations in the temp_node_record */
      curr_node = node_head;
      while (curr_node) {
         curr_node->temp_loc.x = curr_node->loc.x;
         curr_node->temp_loc.y = curr_node->loc.y;
         curr_node->temp_loc.z = curr_node->loc.z;
         curr_node = curr_node->next_node;
      }

      /* move all nodes a short distance, based on the location
       * of its neighbor nodes */
      curr_node = node_head;
      while (curr_node) {

         /* Any algorithm for this needs to use the same data that finding the
          * surface normal would give. Is it pointless to just take a blind
          * average of adjoining nodes? Probably not. */
         // fprintf(stderr,"   node z=%g",curr_node->loc.z); fflush(stderr);

         /* Until then, let's just use a basic Laplace-like operation */
         sum.x = 0.0;
         sum.y = 0.0;
         sum.z = 0.0;
         for (j=0; j<curr_node->num_adj_nodes; j++) {
            sum.x += curr_node->adj_node[j]->loc.x;
            sum.y += curr_node->adj_node[j]->loc.y;
            sum.z += curr_node->adj_node[j]->loc.z;
         }
         // fprintf(stderr,"   node x=%g",curr_node->loc.x); fflush(stderr);
         curr_node->temp_loc.x = (curr_node->temp_loc.x + 0.1*sum.x/curr_node->num_adj_nodes) / 1.1;
         curr_node->temp_loc.y = (curr_node->temp_loc.y + 0.1*sum.y/curr_node->num_adj_nodes) / 1.1;
         curr_node->temp_loc.z = (curr_node->temp_loc.z + 0.1*sum.z/curr_node->num_adj_nodes) / 1.1;
         // fprintf(stderr," to x=%g\n",curr_node->loc.x); fflush(stderr);

         /* for testing, if this is the last pass, save the delta as a normal */
         if (1==0) {		/* if (i==num_cycles-1) */
            sum.x = sum.x/curr_node->num_adj_nodes - curr_node->temp_loc.x;
            sum.y = sum.y/curr_node->num_adj_nodes - curr_node->temp_loc.y;
            sum.z = sum.z/curr_node->num_adj_nodes - curr_node->temp_loc.z;
            sum = norm(sum);

            /* make sure it is pointing outwards...how?
             * compare the normal from above with the normal of any attached tri */
            test_tri = curr_node->conn_tri[0];
            tri_norm = cross(from(test_tri->node[0]->loc,test_tri->node[1]->loc),from(test_tri->node[0]->loc,test_tri->node[2]->loc));
            if (dot(sum,tri_norm) < 0.0) {
               sum.x = -sum.x;
               sum.y = -sum.y;
               sum.z = -sum.z;
            }

            /* loop over all connected triangles, save normal as norm[] */
            for (j=0; j<curr_node->num_conn; j++) {
               test_tri = curr_node->conn_tri[j];
               test_tri->use_norm = 1;
               index = curr_node->conn_tri_node[j];
               test_tri->norm[index].x = sum.x;
               test_tri->norm[index].y = sum.y;
               test_tri->norm[index].z = sum.z;
            }
         }	/* end retired normals block */

         curr_node = curr_node->next_node;
      }

      /* apply the new locations to the nodes */
      curr_node = node_head;
      while (curr_node) {
         curr_node->loc.x = curr_node->temp_loc.x;
         curr_node->loc.y = curr_node->temp_loc.y;
         curr_node->loc.z = curr_node->temp_loc.z;
         curr_node = curr_node->next_node;
      }
   }
   if (num_cycles > 0) fprintf(stderr,"\n");

   /* return some metric */
   return(0);
}


/*
 * Run a surface-tension-like function over the entire surface to
 * smooth the nodes; ref Tryggvason, JCP 169 (2) p708.
 *
 * apply one second's worth of surface tension with the given coefficient
 */
int three_d_surface_tension(tri_pointer tri_head,double coeff) {

   int i,j;//,index;
   //int num = 0;
   int num_cycles = 10;
   VEC force,tri_norm,oneedge;
   node_ptr curr_node;
   node_ptr startn,endn;
   tri_pointer curr_tri;

   // if coefficient is too high, increase the number of cycles
   // if (coeff > 0.1) num_cycles = 2+(int)(coeff*10.);

   /* Loop this routine a number of times */
   if (num_cycles > 0) fprintf(stderr,"Smoothing surface");
   for (i=0; i<num_cycles; i++) {

      fprintf(stderr,".");
      fflush(stderr);

      /* save all current node locations in the temp_node_record */
      curr_node = node_head;
      while (curr_node) {
         curr_node->temp_loc.x = curr_node->loc.x;
         curr_node->temp_loc.y = curr_node->loc.y;
         curr_node->temp_loc.z = curr_node->loc.z;
         curr_node = curr_node->next_node;
      }

      /* move all nodes a short distance, based on the location
       * of its neighbor nodes */
      curr_tri = tri_head;
      while (curr_tri) {

         // the triangle norm (make it the unit-length negative norm!)
         tri_norm = norm(cross(from(curr_tri->node[0]->loc,curr_tri->node[1]->loc),from(curr_tri->node[2]->loc,curr_tri->node[1]->loc)));
         // tri_norm = cross(from(curr_tri->node[0]->loc,curr_tri->node[1]->loc),from(curr_tri->node[2]->loc,curr_tri->node[1]->loc));

         // for each edge
         for (j=0;j<3;j++) {

            startn = curr_tri->node[j];
            endn = curr_tri->node[(j+1)%3];
            oneedge = from(startn->loc,endn->loc);
            // oneedge = norm(from(startn->loc,endn->loc));

            // force is coefficient times edge cross normal
            force = vscale(0.5*coeff/num_cycles,cross(oneedge,tri_norm));

            // half of the force goes on each end node

            // fprintf(stderr,"   node x=%g",curr_node->loc.x); fflush(stderr);
            startn->temp_loc.x += force.x;
            startn->temp_loc.y += force.y;
            startn->temp_loc.z += force.z;
            endn->temp_loc.x += force.x;
            endn->temp_loc.y += force.y;
            endn->temp_loc.z += force.z;
            // fprintf(stderr," to x=%g\n",curr_node->loc.x); fflush(stderr);
         }

         curr_tri = curr_tri->next_tri;
      }

      /* apply the new locations to the nodes */
      curr_node = node_head;
      while (curr_node) {
         curr_node->loc.x = curr_node->temp_loc.x;
         curr_node->loc.y = curr_node->temp_loc.y;
         curr_node->loc.z = curr_node->temp_loc.z;
         curr_node = curr_node->next_node;
      }
   }
   if (num_cycles > 0) fprintf(stderr,"\n");

   /* return some metric */
   return(0);
}


/*
 * Find the normal vectors of the nodes by averaging the normals
 * of ALL connecting triangles.
 *
 * method 1 uses the mean of the normals of all adjacent elements
 * method 2 uses the area-weighted normals of adj elems
 * method 3 uses the angle-weighted normals of adj elems
 */
int compute_normals_2(tri_pointer tri_head, int method) {

   int index,j;
   double weight;
   VEC sum,tri_norm;
   VEC e1,e2;
   node_ptr curr_node;
   tri_pointer test_tri;

   fprintf(stderr,"Computing normal vectors.");
   fflush(stderr);

   curr_node = node_head;
   while (curr_node) {
      //fprintf(stderr,"\nnode %d\n",curr_node->index); fflush(stderr);
      //fprintf(stderr,"  num_conn %d\n",curr_node->num_conn); fflush(stderr);
      //for (j=0; j<curr_node->num_conn; j++) {
      //   fprintf(stderr,"  tri %d\n",curr_node->conn_tri[j]->index); fflush(stderr);
      //   for (i=0; i<3; i++) fprintf(stderr,"    node %d\n",curr_node->conn_tri[j]->node[i]->index); fflush(stderr);
      //}

      // loop over all connected triangles, find normal as average
      sum.x = 0;
      sum.y = 0;
      sum.z = 0;
      weight = 0.;
      for (j=0; j<curr_node->num_conn; j++) {
         test_tri = curr_node->conn_tri[j];
         //fprintf(stderr,"  tri %d\n",test_tri->index); fflush(stderr);
         //fprintf(stderr,"  nodes %d %d %d\n",test_tri->node[0]->index,test_tri->node[1]->index,test_tri->node[2]->index); fflush(stderr);
         //for (i=0; i<3; i++) fprintf(stderr,"  %d  %g %g %g\n",i,test_tri->node[i]->loc.x,test_tri->node[i]->loc.y,test_tri->node[i]->loc.z); fflush(stderr);
         if (method == 2) {
           // norm can use any nodes
           e1 = from(test_tri->node[0]->loc,test_tri->node[1]->loc);
           e2 = from(test_tri->node[0]->loc,test_tri->node[2]->loc);
           tri_norm = cross(e1,e2);
           // weight based on area of tri
           weight = 0.5*length(tri_norm);
           // must normalize *after* finding area
           tri_norm = norm(tri_norm);
         } else if (method == 3) {
           // norm must use vectors from key node
           index = curr_node->conn_tri_node[j];
           e1 = from(test_tri->node[index]->loc,test_tri->node[(index+1)%3]->loc);
           e2 = from(test_tri->node[index]->loc,test_tri->node[(index+2)%3]->loc);
           tri_norm = cross(e1,e2);
           tri_norm = norm(tri_norm);
           // weight based on angle swept by tri from given node
           weight = acos(dot(e1,e2)/sqrt(lengthsq(e1)*lengthsq(e2)));
           //fprintf(stderr,"  %d wt %g  e1,e2 %g %g  norm %g %g %g\n",j,weight,lengthsq(e1),lengthsq(e2),tri_norm.x,tri_norm.y,tri_norm.z);
         } else {
           // all tris weighted evenly
           weight = 1.;
           // must normalize normal vector nonetheless
           e1 = from(test_tri->node[0]->loc,test_tri->node[1]->loc);
           e2 = from(test_tri->node[0]->loc,test_tri->node[2]->loc);
           tri_norm = cross(e1,e2);
           tri_norm = norm(tri_norm);
         }
         if (!isnan(tri_norm.x)) {
           if (weight == 0.) weight = 1.e-6;
           sum.x += tri_norm.x*weight;
           sum.y += tri_norm.y*weight;
           sum.z += tri_norm.z*weight;
         }
         //wsum += weight;
      }
      //if (sum.x+sum.y+sum.z == 0.) sum.z = 1.;
      //sum.x /= weight;
      //sum.y /= weight;
      //sum.z /= weight;
      sum = norm(sum);
      //fprintf(stderr,"  norm %g %g %g\n",sum.x,sum.y,sum.z);
      if (isnan(sum.x) || fabs(sum.x+sum.y+sum.z) < 1.e-6) {
         //fprintf(stderr,"node %d\n",curr_node->index);
         //fprintf(stderr,"  weight %g\n",weight);
         //fprintf(stderr,"  norm %g %g %g\n",sum.x,sum.y,sum.z);
         //exit(0);
         sum.x = 0.;
         sum.y = 0.;
         sum.z = 1.;
      }

      // loop over all connected triangles, save normal as norm[]
      for (j=0; j<curr_node->num_conn; j++) {
         test_tri = curr_node->conn_tri[j];
         //if (test_tri->use_norm == FALSE) {
            // cannot set the flag here anymore!
            //test_tri->use_norm = TRUE;
            index = curr_node->conn_tri_node[j];
            test_tri->norm[index].x = sum.x;
            test_tri->norm[index].y = sum.y;
            test_tri->norm[index].z = sum.z;
         //}
      }

      curr_node = curr_node->next_node;
   }

   // now, we have all of the normals, set the flags
   test_tri = tri_head;
   while (test_tri) {
     test_tri->use_norm = TRUE;
     test_tri = test_tri->next_tri;
   }

   fprintf(stderr,"\n");

   return(0);
}


/*
 * Find the normal vectors on each tri's nodes by averaging the normals
 * of any connecting triangles that are adjacent and not pointing too
 * sharply away.
 *
 * method 1 uses the mean of the normals of all adjacent elements
 * method 2 uses the area-weighted normals of adj elems
 * method 3 uses the angle-weighted normals of adj elems
 *
 * this method is a copy of compute_normals_2
 */
int compute_normals_3(tri_pointer tri_head, int method, int use_sharp, double sharp_thresh) {

   int index,i,j,icnt,newj;
   double weight,thresh;
   VEC sum,tri_norm,curr_norm;
   VEC e1,e2;
   node_ptr curr_node;
   tri_pointer curr_tri,test_tri,last_tri;

   // ----------------------------------------------------------------
   // first step, make sure all triangles have real area,
   // and not "nan" for a normal!

   fprintf(stderr,"Checking for degenerate tris.");
   fflush(stderr);

   icnt = 0;
   last_tri = NULL;
   curr_tri = tri_head;
   while (curr_tri) {

     // find the normal of this tri
     e1 = from(curr_tri->node[0]->loc,curr_tri->node[1]->loc);
     e2 = from(curr_tri->node[0]->loc,curr_tri->node[2]->loc);
     curr_norm = cross(e1,e2);
     curr_norm = norm(curr_norm);

     // if the normal for this tri is bad, then throw away the triangle!
     // Whoah, that would affect connectivity! Shit.
     // But we have to, or else this bad tri will mess up other tris!
     if (curr_norm.x < 2. && curr_norm.x > -2) {

       // if true, then this is a good tri! go on to the next.
       last_tri = curr_tri;
       curr_tri = curr_tri->next_tri;

     } else {

       // it's a bad tri, remove it CAREFULLY!
       //fprintf(stderr,"\ntri %d IS BAD\n",curr_tri->index);
       //fprintf(stderr,"  norm %g %g %g\n",curr_norm.x,curr_norm.y,curr_norm.z);
       //fprintf(stderr,"  nodes %d %d %d\n",curr_tri->node[0]->index,curr_tri->node[1]->index,curr_tri->node[2]->index);

       // remove the tri from the connectivity list for each node
       for (i=0; i<3; i++) {
         curr_node = curr_tri->node[i];
         //fprintf(stderr,"  node %d has %d connecting tris\n",i,curr_node->num_conn);
         // loop through the list, compacting it as you go
         newj = 0;
         for (j=0; j<curr_node->num_conn; j++) {
           if (curr_node->conn_tri[j] != curr_tri) {
             curr_node->conn_tri[newj] = curr_node->conn_tri[j];
             curr_node->conn_tri_node[newj] = curr_node->conn_tri_node[j];
             newj++;
           }
         }
         // change num_conn as well
         curr_node->num_conn = newj;
         //fprintf(stderr,"      now has %d\n",curr_node->num_conn);
       }

       // now, remove the triangle from the linked list!
       test_tri = curr_tri;
       curr_tri = curr_tri->next_tri;
       free(test_tri);
       last_tri->next_tri = curr_tri;
       // keep last_tri unchanged
     }


     if (icnt/DOTPER == (icnt+DPMO)/DOTPER) fprintf(stderr,".");
     icnt++;
   }
   fprintf(stderr,"\n");

   // ----------------------------------------------------------------
   fprintf(stderr,"Computing normal vectors.");
   fflush(stderr);

   // convert threshhold in degrees to a dot product threshhold
   thresh = cos(sharp_thresh*M_PI/180.);

   icnt = 0;
   curr_tri = tri_head;
   while (curr_tri) {

    // find the normal of this tri
    e1 = from(curr_tri->node[0]->loc,curr_tri->node[1]->loc);
    e2 = from(curr_tri->node[0]->loc,curr_tri->node[2]->loc);
    curr_norm = cross(e1,e2);
    curr_norm = norm(curr_norm);

    // loop over each node in this tri
    for (i=0; i<3; i++) {

      curr_node = curr_tri->node[i];

      //fprintf(stderr,"\nnode %d\n",curr_node->index);
      //fprintf(stderr,"  norm %g %g %g\n",curr_norm.x,curr_norm.y,curr_norm.z);
      //fprintf(stderr,"  node %g %g %g\n",curr_node->loc.x,curr_node->loc.y,curr_node->loc.z);
      //fprintf(stderr,"  num_conn %d\n",curr_node->num_conn); fflush(stderr);
      //for (j=0; j<curr_node->num_conn; j++) {
      //   fprintf(stderr,"  tri %d\n",curr_node->conn_tri[j]->index); fflush(stderr);
         // for (i=0; i<3; i++) fprintf(stderr,"    node %d\n",curr_node->conn_tri[j]->node[i]->index); fflush(stderr);
      //}

      // loop over all connected triangles, find normal as average
      sum.x = 0;
      sum.y = 0;
      sum.z = 0;
      weight = 0.;
      for (j=0; j<curr_node->num_conn; j++) {

         test_tri = curr_node->conn_tri[j];
         if (method == 2) {
           // norm can use any nodes
           e1 = from(test_tri->node[0]->loc,test_tri->node[1]->loc);
           e2 = from(test_tri->node[0]->loc,test_tri->node[2]->loc);
           tri_norm = cross(e1,e2);
           // weight based on area of tri
           weight = 0.5*length(tri_norm);
           // must normalize *after* finding area
           tri_norm = norm(tri_norm);

         } else if (method == 3) {
           // norm must use vectors from key node
           index = curr_node->conn_tri_node[j];
           //fprintf(stderr,"  node %d\n",test_tri->node[index]->index);
           //if (curr_node != test_tri->node[index]) fprintf(stderr,"  whoops");
           e1 = from(test_tri->node[index]->loc,test_tri->node[(index+1)%3]->loc);
           e2 = from(test_tri->node[index]->loc,test_tri->node[(index+2)%3]->loc);
           tri_norm = cross(e1,e2);
           tri_norm = norm(tri_norm);
           // weight based on angle swept by tri from given node
           weight = acos(dot(e1,e2)/sqrt(lengthsq(e1)*lengthsq(e2)));
           //if (weight < 0. || weight > 3.15) {
           //fprintf(stderr,"  %d wt %g norm %g %g %g\n",j,weight,tri_norm.x,tri_norm.y,tri_norm.z);
           //}

         } else {
           // all tris weighted evenly
           weight = 1.;
           // must normalize normal vector nonetheless
           e1 = from(test_tri->node[0]->loc,test_tri->node[1]->loc);
           e2 = from(test_tri->node[0]->loc,test_tri->node[2]->loc);
           tri_norm = cross(e1,e2);
           tri_norm = norm(tri_norm);
         }
         if (use_sharp) {
           // check vs threshhold
           // is the normal of the test tri significantly different from curr_tri?
           if (dot(curr_norm,tri_norm) < thresh) weight = 0.;
           //fprintf(stderr,"     dotp %g thresh %g\n",dot(curr_norm,tri_norm),thresh);
         }
         sum.x += tri_norm.x*weight;
         sum.y += tri_norm.y*weight;
         sum.z += tri_norm.z*weight;
      }
      sum = norm(sum);

      // check for inconclusive normals and assign triangle normal
      if (!(sum.x < 2. && sum.x > -2)) {
        fprintf(stderr,"    normal is bad! %g %g %g\n",sum.x,sum.y,sum.z);
        sum.x = curr_norm.x;
        sum.y = curr_norm.y;
        sum.z = curr_norm.z;
      }
      //if (!(sum.x < 2. && sum.x > -2)) {
      //  fprintf(stderr,"    normal is bad again! %g %g %g\n",sum.x,sum.y,sum.z);
      //  fprintf(stderr,"    %g %g %g\n",curr_tri->node[0]->loc.x,curr_tri->node[0]->loc.y,curr_tri->node[0]->loc.z);
      //  fprintf(stderr,"    %g %g %g\n",curr_tri->node[1]->loc.x,curr_tri->node[1]->loc.y,curr_tri->node[1]->loc.z);
      //  fprintf(stderr,"    %g %g %g\n",curr_tri->node[2]->loc.x,curr_tri->node[2]->loc.y,curr_tri->node[2]->loc.z);
      //}

      // is this normal very different from the tri's normal?
      //if (dot(curr_norm,sum) < thresh) {
      //  fprintf(stderr,"  tri norm %g %g %g\n",curr_norm.x,curr_norm.y,curr_norm.z);
      //  fprintf(stderr,"    node norm %g %g %g\n",sum.x,sum.y,sum.z);
      //  fprintf(stderr,"    dot prod %g\n",dot(curr_norm,sum));
      //}

      // save normal
      curr_tri->norm[i].x = sum.x;
      curr_tri->norm[i].y = sum.y;
      curr_tri->norm[i].z = sum.z;
      //fprintf(stderr,"  final norm %g %g %g\n",sum.x,sum.y,sum.z);

    }
    curr_tri->use_norm = TRUE;
    curr_tri = curr_tri->next_tri;

    if (icnt/DOTPER == (icnt+DPMO)/DOTPER) fprintf(stderr,".");
    icnt++;
   }

   fprintf(stderr,"\n");

   return(0);
}


/*
 * Grow a surface a fixed distance along each nodes' normal
 */
int grow_surface_along_normal(tri_pointer tri_head,double grow_distance) {

   int corner;
   node_ptr curr_node;
   tri_pointer test_tri;

   fprintf(stderr,"Growing boundary...");
   fflush(stderr);

   curr_node = node_head;
   while (curr_node) {

      // find the normal
      test_tri = curr_node->conn_tri[0];
      corner = curr_node->conn_tri_node[0];

      curr_node->loc.x += grow_distance*test_tri->norm[corner].x;
      curr_node->loc.y += grow_distance*test_tri->norm[corner].y;
      curr_node->loc.z += grow_distance*test_tri->norm[corner].z;

      curr_node = curr_node->next_node;
   }

   fprintf(stderr,"\n");

   return(0);
}
