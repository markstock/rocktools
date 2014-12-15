/*************************************************************
 *
 *  detailutil.c - Utility subroutines for use primarily 
 *	with rockdetail
 *
 *  Mark J. Stock, mstock@umich.edu
 *
 *
 * rocktools - Tools for creating and manipulating triangular meshes
 * Copyright (C) 1999,2002-4,6,14  Mark J. Stock
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


tri_pointer split_tri(int,tri_pointer);
tri_pointer split_tri_hex(int,tri_pointer);
tri_pointer split_tri_5(int,tri_pointer);
int make_sphere(tri_pointer);
int find_adjacent_child(tri_pointer,tri_pointer,node_ptr,node_ptr);
node_ptr create_midpoint(node_ptr,node_ptr);
node_ptr create_midpoint_2(int,node_ptr,node_ptr);
node_ptr create_midpoint_3(int,node_ptr,node_ptr,node_ptr);
node_ptr create_midpoint_4(int,node_ptr,node_ptr,node_ptr,node_ptr);
node_ptr create_midpoint_5(int,node_ptr,node_ptr,node_ptr,node_ptr);
void move_existing_node_5( int, node_ptr);
void perturb_node_5 (VEC*, int, VEC);
node_ptr create_midpoint_spline(node_ptr,node_ptr,int*);
node_ptr create_center_point(int,node_ptr,node_ptr,node_ptr);

// from smoothutil.c
extern void compute_normals_2 (tri_pointer,int);

extern double normal_shake;
extern double normal_exponent;
extern double normal_bias;
extern double base_shake;
extern double base_exponent;
extern int use_spline;
extern int use_gaussian_random;
extern int clamp_edges;
extern int use_thresh;
extern double area_thresh;
extern int use_dist;
extern double distance_thresh;
extern VEC viewp;
extern int force_sphere;
extern double sphere_rad;
extern double find_tri_dist(tri_pointer,VEC);

/*
 * split_tri_hex takes a linked list of triangles and splits each
 * triangle into 6/4/2/0 new triangles. It adds 1 new node to
 * the node list, and places it close to, but not exactly on, the 
 * the center of the original triangle.
 *
 * This routine mimics the hexagonal subdivision scheme originally
 * proposed by Mandelbrot, and should reduce the effect of edging
 * over successive generations.
 */
tri_pointer split_tri_hex(int depth,tri_pointer tri_head) {

   int i,j,k;
   int adj_side = -1;
   int local_side[6];
   //int tindex = 0;
   //int nindex = 100;
   double temp_area;
   tri_pointer this_tri;
   tri_pointer new_tri_head = NULL;
   tri_pointer new_tri[2];
   tri_pointer local_tri[6];
   node_ptr new_node,new_node2;

   fprintf(stderr,"Method 2, depth = %d\n",depth); fflush(stderr);

   /* for each triangle in the old list */
   this_tri = tri_head;
   while (this_tri) {

      // fprintf(stderr,"Checking tri %d\n",this_tri->index);
      // fprintf(stdout,"  this tri has nodes %d %d %d\n",this_tri->node[0]->index,this_tri->node[1]->index,this_tri->node[2]->index);

      // if the tri is too small, do not split it, just copy it to the new list
      // For first take, split all tris
      if (!this_tri->splittable) {
         new_tri[0] = alloc_new_tri();
         new_tri[0]->splittable = FALSE;
         new_tri[0]->node[0] = this_tri->node[0];
         new_tri[0]->node[1] = this_tri->node[1];
         new_tri[0]->node[2] = this_tri->node[2];
         new_tri[0]->adjacent[0] = this_tri->adjacent[0];
         new_tri[0]->adjacent[1] = this_tri->adjacent[1];
         new_tri[0]->adjacent[2] = this_tri->adjacent[2];
         new_tri[0]->next_tri = new_tri_head;
         new_tri_head = new_tri[0];
         this_tri = this_tri->next_tri;
         continue;
      }

      // choose the one new node location
      new_node = create_center_point(depth,this_tri->node[0],this_tri->node[1],this_tri->node[2]);

      // for each side, check to see if two triangles were already made
      for (i=0; i<3; i++) {
         j = (i+1)%3;

         // fprintf(stdout,"adjacent tri %d is %d\n",i,this_tri->adjacent[i]->index);
         // fprintf(stdout,"  has nodes %d %d %d\n",this_tri->adjacent[i]->node[0]->index,this_tri->adjacent[i]->node[1]->index,this_tri->adjacent[i]->node[2]->index);
         // fflush(stdout);

         // has side i been split?
         if (this_tri->midpoint[i]) {

            // if so, tell this side's two tris what their new node pointer is

            // which of the neighbor's children need this info?
            // find the index of the shared side from the adjacent tri's entry
            for (k=0; k<3; k++) {
               if (this_tri->node[i] == this_tri->adjacent[i]->node[k]) {

                  if (this_tri->node[j] == this_tri->adjacent[i]->node[(k+2)%3]) {
                     // both tris are oriented the same
                     // fprintf(stdout,"  side has been split, k=%d, tris are oriented\n",k);
                     // fflush(stdout);
                     adj_side = mod(k+2,3);
                     local_tri[i*2] = this_tri->adjacent[i]->adjacent[adj_side]->next_tri;
                     local_side[i*2] = 0;
                     local_tri[i*2+1] = this_tri->adjacent[i]->adjacent[adj_side];
                     local_side[i*2+1] = 2;
                     break;
                  } else {
                     // tris are oriented opposite
                     // fprintf(stdout,"  side has been split, k=%d, tris are opposite\n",k);
                     // fflush(stdout);
                     adj_side = k;
                     local_tri[i*2] = this_tri->adjacent[i]->adjacent[adj_side];
                     local_side[i*2] = 2;
                     local_tri[i*2+1] = this_tri->adjacent[i]->adjacent[adj_side]->next_tri;
                     local_side[i*2+1] = 0;
                     break;
                  }
               }
            }
            // fprintf(stdout,"  index of the shared side %d\n",adj_side);
            // fflush(stdout);

            // set the existing side tris' empty node to new_node
            this_tri->adjacent[i]->adjacent[adj_side]->node[0] = new_node;
            this_tri->adjacent[i]->adjacent[adj_side]->next_tri->node[0] = new_node;

            // let the parent point to the first of the two children
            //   on this side, regardless of the orientation
            this_tri->adjacent[i] = this_tri->adjacent[i]->adjacent[adj_side];

            // calculate triangle's area, if below either threshhold, set it up so it will not split
            for (k=0; k<2; k++) local_tri[i*2+k]->splittable = TRUE;
            if (use_thresh || use_dist) {
               for (k=0; k<2; k++) {
                  temp_area = find_area(local_tri[i*2+k]);
                  if (use_thresh) {
                     if (temp_area < area_thresh) {
                        // this is a flag to the splitter, do not split further
                        local_tri[i*2+k]->splittable = FALSE;
                     }
                  }
                  if (use_dist) {
                     // if (sqrt(temp_area)/length(from(viewp,local_tri[i*2+k]->node[0]->loc)) < distance_thresh) {
                     if (sqrt(temp_area)/find_tri_dist(local_tri[i*2+k],viewp) < distance_thresh) {
                        // this is a flag to the splitter, do not split further
                        local_tri[i*2+k]->splittable = FALSE;
                     }
                  }
               }
            }


         } else {	// if side i has not been split

            // if not, create two new elements, fill in as much data as possible
            // fprintf(stdout,"  side has not been split\n");
            // fflush(stdout);

            // initialize two new elements
            new_tri[0] = alloc_new_tri();
            new_tri[1] = alloc_new_tri();

            // define two new elements
            new_tri[0]->node[0] = NULL;
            new_tri[0]->node[1] = new_node;
            new_tri[0]->node[2] = this_tri->node[i];
            new_tri[0]->adjacent[0] = new_tri[1];
            new_tri[0]->adjacent[1] = NULL;
            new_tri[0]->adjacent[2] = NULL;
            new_tri[1]->node[0] = NULL;
            new_tri[1]->node[1] = this_tri->node[j];
            new_tri[1]->node[2] = new_node;
            new_tri[1]->adjacent[0] = NULL;
            new_tri[1]->adjacent[1] = NULL;
            new_tri[1]->adjacent[2] = new_tri[0];

            // new tri's midpoints will not be set this recursion level, set to NULL
            for (k=0; k<2; k++) for (j=0; j<3; j++) new_tri[k]->midpoint[j] = NULL;
            j = (i+1)%3;

            local_tri[i*2] = new_tri[0];
            local_side[i*2] = 1;
            local_tri[i*2+1] = new_tri[1];
            local_side[i*2+1] = 1;

            if (this_tri->adjacent[i]) {

               // if there exists a neighboring parent, tell it that we've
               //    got a new node over here

               // unless, of course, that element is to never be split again
               if (!this_tri->adjacent[i]->splittable) {

                  // in that case, just split the edge at the midpoint
                  //    and never deal with it again.
                  j = (i+1)%3;
                  new_node2 = create_midpoint(this_tri->node[i],this_tri->node[j]);
                  new_tri[0]->node[0] = new_node2;
                  new_tri[1]->node[0] = new_node2;

                  // now, must check the triangles created here for smallness
                  // calculate triangle's area, if below either threshhold, set it up so it will not split
                  for (k=0; k<2; k++) new_tri[k]->splittable = TRUE;
                  if (use_thresh || use_dist) {
                     for (k=0; k<2; k++) {
                        temp_area = find_area(new_tri[k]);
                        if (use_thresh) {
                           if (temp_area < area_thresh) {
                              // this is a flag to the splitter, do not split further
                              new_tri[k]->splittable = FALSE;
                           }
                        }
                        if (use_dist) {
                           // if (sqrt(temp_area)/length(from(viewp,new_tri[k]->node[0]->loc)) < distance_thresh) {
                           if (sqrt(temp_area)/find_tri_dist(new_tri[k],viewp) < distance_thresh) {
                              // this is a flag to the splitter, do not split further
                              new_tri[k]->splittable = FALSE;
                           }
                        }
                     }
                  }


               } else {
                  // the neighboring parent exists and isn't too small
            
                  // find the index of the shared side from the adjacent tri's entry
                  for (k=0; k<3; k++) {
                     // fprintf(stdout,"    a neighbor parent exists %d\n",k);
                     // fflush(stdout);
                     if (this_tri->node[i] == this_tri->adjacent[i]->node[k]) {
   
                        if (this_tri->node[j] == this_tri->adjacent[i]->node[(k+2)%3]) {
                           // both tris are oriented the same
                           adj_side = (k+2)%3;
                           break;
                        } else {
                           // tris are oriented opposite
                           adj_side = k;
                           break;
                        }
                     }
                  }
               }
               // fprintf(stdout,"  index of the shared side %d\n",adj_side);
               // fflush(stdout);
               this_tri->adjacent[i]->midpoint[adj_side] = new_node;

            } else {

               // if there doesn't even *exist* an adjacent triangle, set a new edge point
               // ALWAYS use the true midpoint, because there may be
               //   a neighbor triangle that is too small to split
               //   Well, not any more, with the if above

               // use the better method for determining a new midpoint
               k = (j+1)%3;

               if (clamp_edges) {
                  // clamp the edge, use the mathematical midpoint
                  new_node2 = create_midpoint(this_tri->node[i],this_tri->node[j]);
               } else {
                  // let the edge get jagged, use the current tri's 3 points
                  new_node2 = create_midpoint_3(depth,this_tri->node[i],
                                                   this_tri->node[j],
                                                   this_tri->node[k]);
               }

               // and set the appropriate node pointers
               new_tri[0]->node[0] = new_node2;
               new_tri[1]->node[0] = new_node2;
               // new_tri->adjacents already point to NULL

               // calculate triangle's area, if below either threshhold, set it up so it will not split
               for (k=0; k<2; k++) new_tri[k]->splittable = TRUE;
               if (use_thresh || use_dist) {
                  for (k=0; k<2; k++) {
                     temp_area = find_area(new_tri[k]);
                     if (use_thresh) {
                        if (temp_area < area_thresh) {
                           // this is a flag to the splitter, do not split further
                           new_tri[k]->splittable = FALSE;
                        }
                     }
                     if (use_dist) {
                        // if (sqrt(temp_area)/length(from(viewp,new_tri[k]->node[0]->loc)) < distance_thresh) {
                        if (sqrt(temp_area)/find_tri_dist(new_tri[k],viewp) < distance_thresh) {
                           // this is a flag to the splitter, do not split further
                           new_tri[k]->splittable = FALSE;
                        }
                     }
                  }
               }

            }

            // let the parent point to the first of the two children
            this_tri->adjacent[i] = new_tri[0];

            // add the 2 new tris to the new list
            new_tri[1]->next_tri = new_tri_head;
            new_tri[0]->next_tri = new_tri[1];
            new_tri_head = new_tri[0];

         }   // end if side has been split

      }   // end loop for i=0 to 2


      // now that each of the three sides has been created/updated, we need to
      // update the adjacent pointers within the parent tri
      /*
      for (i=0; i<6; i++) fprintf(stdout,"    local_side[%d] = %d\n",i,local_side[i]);
      for (i=0; i<6; i++) {
         if (local_tri[i]) {
            fprintf(stdout,"    local_tri[%d] exists\n",i);
            fflush(stdout);
         }
      }
      */

      local_tri[0]->adjacent[local_side[0]] = local_tri[5];
      local_tri[5]->adjacent[local_side[5]] = local_tri[0];

      local_tri[2]->adjacent[local_side[2]] = local_tri[1];
      local_tri[1]->adjacent[local_side[1]] = local_tri[2];

      local_tri[4]->adjacent[local_side[4]] = local_tri[3];
      local_tri[3]->adjacent[local_side[3]] = local_tri[4];

      // for my own amusement, see if all of this parent->adjacent
      //   i.e. children, have a next_tri
      /*
      for (i=0; i<3; i++) {
         if (this_tri->adjacent[i]->next_tri) {
            fprintf(stdout,"    this_tri->adjacent[%d]->next_tri exists\n",i);
            fflush(stdout);
         } else {
            fprintf(stdout,"    this_tri->adjacent[%d]->next_tri DOES NOT exist\n",i);
            fflush(stdout);
         }
      }
      */


      // we are now done with this parent tri, we may reference it later, though
      this_tri = this_tri->next_tri;
   }

   // now, call the sphericalizing routine
   if (force_sphere) make_sphere(new_tri_head);

   /* replace the old list with the new list */
   return new_tri_head;
}


/*
 * split_tri takes a linked list of triangles and splits each
 * triangle into 4 new triangles. It adds 3 new nodes to the node
 * list, and places them close to, but not exactly on, the 
 * the midpoint between any two existing nodes.
 */
tri_pointer split_tri(int depth,tri_pointer tri_head) {

   int h,i,j,k;
   int adj_side = -1;
   int far_corner = -1;
   int num_new_tri;
   int tri_cnt = 0;
   int has_adjacent_been_split[3];
   double temp_area;
   tri_pointer this_tri;
   tri_pointer new_tri_head = NULL;
   tri_pointer new_tri[4];
   node_ptr new_node[3];

   fprintf(stderr,"Method 1, depth = %d\n",depth);
   j = -1; k = -1;

   // for each triangle in the old list
   this_tri = tri_head;
   while (this_tri) {

      //fprintf(stderr,"Checking tri...\n");

      // if the tri is too small, do not split it, just copy it to the new list
      // new logic: if the tri is flagged to not split, split it minimally to
      //    account for neighboring splittable triangles!
      if (!this_tri->splittable) {

         //fprintf(stderr,"Unsplittable tri %d\n",this_tri->index); fflush(stderr);

         // first, see how many neighbors will be split
         num_new_tri = 1;
         for (i=0; i<3; i++) new_node[i] = NULL;
// NOTE: to make it work, i<3 on next line!
         for (i=0; i<3; i++) {

            // check to see if a midpoint split already exists, use it if it does
            if (this_tri->midpoint[i]) {

               // that edge has already been split, use the point
               has_adjacent_been_split[i] = TRUE;
               new_node[i] = this_tri->midpoint[i];
               num_new_tri++;

               // if (this_tri->adjacent[i]->splittable) {
                  // fprintf(stderr,"     corresponding to splittable tri %d\n",this_tri->adjacent[i]->index); fflush(stderr);
               // } else {
                  // fprintf(stderr,"     corresponding to unsplittable tri %d\n",this_tri->adjacent[i]->index); fflush(stderr);
               // }

            } else {

               // if it doesn't, see if it eventually will
               has_adjacent_been_split[i] = FALSE;
               j = mod(i+1,3);
               k = mod(j+1,3);

               // fprintf(stderr,"   creating midpoint for side %d\n",i); fflush(stderr);
               if (!this_tri->adjacent[i]) {

                  // there is no adjacent triangle, don't split this edge

               } else if (this_tri->adjacent[i]->splittable) {

                  // fprintf(stderr,"     corresponding to splittable tri %d\n",this_tri->adjacent[i]->index); fflush(stderr);

                  // there is an adjacent triangle that can be split
                  // find the index of the adjacent tri's farthest node
                  for (h=0; h<3; h++) {
                     if (this_tri->node[i] == this_tri->adjacent[i]->node[h]) {
                        far_corner = mod(h+1,3);
                        break;
                     }
                  }
                  new_node[i] = create_midpoint_5(depth,this_tri->node[i],
                                                  this_tri->node[j],
                                                  this_tri->node[k],
                                                  this_tri->adjacent[i]->node[far_corner]);

                  // find the index of the shared side from the adjacent tri's entry
                  for (j=0; j<3; j++) {
                     if (this_tri->node[i] == this_tri->adjacent[i]->node[j]) {
                        adj_side = mod(j+2,3);
                        break;
                     }
                  }
                  j = mod(adj_side+1,3);
                  this_tri->adjacent[i]->midpoint[adj_side] = new_node[i];
                  // fprintf(stderr,"      and told adjacent tri %d that its side %d has a midpoint already\n",i,adj_side);
                  num_new_tri++;

               } else {

                  // fprintf(stderr,"     corresponding to unsplittable tri %d\n",this_tri->adjacent[i]->index); fflush(stderr);

                  // there is an adjacent triangle, but it's flagged to not split,
                  //    so don't split this side at all

               }
            }
         }

         // now, based on how many tris we have, create them
         // fprintf(stderr,"  splitting into %d tris\n",num_new_tri); fflush(stderr);

         // create the four new triangles
         for (i=0; i<num_new_tri; i++) {
            new_tri[i] = alloc_new_tri();
            new_tri[i]->splittable = FALSE;
         }

         if (num_new_tri == 4) {

            // treat it just like a regular split element, but flag all 4
            //    children as non-splittable

            // no change to old method
            // by default, set all external adjacent pointers to NULL, fix later
            new_tri[0]->node[0] = this_tri->node[0];
            new_tri[0]->node[1] = new_node[0];
            new_tri[0]->node[2] = new_node[2];
            new_tri[0]->adjacent[0] = NULL;
            new_tri[0]->adjacent[1] = new_tri[3];
            new_tri[0]->adjacent[2] = NULL;
            new_tri[0]->index = tri_cnt++;
            new_tri[1]->node[0] = new_node[0];
            new_tri[1]->node[1] = this_tri->node[1];
            new_tri[1]->node[2] = new_node[1];
            new_tri[1]->adjacent[0] = NULL;
            new_tri[1]->adjacent[1] = NULL;
            new_tri[1]->adjacent[2] = new_tri[3];
            new_tri[1]->index = tri_cnt++;
            new_tri[2]->node[0] = new_node[2];
            new_tri[2]->node[1] = new_node[1];
            new_tri[2]->node[2] = this_tri->node[2];
            new_tri[2]->adjacent[0] = new_tri[3];
            new_tri[2]->adjacent[1] = NULL;
            new_tri[2]->adjacent[2] = NULL;
            new_tri[2]->index = tri_cnt++;
            new_tri[3]->node[0] = new_node[1];
            new_tri[3]->node[1] = new_node[2];
            new_tri[3]->node[2] = new_node[0];
            new_tri[3]->adjacent[0] = new_tri[2];
            new_tri[3]->adjacent[1] = new_tri[0];
            new_tri[3]->adjacent[2] = new_tri[1];
            new_tri[3]->index = tri_cnt++;

            // If parent had no adjacent tri, children on that edge have none, either.
            // BUT, if the parent DID have an adjacent tri, children will also.
            // So, if the parent had an adjacent tri, and it has been split, find the
            // children's adjacent tris among the parent's adjacent's children
            if (this_tri->adjacent[0] && has_adjacent_been_split[0]) {

               // then use the child tri pointed to by this_tri->adjacent[0]->adjacent[0]
               // to begin searching for the tris adjacent to the 2 new ones on this side

               if (!find_adjacent_child(this_tri->adjacent[0]->adjacent[0],
                     new_tri[0],this_tri->node[0],new_node[0]) )
                  fprintf(stderr,"Could not find adjacent child 1.\n");
               if (!find_adjacent_child(this_tri->adjacent[0]->adjacent[0],
                     new_tri[1],new_node[0],this_tri->node[1]) )
                  fprintf(stderr,"Could not find adjacent child 2.\n");
            }

            // do the same for parent's side 1
            if (this_tri->adjacent[1] && has_adjacent_been_split[1]) {
               if (!find_adjacent_child(this_tri->adjacent[1]->adjacent[0],
                     new_tri[1],this_tri->node[1],new_node[1]) )
                  fprintf(stderr,"Could not find adjacent child 3.\n");
               if (!find_adjacent_child(this_tri->adjacent[1]->adjacent[0],
                     new_tri[2],new_node[1],this_tri->node[2]) )
                  fprintf(stderr,"Could not find adjacent child 4.\n");
            }

            // and for parent's side 2
            if (this_tri->adjacent[2] && has_adjacent_been_split[2]) {
               if (!find_adjacent_child(this_tri->adjacent[2]->adjacent[0],
                     new_tri[2],this_tri->node[2],new_node[2]) )
                  fprintf(stderr,"Could not find adjacent child 5.\n");
               if (!find_adjacent_child(this_tri->adjacent[2]->adjacent[0],
                     new_tri[0],new_node[2],this_tri->node[0]) )
                  fprintf(stderr,"Could not find adjacent child 6.\n");
            }

            // new tri's midpoints will not be set this recursion level, set to NULL
            for (i=0; i<4; i++) for (j=0; j<3; j++) new_tri[i]->midpoint[j] = NULL;

            // Important: set this now-split parent triangle's first adjacent pointer to the
            // first of the 4 new triangles created, this information will be useful later
            this_tri->adjacent[0] = new_tri[0];

         } else if (num_new_tri == 3) {

            // first, which edge doesn't have a new_node?
            for (h=0; h<3; h++) if (!new_node[h]) {
               i = (h+1)%3;
               j = (h+2)%3;
               k = h;
               break;
            }

            // geometry of the three triangles
            new_tri[0]->node[0] = this_tri->node[i];
            new_tri[0]->node[1] = new_node[i];
            new_tri[0]->node[2] = this_tri->node[k];
            new_tri[0]->adjacent[0] = NULL;
            new_tri[0]->adjacent[1] = new_tri[2];
            new_tri[0]->adjacent[2] = NULL;
            new_tri[0]->index = tri_cnt++;
            new_tri[1]->node[0] = new_node[i];
            new_tri[1]->node[1] = this_tri->node[j];
            new_tri[1]->node[2] = new_node[j];
            new_tri[1]->adjacent[0] = NULL;
            new_tri[1]->adjacent[1] = NULL;
            new_tri[1]->adjacent[2] = new_tri[2];
            new_tri[1]->index = tri_cnt++;
            new_tri[2]->node[0] = new_node[i];
            new_tri[2]->node[1] = new_node[j];
            new_tri[2]->node[2] = this_tri->node[k];
            new_tri[2]->adjacent[0] = new_tri[1];
            new_tri[2]->adjacent[1] = NULL;
            new_tri[2]->adjacent[2] = new_tri[0];
            new_tri[2]->index = tri_cnt++;

            // now, set those adjacent pointers
            if (this_tri->adjacent[i] && has_adjacent_been_split[i]) {
               if (!find_adjacent_child(this_tri->adjacent[i]->adjacent[0],
                     new_tri[0],this_tri->node[i],new_node[i]) )
                  fprintf(stderr,"Could not find adjacent child 1.\n");
               if (!find_adjacent_child(this_tri->adjacent[i]->adjacent[0],
                     new_tri[1],new_node[i],this_tri->node[j]) )
                  fprintf(stderr,"Could not find adjacent child 2.\n");
            }
            if (this_tri->adjacent[j] && has_adjacent_been_split[j]) {
               if (!find_adjacent_child(this_tri->adjacent[j]->adjacent[0],
                     new_tri[1],this_tri->node[j],new_node[j]) )
                  fprintf(stderr,"Could not find adjacent child 3.\n");
               if (!find_adjacent_child(this_tri->adjacent[j]->adjacent[0],
                     new_tri[2],new_node[j],this_tri->node[k]) )
                  fprintf(stderr,"Could not find adjacent child 4.\n");
            }
            if (this_tri->adjacent[k] && has_adjacent_been_split[k]) {
               fprintf(stderr,"Side k should not have been split!\n");
            }

            // set that other stuff
            for (i=0; i<3; i++) for (j=0; j<3; j++) new_tri[i]->midpoint[j] = NULL;
            this_tri->adjacent[0] = new_tri[0];

         } else if (num_new_tri == 2) {

            // first, which edge has the new_node?
            for (h=0; h<3; h++) if (new_node[h]) {
               i = h;
               j = (h+1)%3;
               k = (h+2)%3;
               break;
            }

            // geometry of the two triangles
            new_tri[0]->node[0] = this_tri->node[i];
            new_tri[0]->node[1] = new_node[i];
            new_tri[0]->node[2] = this_tri->node[k];
            new_tri[0]->adjacent[0] = NULL;
            new_tri[0]->adjacent[1] = new_tri[1];
            new_tri[0]->adjacent[2] = NULL;
            new_tri[0]->index = tri_cnt++;
            new_tri[1]->node[0] = new_node[i];
            new_tri[1]->node[1] = this_tri->node[j];
            new_tri[1]->node[2] = this_tri->node[k];
            new_tri[1]->adjacent[0] = NULL;
            new_tri[1]->adjacent[1] = NULL;
            new_tri[1]->adjacent[2] = new_tri[0];
            new_tri[1]->index = tri_cnt++;

            // now, set those adjacent pointers
            if (this_tri->adjacent[i] && has_adjacent_been_split[i]) {
               if (!find_adjacent_child(this_tri->adjacent[i]->adjacent[0],
                     new_tri[0],this_tri->node[i],new_node[i]) )
                  fprintf(stderr,"Could not find adjacent child 1.\n");
               if (!find_adjacent_child(this_tri->adjacent[i]->adjacent[0],
                     new_tri[1],new_node[i],this_tri->node[j]) )
                  fprintf(stderr,"Could not find adjacent child 2.\n");
            }
            if (this_tri->adjacent[j] && has_adjacent_been_split[j]) {
               fprintf(stderr,"Side j should not have been split!\n");
            }
            if (this_tri->adjacent[k] && has_adjacent_been_split[k]) {
               fprintf(stderr,"Side k should not have been split!\n");
            }


            // set that other stuff
            for (i=0; i<2; i++) for (j=0; j<3; j++) new_tri[i]->midpoint[j] = NULL;
            this_tri->adjacent[0] = new_tri[0];

         } else {	// if num_new_tri==1

            // then just pass this tri on to the next stage
            // new_tri[0] = alloc_new_tri();
            // new_tri[0]->splittable = FALSE;
            new_tri[0]->node[0] = this_tri->node[0];
            new_tri[0]->node[1] = this_tri->node[1];
            new_tri[0]->node[2] = this_tri->node[2];
            new_tri[0]->adjacent[0] = NULL;
            new_tri[0]->adjacent[1] = NULL;
            new_tri[0]->adjacent[2] = NULL;
            new_tri[0]->index = tri_cnt++;
            // new_tri[0]->next_tri = new_tri_head;
            // new_tri_head = new_tri[0];
            // this_tri = this_tri->next_tri;
            // continue;

         }

         // add all of the tris to the list
         for (i=num_new_tri-1; i>-1; i--) {
            new_tri[i]->next_tri = new_tri_head;
            new_tri_head = new_tri[i];
         }

         // and jump to the next parent tri
         this_tri = this_tri->next_tri;
         continue;
      }

      //--------------------------------------------------------------------------
      // if this is a normal, splittable triangle, create 4 child triangles

      //fprintf(stderr,"Splittable tri %d\n",this_tri->index); fflush(stderr);

      // choose the three edge split points
      for (i=0; i<3; i++) new_node[i] = NULL;
      for (i=0; i<3; i++) {

         // check to see if a midpoint split already exists, use it if it does
         if (this_tri->midpoint[i]) {

            has_adjacent_been_split[i] = TRUE;
            //fprintf(stderr,"   already have midpoint on side %d\n",i); fflush(stderr);
            new_node[i] = this_tri->midpoint[i];

            // if (this_tri->adjacent[i]->splittable) {
               // fprintf(stderr,"     corresponding to splittable tri %d\n",this_tri->adjacent[i]->index); fflush(stderr);
            // } else {
               // fprintf(stderr,"     corresponding to unsplittable tri %d\n",this_tri->adjacent[i]->index); fflush(stderr);
            // }

         } else {

            /* if it doesn't, choose one and send the information to the adjacent tri, if there is one */
            has_adjacent_been_split[i] = FALSE;
            j = mod(i+1,3);

            //fprintf(stderr,"   creating midpoint for side %d\n",i); fflush(stderr);
            /* use the better method for determining a new midpoint */
            k = mod(j+1,3);

            if (!this_tri->adjacent[i]) {

               //fprintf(stderr,"     there is no adjacent tri %d\n"); fflush(stderr);

               /* there is no adjacent triangle */
               if (clamp_edges) {
                  /* clamp the edge, use the mathematical midpoint */
                  new_node[i] = create_midpoint(this_tri->node[i],this_tri->node[j]);
               } else {
                  /* let the edge get jagged, use the current tri's 3 points */
                  new_node[i] = create_midpoint_3(depth,this_tri->node[i],
                                                  this_tri->node[j],
                                                  this_tri->node[k]);
               }

            } else if (this_tri->adjacent[i]->splittable) {

               //fprintf(stderr,"     corresponding to splittable tri %d\n",this_tri->adjacent[i]->index); fflush(stderr);

               // there is an adjacent triangle that can be split
               // find the index of the adjacent tri's farthest node
               for (h=0; h<3; h++) {
                  if (this_tri->node[i] == this_tri->adjacent[i]->node[h]) {
                     far_corner = mod(h+1,3);
                     break;
                  }
               }
               new_node[i] = create_midpoint_5(depth,this_tri->node[i],
                                               this_tri->node[j],
                                               this_tri->node[k],
                                               this_tri->adjacent[i]->node[far_corner]);

               // find the index of the shared side from the adjacent tri's entry
               for (j=0; j<3; j++) {
                  if (this_tri->node[i] == this_tri->adjacent[i]->node[j]) {
                     adj_side = mod(j+2,3);
                     break;
                  }
               }
               j = mod(adj_side+1,3);
               this_tri->adjacent[i]->midpoint[adj_side] = new_node[i];
               //fprintf(stderr,"      and told adjacent tri %d that its side %d has a midpoint already\n",i,adj_side); fflush(stderr);

            } else {

               //fprintf(stderr,"     corresponding to unsplittable tri %d\n",this_tri->adjacent[i]->index); fflush(stderr);

               // there is an adjacent triangle, but it's set to not split;
               // clamp the edge, use the mathematical midpoint
               // new_node[i] = create_midpoint(this_tri->node[i],this_tri->node[j]);

               // new method: split that edge normally, anyways! AND tell the
               //    adjacent tri that we split

               // find the index of the adjacent tri's farthest node
               for (h=0; h<3; h++) {
                  if (this_tri->node[i] == this_tri->adjacent[i]->node[h]) {
                     far_corner = mod(h+1,3);
                     break;
                  }
               }
               new_node[i] = create_midpoint_5(depth,this_tri->node[i],
                                               this_tri->node[j],
                                               this_tri->node[k],
                                               this_tri->adjacent[i]->node[far_corner]);
               // find the index of the shared side from the adjacent tri's entry
               for (j=0; j<3; j++) {
                  if (this_tri->node[i] == this_tri->adjacent[i]->node[j]) {
                     adj_side = mod(j+2,3);
                     break;
                  }
               }
               j = mod(adj_side+1,3);
               this_tri->adjacent[i]->midpoint[adj_side] = new_node[i];
            }

         }
      }

      // create the four new triangles
      for (i=0; i<4; i++) new_tri[i] = alloc_new_tri();

      // no change to old method
      // by default, set all external adjacent pointers to NULL, fix later
      new_tri[0]->node[0] = this_tri->node[0];
      new_tri[0]->node[1] = new_node[0];
      new_tri[0]->node[2] = new_node[2];
      new_tri[0]->adjacent[0] = NULL;
      new_tri[0]->adjacent[1] = new_tri[3];
      new_tri[0]->adjacent[2] = NULL;
      new_tri[0]->index = tri_cnt++;
      new_tri[1]->node[0] = new_node[0];
      new_tri[1]->node[1] = this_tri->node[1];
      new_tri[1]->node[2] = new_node[1];
      new_tri[1]->adjacent[0] = NULL;
      new_tri[1]->adjacent[1] = NULL;
      new_tri[1]->adjacent[2] = new_tri[3];
      new_tri[1]->index = tri_cnt++;
      new_tri[2]->node[0] = new_node[2];
      new_tri[2]->node[1] = new_node[1];
      new_tri[2]->node[2] = this_tri->node[2];
      new_tri[2]->adjacent[0] = new_tri[3];
      new_tri[2]->adjacent[1] = NULL;
      new_tri[2]->adjacent[2] = NULL;
      new_tri[2]->index = tri_cnt++;
      new_tri[3]->node[0] = new_node[1];
      new_tri[3]->node[1] = new_node[2];
      new_tri[3]->node[2] = new_node[0];
      new_tri[3]->adjacent[0] = new_tri[2];
      new_tri[3]->adjacent[1] = new_tri[0];
      new_tri[3]->adjacent[2] = new_tri[1];
      new_tri[3]->index = tri_cnt++;

      /* If parent had no adjacent tri, children on that edge have none, either.
       * BUT, if the parent DID have an adjacent tri, children will also.
       * So, if the parent had an adjacent tri, and it has been split, find the children's
       * adjacent tris among the parent's adjacent's children */
      /* but what about the actual node locations? NO, this is just setup for the next step */
      if (this_tri->adjacent[0] && has_adjacent_been_split[0]) {

         /* then use the child tri pointed to by this_tri->adjacent[0]->adjacent[0]
          * to begin searching for the triangles adjacent to the 2 new ones on this side */

         if (!find_adjacent_child(this_tri->adjacent[0]->adjacent[0],new_tri[0],this_tri->node[0],new_node[0]) )
            fprintf(stderr,"Could not find adjacent child 1.\n");
         if (!find_adjacent_child(this_tri->adjacent[0]->adjacent[0],new_tri[1],new_node[0],this_tri->node[1]) )
            fprintf(stderr,"Could not find adjacent child 2.\n");
      }

      /* do the same for parent's side 1 */
      if (this_tri->adjacent[1] && has_adjacent_been_split[1]) {
         if (!find_adjacent_child(this_tri->adjacent[1]->adjacent[0],new_tri[1],this_tri->node[1],new_node[1]) )
            fprintf(stderr,"Could not find adjacent child 3.\n");
         if (!find_adjacent_child(this_tri->adjacent[1]->adjacent[0],new_tri[2],new_node[1],this_tri->node[2]) )
            fprintf(stderr,"Could not find adjacent child 4.\n");
      }

      /* and for parent's side 2 */
      if (this_tri->adjacent[2] && has_adjacent_been_split[2]) {
         if (!find_adjacent_child(this_tri->adjacent[2]->adjacent[0],new_tri[2],this_tri->node[2],new_node[2]) )
            fprintf(stderr,"Could not find adjacent child 5.\n");
         if (!find_adjacent_child(this_tri->adjacent[2]->adjacent[0],new_tri[0],new_node[2],this_tri->node[0]) )
            fprintf(stderr,"Could not find adjacent child 6.\n");
      }

      /* new tri's midpoints will not be set this recursion level, set to NULL */
      for (i=0; i<4; i++) for (j=0; j<3; j++) new_tri[i]->midpoint[j] = NULL;

      /* calculate triangle's area, if below either threshhold, set it up so it will not split */
      for (i=0; i<4; i++) new_tri[i]->splittable = TRUE;
      if (use_thresh || use_dist) {
         for (i=0; i<4; i++) {
            temp_area = find_area(new_tri[i]);
            if (use_thresh) {
               if (temp_area < area_thresh) {
                  /* this is a flag to the splitter, do not split further */
                  new_tri[i]->splittable = FALSE;
               }
            }
            if (use_dist) {
               // if (sqrt(temp_area)/length(from(viewp,new_tri[i]->node[0]->loc)) < distance_thresh) {
               if (sqrt(temp_area)/find_tri_dist(new_tri[i],viewp) < distance_thresh) {
                  /* this is a flag to the splitter, do not split further */
                  new_tri[i]->splittable = FALSE;
               }
            }
         }
      }

      /* add the 4 new tris to the new list */
      new_tri[3]->next_tri = new_tri_head;
      new_tri[2]->next_tri = new_tri[3];
      new_tri[1]->next_tri = new_tri[2];
      new_tri[0]->next_tri = new_tri[1];
      new_tri_head = new_tri[0];

      /* Important: set this now-split parent triangle's first adjacent pointer to the
       * first of the 4 new triangles created, this information will be useful later */
      this_tri->adjacent[0] = new_tri[0];

      /* we are now done with this parent tri, we may reference it later, though */
      this_tri = this_tri->next_tri;
   }

   // finally, call the sphericalizing routine
   if (force_sphere) make_sphere(new_tri_head);

   return(new_tri_head);
}


/*
 * split_tri_5
 *
 * Takes a linked list of triangles and splits any number
 * of triangles into 4 new triangles.
 * It only adds one new node for each edge, regardless of how many triangles
 * share that edge.
 * All new nodes are placed either at the midpoint of their edges, or 
 * at a location corresponding to the spline surface through the close
 * nodes.
 * Only after the new topology is created are the nodes moved from their
 * initial positions.
 */
tri_pointer split_tri_5 (int depth, tri_pointer tri_head) {

   int h,i,j,k;
   int adj_side = -1;
   int num_new_tri;
   int tri_cnt = 0;
   int node_cnt = 0;
   int has_adjacent_been_split[3];
   double temp_area;
   tri_pointer this_tri;
   tri_pointer new_tri_head = NULL;
   tri_pointer new_tri[4];
   node_ptr this_node;
   node_ptr new_node[3];

   fprintf(stderr,"Method 3, depth = %d\n",depth);
   j = -1; k = -1;

   // count the nodes (we need this to set the index)
   node_cnt = count_nodes();

   // first, compute normals for all triangles (use best method)
   if (depth == 0) (void) compute_normals_2 (tri_head,3);

   // for each triangle in the old list
   this_tri = tri_head;
   while (this_tri) {

      //fprintf(stderr,"Checking tri...%d\n",this_tri->index);

      // if the tri is too small, do not split it, just copy it to the new list
      // new logic: if the tri is flagged to not split, split it minimally to
      //    account for neighboring splittable triangles!
      if (!this_tri->splittable) {

         //fprintf(stderr,"Unsplittable tri %d\n",this_tri->index); fflush(stderr);

         // first, see how many neighbors will be split
         num_new_tri = 1;
         for (i=0; i<3; i++) new_node[i] = NULL;
         for (i=0; i<3; i++) {

            // check to see if a midpoint split already exists, use it if it does
            if (this_tri->midpoint[i]) {

               // that edge has already been split, use the point
               has_adjacent_been_split[i] = TRUE;
               new_node[i] = this_tri->midpoint[i];
               num_new_tri++;

               // if (this_tri->adjacent[i]->splittable) {
                  // fprintf(stderr,"     corresponding to splittable tri %d\n",this_tri->adjacent[i]->index); fflush(stderr);
               // } else {
                  // fprintf(stderr,"     corresponding to unsplittable tri %d\n",this_tri->adjacent[i]->index); fflush(stderr);
               // }

            } else {

               // if it doesn't, see if it eventually will
               has_adjacent_been_split[i] = FALSE;
               j = mod(i+1,3);
               k = mod(j+1,3);

               // fprintf(stderr,"   creating midpoint for side %d\n",i); fflush(stderr);
               if (!this_tri->adjacent[i]) {

                  // there is no adjacent triangle, don't split this edge

               } else if (this_tri->adjacent[i]->splittable) {

                  // fprintf(stderr,"     corresponding to splittable tri %d\n",this_tri->adjacent[i]->index); fflush(stderr);

                  // there is an adjacent triangle that can be split
                  if (use_spline)
                     new_node[i] = create_midpoint_spline(this_tri->node[i],this_tri->node[j],&node_cnt);
                  else
                     new_node[i] = create_midpoint(this_tri->node[i],this_tri->node[j]);

                  // find the index of the shared side from the adjacent tri's entry
                  for (j=0; j<3; j++) {
                     if (this_tri->node[i] == this_tri->adjacent[i]->node[j]) {
                        adj_side = mod(j+2,3);
                        break;
                     }
                  }
                  j = mod(adj_side+1,3);
                  this_tri->adjacent[i]->midpoint[adj_side] = new_node[i];
                  // fprintf(stderr,"      and told adjacent tri %d that its side %d has a midpoint already\n",i,adj_side);
                  num_new_tri++;

               } else {

                  // fprintf(stderr,"     corresponding to unsplittable tri %d\n",this_tri->adjacent[i]->index); fflush(stderr);

                  // there is an adjacent triangle, but it's flagged to not split,
                  //    so don't split this side at all

               }
            }
         }

         // now, based on how many tris we have, create them
         // fprintf(stderr,"  splitting into %d tris\n",num_new_tri); fflush(stderr);

         // create the N new triangles
         for (i=0; i<num_new_tri; i++) {
            new_tri[i] = alloc_new_tri();
            new_tri[i]->splittable = FALSE;
         }

         if (num_new_tri == 4) {

            // treat it just like a regular split element, but flag all 4
            //    children as non-splittable

            // no change to old method
            // by default, set all external adjacent pointers to NULL, fix later
            new_tri[0]->node[0] = this_tri->node[0];
            new_tri[0]->node[1] = new_node[0];
            new_tri[0]->node[2] = new_node[2];
            new_tri[0]->adjacent[0] = NULL;
            new_tri[0]->adjacent[1] = new_tri[3];
            new_tri[0]->adjacent[2] = NULL;
            new_tri[0]->index = tri_cnt++;
            new_tri[1]->node[0] = new_node[0];
            new_tri[1]->node[1] = this_tri->node[1];
            new_tri[1]->node[2] = new_node[1];
            new_tri[1]->adjacent[0] = NULL;
            new_tri[1]->adjacent[1] = NULL;
            new_tri[1]->adjacent[2] = new_tri[3];
            new_tri[1]->index = tri_cnt++;
            new_tri[2]->node[0] = new_node[2];
            new_tri[2]->node[1] = new_node[1];
            new_tri[2]->node[2] = this_tri->node[2];
            new_tri[2]->adjacent[0] = new_tri[3];
            new_tri[2]->adjacent[1] = NULL;
            new_tri[2]->adjacent[2] = NULL;
            new_tri[2]->index = tri_cnt++;
            new_tri[3]->node[0] = new_node[1];
            new_tri[3]->node[1] = new_node[2];
            new_tri[3]->node[2] = new_node[0];
            new_tri[3]->adjacent[0] = new_tri[2];
            new_tri[3]->adjacent[1] = new_tri[0];
            new_tri[3]->adjacent[2] = new_tri[1];
            new_tri[3]->index = tri_cnt++;

            // If parent had no adjacent tri, children on that edge have none, either.
            // BUT, if the parent DID have an adjacent tri, children will also.
            // So, if the parent had an adjacent tri, and it has been split, find the
            // children's adjacent tris among the parent's adjacent's children
            if (this_tri->adjacent[0] && has_adjacent_been_split[0]) {

               // then use the child tri pointed to by this_tri->adjacent[0]->adjacent[0]
               // to begin searching for the tris adjacent to the 2 new ones on this side

               if (!find_adjacent_child(this_tri->adjacent[0]->adjacent[0],
                     new_tri[0],this_tri->node[0],new_node[0]) )
                  fprintf(stderr,"Could not find adjacent child 1.\n");
               if (!find_adjacent_child(this_tri->adjacent[0]->adjacent[0],
                     new_tri[1],new_node[0],this_tri->node[1]) )
                  fprintf(stderr,"Could not find adjacent child 2.\n");
            }

            // do the same for parent's side 1
            if (this_tri->adjacent[1] && has_adjacent_been_split[1]) {
               if (!find_adjacent_child(this_tri->adjacent[1]->adjacent[0],
                     new_tri[1],this_tri->node[1],new_node[1]) )
                  fprintf(stderr,"Could not find adjacent child 3.\n");
               if (!find_adjacent_child(this_tri->adjacent[1]->adjacent[0],
                     new_tri[2],new_node[1],this_tri->node[2]) )
                  fprintf(stderr,"Could not find adjacent child 4.\n");
            }

            // and for parent's side 2
            if (this_tri->adjacent[2] && has_adjacent_been_split[2]) {
               if (!find_adjacent_child(this_tri->adjacent[2]->adjacent[0],
                     new_tri[2],this_tri->node[2],new_node[2]) )
                  fprintf(stderr,"Could not find adjacent child 5.\n");
               if (!find_adjacent_child(this_tri->adjacent[2]->adjacent[0],
                     new_tri[0],new_node[2],this_tri->node[0]) )
                  fprintf(stderr,"Could not find adjacent child 6.\n");
            }

         // reset node parameters here, too; new nodes first
#ifdef CONN
         add_conn_tri (new_node[0], new_tri[0], 1);
         add_conn_tri (new_node[0], new_tri[1], 0);
         add_conn_tri (new_node[0], new_tri[3], 2);
         //nc = new_node[0]->num_conn;
         //new_node[0]->conn_tri[nc] = new_tri[0];
         //new_node[0]->conn_tri_node[nc++] = 1;
         //new_node[0]->conn_tri[nc] = new_tri[1];
         //new_node[0]->conn_tri_node[nc++] = 0;
         //new_node[0]->conn_tri[nc] = new_tri[3];
         //new_node[0]->conn_tri_node[nc++] = 2;
         //new_node[0]->num_conn = nc;

         add_conn_tri (new_node[1], new_tri[1], 2);
         add_conn_tri (new_node[1], new_tri[2], 1);
         add_conn_tri (new_node[1], new_tri[3], 0);
         //nc = new_node[1]->num_conn;
         //new_node[1]->conn_tri[nc] = new_tri[1];
         //new_node[1]->conn_tri_node[nc++] = 2;
         //new_node[1]->conn_tri[nc] = new_tri[2];
         //new_node[1]->conn_tri_node[nc++] = 1;
         //new_node[1]->conn_tri[nc] = new_tri[3];
         //new_node[1]->conn_tri_node[nc++] = 0;
         //new_node[1]->num_conn = nc;

         add_conn_tri (new_node[2], new_tri[2], 0);
         add_conn_tri (new_node[2], new_tri[0], 2);
         add_conn_tri (new_node[2], new_tri[3], 1);
         //nc = new_node[2]->num_conn;
         //new_node[2]->conn_tri[nc] = new_tri[2];
         //new_node[2]->conn_tri_node[nc++] = 0;
         //new_node[2]->conn_tri[nc] = new_tri[0];
         //new_node[2]->conn_tri_node[nc++] = 2;
         //new_node[2]->conn_tri[nc] = new_tri[3];
         //new_node[2]->conn_tri_node[nc++] = 1;
         //new_node[2]->num_conn = nc;

         // and the old nodes; they need to have their pointer to the old parent tri
         //   removed and replaced with a pointer to the new child tri
         for (i=0; i<3; i++) {
            this_node = this_tri->node[i];
            for (j=0; j<this_node->num_conn; j++) {
               if (this_node->conn_tri[j] == this_tri) {
                  this_node->conn_tri[j] = new_tri[i];
                  this_node->conn_tri_node[j] = i;
               }
            }
         }
#endif

            // new tri's midpts will not be set this recursion level, set to NULL
            for (i=0; i<4; i++)
               for (j=0; j<3; j++)
                  new_tri[i]->midpoint[j] = NULL;

            // Important: set this now-split parent triangle's first adjacent
            //   pointer to the first of the 4 new triangles created, this
            //   information will be useful later
            this_tri->adjacent[0] = new_tri[0];

         } else if (num_new_tri == 3) {

            // first, which edge doesn't have a new_node?
            for (h=0; h<3; h++) if (!new_node[h]) {
               i = (h+1)%3;
               j = (h+2)%3;
               k = h;
               break;
            }

            // geometry of the three triangles
            new_tri[0]->node[0] = this_tri->node[i];
            new_tri[0]->node[1] = new_node[i];
            new_tri[0]->node[2] = this_tri->node[k];
            new_tri[0]->adjacent[0] = NULL;
            new_tri[0]->adjacent[1] = new_tri[2];
            new_tri[0]->adjacent[2] = NULL;
            new_tri[0]->index = tri_cnt++;
            new_tri[1]->node[0] = new_node[i];
            new_tri[1]->node[1] = this_tri->node[j];
            new_tri[1]->node[2] = new_node[j];
            new_tri[1]->adjacent[0] = NULL;
            new_tri[1]->adjacent[1] = NULL;
            new_tri[1]->adjacent[2] = new_tri[2];
            new_tri[1]->index = tri_cnt++;
            new_tri[2]->node[0] = new_node[i];
            new_tri[2]->node[1] = new_node[j];
            new_tri[2]->node[2] = this_tri->node[k];
            new_tri[2]->adjacent[0] = new_tri[1];
            new_tri[2]->adjacent[1] = NULL;
            new_tri[2]->adjacent[2] = new_tri[0];
            new_tri[2]->index = tri_cnt++;

            // now, set those adjacent pointers
            if (this_tri->adjacent[i] && has_adjacent_been_split[i]) {
               if (!find_adjacent_child(this_tri->adjacent[i]->adjacent[0],
                     new_tri[0],this_tri->node[i],new_node[i]) )
                  fprintf(stderr,"Could not find adjacent child 1.\n");
               if (!find_adjacent_child(this_tri->adjacent[i]->adjacent[0],
                     new_tri[1],new_node[i],this_tri->node[j]) )
                  fprintf(stderr,"Could not find adjacent child 2.\n");
            }
            if (this_tri->adjacent[j] && has_adjacent_been_split[j]) {
               if (!find_adjacent_child(this_tri->adjacent[j]->adjacent[0],
                     new_tri[1],this_tri->node[j],new_node[j]) )
                  fprintf(stderr,"Could not find adjacent child 3.\n");
               if (!find_adjacent_child(this_tri->adjacent[j]->adjacent[0],
                     new_tri[2],new_node[j],this_tri->node[k]) )
                  fprintf(stderr,"Could not find adjacent child 4.\n");
            }
            if (this_tri->adjacent[k] && has_adjacent_been_split[k]) {
               fprintf(stderr,"Side k should not have been split!\n");
            }

#ifdef CONN
            // reset node parameters here, too;
            // first, the new nodes
            add_conn_tri (new_node[i], new_tri[0], 1);
            add_conn_tri (new_node[i], new_tri[1], 0);
            add_conn_tri (new_node[i], new_tri[2], 0);
            //nc = new_node[i]->num_conn;
            //new_node[i]->conn_tri[nc] = new_tri[0];
            //new_node[i]->conn_tri_node[nc++] = 1;
            //new_node[i]->conn_tri[nc] = new_tri[1];
            //new_node[i]->conn_tri_node[nc++] = 0;
            //new_node[i]->conn_tri[nc] = new_tri[2];
            //new_node[i]->conn_tri_node[nc++] = 0;
            //new_node[i]->num_conn = nc;

            add_conn_tri (new_node[j], new_tri[1], 2);
            add_conn_tri (new_node[j], new_tri[2], 1);
            //nc = new_node[j]->num_conn;
            //new_node[j]->conn_tri[nc] = new_tri[1];
            //new_node[j]->conn_tri_node[nc++] = 2;
            //new_node[j]->conn_tri[nc] = new_tri[2];
            //new_node[j]->conn_tri_node[nc++] = 1;
            //new_node[j]->num_conn = nc;

            // then the node across from the new node (add 2 tris, remove old one)
            this_node = this_tri->node[k];
            add_conn_tri (this_node, new_tri[0], 2);
            //nc = this_node->num_conn;
            //this_node->conn_tri[nc] = new_tri[0];
            //this_node->conn_tri_node[nc++] = 2;
            //this_node->num_conn = nc;
            for (h=0; h<this_node->num_conn; h++) {
               if (this_node->conn_tri[h] == this_tri) {
                  this_node->conn_tri[h] = new_tri[2];
                  this_node->conn_tri_node[h] = 2;
               }
            }

            // the last 2 nodes need to have their pointer to the old parent tri
            //   removed and replaced with a pointer to the new child tri
            this_node = this_tri->node[i];
            for (h=0; h<this_node->num_conn; h++) {
               if (this_node->conn_tri[h] == this_tri) {
                  this_node->conn_tri[h] = new_tri[0];
                  this_node->conn_tri_node[h] = 0;
               }
            }
            this_node = this_tri->node[j];
            for (h=0; h<this_node->num_conn; h++) {
               if (this_node->conn_tri[h] == this_tri) {
                  this_node->conn_tri[h] = new_tri[1];
                  this_node->conn_tri_node[h] = 1;
               }
            }
#endif

            // set that other stuff
            for (i=0; i<3; i++)
               for (j=0; j<3; j++)
                  new_tri[i]->midpoint[j] = NULL;

            this_tri->adjacent[0] = new_tri[0];

         } else if (num_new_tri == 2) {

            // first, which edge has the new_node?
            for (h=0; h<3; h++) if (new_node[h]) {
               i = h;
               j = (h+1)%3;
               k = (h+2)%3;
               break;
            }

            // geometry of the two triangles
            new_tri[0]->node[0] = this_tri->node[i];
            new_tri[0]->node[1] = new_node[i];
            new_tri[0]->node[2] = this_tri->node[k];
            new_tri[0]->adjacent[0] = NULL;
            new_tri[0]->adjacent[1] = new_tri[1];
            new_tri[0]->adjacent[2] = NULL;
            new_tri[0]->index = tri_cnt++;
            new_tri[1]->node[0] = new_node[i];
            new_tri[1]->node[1] = this_tri->node[j];
            new_tri[1]->node[2] = this_tri->node[k];
            new_tri[1]->adjacent[0] = NULL;
            new_tri[1]->adjacent[1] = NULL;
            new_tri[1]->adjacent[2] = new_tri[0];
            new_tri[1]->index = tri_cnt++;

            // now, set those adjacent pointers
            if (this_tri->adjacent[i] && has_adjacent_been_split[i]) {
               if (!find_adjacent_child(this_tri->adjacent[i]->adjacent[0],
                     new_tri[0],this_tri->node[i],new_node[i]) )
                  fprintf(stderr,"Could not find adjacent child 1.\n");
               if (!find_adjacent_child(this_tri->adjacent[i]->adjacent[0],
                     new_tri[1],new_node[i],this_tri->node[j]) )
                  fprintf(stderr,"Could not find adjacent child 2.\n");
            }
            if (this_tri->adjacent[j] && has_adjacent_been_split[j]) {
               fprintf(stderr,"Side j should not have been split!\n");
            }
            if (this_tri->adjacent[k] && has_adjacent_been_split[k]) {
               fprintf(stderr,"Side k should not have been split!\n");
            }

#ifdef CONN
            // reset node parameters here, too;
            // first, the new node
            add_conn_tri (new_node[i], new_tri[0], 1);
            add_conn_tri (new_node[i], new_tri[1], 0);
            //nc = new_node[i]->num_conn;
            //new_node[i]->conn_tri[nc] = new_tri[0];
            //new_node[i]->conn_tri_node[nc++] = 1;
            //new_node[i]->conn_tri[nc] = new_tri[1];
            //new_node[i]->conn_tri_node[nc++] = 0;
            //new_node[i]->num_conn = nc;

            // then the node across from the new node (add 2 tris, remove old one)
            this_node = this_tri->node[k];
            add_conn_tri (this_node, new_tri[0], 2);
            //nc = this_node->num_conn;
            //this_node->conn_tri[nc] = new_tri[0];
            //this_node->conn_tri_node[nc++] = 2;
            //this_node->num_conn = nc;
            for (h=0; h<this_node->num_conn; h++) {
               if (this_node->conn_tri[h] == this_tri) {
                  this_node->conn_tri[h] = new_tri[1];
                  this_node->conn_tri_node[h] = 2;
               }
            }

            // the last 2 nodes need to have their pointer to the old parent tri
            //   removed and replaced with a pointer to the new child tri
            this_node = this_tri->node[i];
            for (h=0; h<this_node->num_conn; h++) {
               if (this_node->conn_tri[h] == this_tri) {
                  this_node->conn_tri[h] = new_tri[0];
                  this_node->conn_tri_node[h] = 0;
               }
            }
            this_node = this_tri->node[j];
            for (h=0; h<this_node->num_conn; h++) {
               if (this_node->conn_tri[h] == this_tri) {
                  this_node->conn_tri[h] = new_tri[1];
                  this_node->conn_tri_node[h] = 1;
               }
            }
#endif

            // set that other stuff
            for (i=0; i<2; i++)
               for (j=0; j<3; j++)
                  new_tri[i]->midpoint[j] = NULL;

            this_tri->adjacent[0] = new_tri[0];

         } else {	// if num_new_tri==1

            // then just pass this tri on to the next stage
            // new_tri[0] = alloc_new_tri();
            // new_tri[0]->splittable = FALSE;
            new_tri[0]->node[0] = this_tri->node[0];
            new_tri[0]->node[1] = this_tri->node[1];
            new_tri[0]->node[2] = this_tri->node[2];
            new_tri[0]->adjacent[0] = NULL;
            new_tri[0]->adjacent[1] = NULL;
            new_tri[0]->adjacent[2] = NULL;
            new_tri[0]->index = tri_cnt++;
            // new_tri[0]->next_tri = new_tri_head;
            // new_tri_head = new_tri[0];
            // this_tri = this_tri->next_tri;
            // continue;

            //  shouldn't we set adjacent pointers?
            new_tri[0]->adjacent[0] = this_tri->adjacent[0];
            new_tri[0]->adjacent[1] = this_tri->adjacent[1];
            new_tri[0]->adjacent[2] = this_tri->adjacent[2];
            this_tri->adjacent[0] = new_tri[0];

#ifdef CONN
            // reset node parameters here, too;
            // the old nodes need to have their pointer to the old parent tri
            //   removed and replaced with a pointer to the new child tri
            for (i=0; i<3; i++) {
               this_node = this_tri->node[i];
               for (j=0; j<this_node->num_conn; j++) {
                  if (this_node->conn_tri[j] == this_tri) {
                     this_node->conn_tri[j] = new_tri[0];
                     this_node->conn_tri_node[j] = i;
                  }
               }
            }
#endif

         }

         // add all of the tris to the list
         // why backwards? oh, well.
         for (i=num_new_tri-1; i>-1; i--) {
            new_tri[i]->next_tri = new_tri_head;
            new_tri_head = new_tri[i];
         }

         // and jump to the next parent tri
         this_tri = this_tri->next_tri;
         continue;
      }

      //--------------------------------------------------------------------------
      // if this is a normal, splittable triangle, create 4 child triangles

      //fprintf(stderr,"Splittable tri %d\n",this_tri->index); fflush(stderr);

      // choose the three edge split points
      for (i=0; i<3; i++) new_node[i] = NULL;

      for (i=0; i<3; i++) {

         //fprintf(stderr,"  when checking side %d...\n",i);

         // check to see if a midpoint split already exists, use it if it does
         if (this_tri->midpoint[i]) {

            has_adjacent_been_split[i] = TRUE;
            //fprintf(stderr,"    already have midpoint on side %d\n",i); fflush(stderr);
            new_node[i] = this_tri->midpoint[i];

         } else {

            // if it doesn't, choose one and send the information to the
            //   adjacent tri, if there is one
            has_adjacent_been_split[i] = FALSE;
            j = mod(i+1,3);

            //fprintf(stderr,"    creating midpoint for side %d\n",i); fflush(stderr);
            // use the better method for determining a new midpoint
            k = mod(j+1,3);

            if (!this_tri->adjacent[i]) {

               //fprintf(stderr,"    there is no adjacent tri %d\n",i); fflush(stderr);

               /* there is no adjacent triangle */
               if (clamp_edges) {
                  // clamp the edge, use the mathematical midpoint
                  new_node[i] = create_midpoint(this_tri->node[i],this_tri->node[j]);
               } else {
                  // let the edge get jagged, use the current tri's 3 points
                  // not for this new version---all nodes get perturbed LATER
                  //new_node[i] = create_midpoint_3(depth,this_tri->node[i],
                  //                                this_tri->node[j],
                  //                                this_tri->node[k]);
                  if (use_spline)
                     new_node[i] = create_midpoint_spline(this_tri->node[i],this_tri->node[j],&node_cnt);
                  else
                     new_node[i] = create_midpoint(this_tri->node[i],this_tri->node[j]);
               }

            } else if (this_tri->adjacent[i]->splittable) {

               //fprintf(stderr,"    corresponding to splittable tri %d\n",this_tri->adjacent[i]->index); fflush(stderr);

               // there is an adjacent triangle that can be split
               // all nodes get perturbed LATER
               if (use_spline)
                  new_node[i] = create_midpoint_spline(this_tri->node[i],this_tri->node[j],&node_cnt);
               else
                  new_node[i] = create_midpoint(this_tri->node[i],this_tri->node[j]);

               // find the index of the shared side from the adjacent tri's entry
               for (h=0; h<3; h++) {
                  if (this_tri->node[i] == this_tri->adjacent[i]->node[h]) {
                     adj_side = mod(h+2,3);
                     break;
                  }
               }
               h = mod(adj_side+1,3);
               this_tri->adjacent[i]->midpoint[adj_side] = new_node[i];
               //fprintf(stderr,"      and told adjacent tri %d that its side %d has a midpoint already\n",
               //   i,adj_side);
               //fflush(stderr);


            } else {

               //fprintf(stderr,"    corresponding to unsplittable tri %d\n",this_tri->adjacent[i]->index); fflush(stderr);

               // there is an adjacent triangle, but it's set to not split;
               // So, split that edge normally, anyways! AND tell the
               //    adjacent tri that we split

               // all nodes get perturbed LATER
               if (use_spline)
                  new_node[i] = create_midpoint_spline(this_tri->node[i],this_tri->node[j],&node_cnt);
               else
                  new_node[i] = create_midpoint(this_tri->node[i],this_tri->node[j]);

               // find the index of the shared side from the adjacent tri's entry
               for (j=0; j<3; j++) {
                  if (this_tri->node[i] == this_tri->adjacent[i]->node[j]) {
                     adj_side = mod(j+2,3);
                     break;
                  }
               }
               j = mod(adj_side+1,3);
               this_tri->adjacent[i]->midpoint[adj_side] = new_node[i];

            }

         }
      }

      // create the four new triangles
      for (i=0; i<4; i++) new_tri[i] = alloc_new_tri();

      // no change to old method
      // by default, set all external adjacent pointers to NULL, fix later
      new_tri[0]->node[0] = this_tri->node[0];
      new_tri[0]->node[1] = new_node[0];
      new_tri[0]->node[2] = new_node[2];
      new_tri[0]->adjacent[0] = NULL;
      new_tri[0]->adjacent[1] = new_tri[3];
      new_tri[0]->adjacent[2] = NULL;
      new_tri[0]->index = tri_cnt++;
      new_tri[1]->node[0] = new_node[0];
      new_tri[1]->node[1] = this_tri->node[1];
      new_tri[1]->node[2] = new_node[1];
      new_tri[1]->adjacent[0] = NULL;
      new_tri[1]->adjacent[1] = NULL;
      new_tri[1]->adjacent[2] = new_tri[3];
      new_tri[1]->index = tri_cnt++;
      new_tri[2]->node[0] = new_node[2];
      new_tri[2]->node[1] = new_node[1];
      new_tri[2]->node[2] = this_tri->node[2];
      new_tri[2]->adjacent[0] = new_tri[3];
      new_tri[2]->adjacent[1] = NULL;
      new_tri[2]->adjacent[2] = NULL;
      new_tri[2]->index = tri_cnt++;
      new_tri[3]->node[0] = new_node[1];
      new_tri[3]->node[1] = new_node[2];
      new_tri[3]->node[2] = new_node[0];
      new_tri[3]->adjacent[0] = new_tri[2];
      new_tri[3]->adjacent[1] = new_tri[0];
      new_tri[3]->adjacent[2] = new_tri[1];
      new_tri[3]->index = tri_cnt++;

      // reset node parameters here, too; new nodes first
#ifdef CONN
      add_conn_tri (new_node[0], new_tri[0], 1);
      add_conn_tri (new_node[0], new_tri[1], 0);
      add_conn_tri (new_node[0], new_tri[3], 2);
      //nc = new_node[0]->num_conn;
      //new_node[0]->conn_tri[nc] = new_tri[0];
      //new_node[0]->conn_tri_node[nc++] = 1;
      //new_node[0]->conn_tri[nc] = new_tri[1];
      //new_node[0]->conn_tri_node[nc++] = 0;
      //new_node[0]->conn_tri[nc] = new_tri[3];
      //new_node[0]->conn_tri_node[nc++] = 2;
      //new_node[0]->num_conn = nc;

      add_conn_tri (new_node[1], new_tri[1], 2);
      add_conn_tri (new_node[1], new_tri[2], 1);
      add_conn_tri (new_node[1], new_tri[3], 0);
      //nc = new_node[1]->num_conn;
      //new_node[1]->conn_tri[nc] = new_tri[1];
      //new_node[1]->conn_tri_node[nc++] = 2;
      //new_node[1]->conn_tri[nc] = new_tri[2];
      //new_node[1]->conn_tri_node[nc++] = 1;
      //new_node[1]->conn_tri[nc] = new_tri[3];
      //new_node[1]->conn_tri_node[nc++] = 0;
      //new_node[1]->num_conn = nc;

      add_conn_tri (new_node[2], new_tri[2], 0);
      add_conn_tri (new_node[2], new_tri[0], 2);
      add_conn_tri (new_node[2], new_tri[3], 1);
      //nc = new_node[2]->num_conn;
      //new_node[2]->conn_tri[nc] = new_tri[2];
      //new_node[2]->conn_tri_node[nc++] = 0;
      //new_node[2]->conn_tri[nc] = new_tri[0];
      //new_node[2]->conn_tri_node[nc++] = 2;
      //new_node[2]->conn_tri[nc] = new_tri[3];
      //new_node[2]->conn_tri_node[nc++] = 1;
      //new_node[2]->num_conn = nc;

      // and the old nodes; they need to have their pointer to the old parent tri
      //   removed and replaced with a pointer to the new child tri
      for (i=0; i<3; i++) {
         this_node = this_tri->node[i];
         for (j=0; j<this_node->num_conn; j++) {
            if (this_node->conn_tri[j] == this_tri) {
               this_node->conn_tri[j] = new_tri[i];
               this_node->conn_tri_node[j] = i;
            }
         }
      }
#endif

      //fprintf(stderr,"Now checking adjacents' adjacents\n");

      /* If parent had no adjacent tri, children on that edge have none, either.
       * BUT, if the parent DID have an adjacent tri, children will also.
       * So, if the parent had an adjacent tri, and it has been split, find the children's
       * adjacent tris among the parent's adjacent's children */
      /* but what about the actual node locations? NO, this is just setup for the next step */
      if (this_tri->adjacent[0] && has_adjacent_been_split[0]) {

         /* then use the child tri pointed to by this_tri->adjacent[0]->adjacent[0]
          * to begin searching for the triangles adjacent to the 2 new ones on this side */

         if (!find_adjacent_child(this_tri->adjacent[0]->adjacent[0],new_tri[0],this_tri->node[0],new_node[0]) )
            fprintf(stderr,"Could not find adjacent child 1.\n");
         if (!find_adjacent_child(this_tri->adjacent[0]->adjacent[0],new_tri[1],new_node[0],this_tri->node[1]) )
            fprintf(stderr,"Could not find adjacent child 2.\n");
      }

      /* do the same for parent's side 1 */
      if (this_tri->adjacent[1] && has_adjacent_been_split[1]) {
         if (!find_adjacent_child(this_tri->adjacent[1]->adjacent[0],new_tri[1],this_tri->node[1],new_node[1]) )
            fprintf(stderr,"Could not find adjacent child 3.\n");
         if (!find_adjacent_child(this_tri->adjacent[1]->adjacent[0],new_tri[2],new_node[1],this_tri->node[2]) )
            fprintf(stderr,"Could not find adjacent child 4.\n");
      }

      /* and for parent's side 2 */
      if (this_tri->adjacent[2] && has_adjacent_been_split[2]) {
         if (!find_adjacent_child(this_tri->adjacent[2]->adjacent[0],new_tri[2],this_tri->node[2],new_node[2]) )
            fprintf(stderr,"Could not find adjacent child 5.\n");
         if (!find_adjacent_child(this_tri->adjacent[2]->adjacent[0],new_tri[0],new_node[2],this_tri->node[0]) )
            fprintf(stderr,"Could not find adjacent child 6.\n");
      }

      /* new tri's midpoints will not be set this recursion level, set to NULL */
      for (i=0; i<4; i++)
         for (j=0; j<3; j++)
            new_tri[i]->midpoint[j] = NULL;

      // calculate triangle's area, if below either threshhold, set it up so
      //   it will not split
      for (i=0; i<4; i++) new_tri[i]->splittable = TRUE;
      if (use_thresh || use_dist) {
         for (i=0; i<4; i++) {
            temp_area = find_area(new_tri[i]);
            if (use_thresh) {
               if (temp_area < area_thresh) {
                  /* this is a flag to the splitter, do not split further */
                  new_tri[i]->splittable = FALSE;
               }
            }
            if (use_dist) {
               // if (sqrt(temp_area)/length(from(viewp,new_tri[i]->node[0]->loc)) < distance_thresh)
               if (sqrt(temp_area)/find_tri_dist(new_tri[i],viewp) < distance_thresh) {
                  /* this is a flag to the splitter, do not split further */
                  new_tri[i]->splittable = FALSE;
               }
            }
         }
      }

      /* add the 4 new tris to the new list */
      new_tri[3]->next_tri = new_tri_head;
      new_tri[2]->next_tri = new_tri[3];
      new_tri[1]->next_tri = new_tri[2];
      new_tri[0]->next_tri = new_tri[1];
      new_tri_head = new_tri[0];

      /* Important: set this now-split parent triangle's first adjacent pointer to the
       * first of the 4 new triangles created, this information will be useful later */
      this_tri->adjacent[0] = new_tri[0];

      /* we are now done with this parent tri, we may reference it later, though */
      this_tri = this_tri->next_tri;
   }

   // -------------- here is where we perturb the nodes -------------

   // first, compute normals for all triangles, use best normal-finding
   (void) compute_normals_2 (new_tri_head,3);

   // dump normals
   //this_tri = new_tri_head;
   //while (this_tri) {
   //   fprintf(stderr,"tri %d\n",this_tri->index);
   //   fprintf(stderr,"  nodes %d %d %d\n",this_tri->node[0]->index,this_tri->node[1]->index,this_tri->node[2]->index);
   //   for (i=0; i<3; i++) {
   //      fprintf(stderr,"  norm %d  %g %g %g\n",i,this_tri->norm[i].x,this_tri->norm[i].y,this_tri->norm[i].z);
   //   }
   //   //fprintf(stderr,"  adj %d %d %d\n",this_tri->adjacent[0]->index,this_tri->adjacent[1]->index,this_tri->adjacent[2]->index);
   //   this_tri = this_tri->next_tri;
   //}

   // then, perturb each according to the normal
   this_node = node_head;
   while (this_node) {
      //fprintf(stderr,"Perturbing node %d\n",this_node->index);
      (void) move_existing_node_5(depth, this_node);
      this_node = this_node->next_node;
   }

   if (force_sphere) {
      // finally, call the sphericalizing routine
      make_sphere(new_tri_head);
   }

   return(new_tri_head);
}


/*
 * make_sphere will move all node locations in the radial direction so that
 * the resulting shape is a sphere
 *
 * It will use the sphere_rad if that number is positive, otherwise,
 * it will get the radius from the maximum radius of the current locations,
 * based on the node-weighted center of the object
 *
 * now, it also resets the normal to be the true normal
 */
int make_sphere (tri_pointer head) {

  int i,j,cnt;
  double center[3],d[3],rad,maxrad;
  tri_pointer currt;

  if (sphere_rad <= 0.) {

    // set all node indexes to -1
    currt = head;
    while (currt) {
      for (i=0; i<3; i++) currt->node[i]->index = -1;
      currt = currt->next_tri;
    }

    // find the center of the object (each node weighs 1)
    currt = head;
    cnt = 0;
    for (i=0; i<3; i++) center[i] = 0.0;
    while (currt) {
      for (i=0; i<3; i++) if (currt->node[i]->index == -1) {
        center[0] += currt->node[i]->loc.x;
        center[1] += currt->node[i]->loc.y;
        center[2] += currt->node[i]->loc.z;
        currt->node[i]->index = 1;
        cnt++;
      }
      currt = currt->next_tri;
    }
    for (i=0; i<3; i++) center[i] /= (double)cnt;
    // fprintf(stderr,"center of %d nodes is at %g %g %g\n",cnt,center[0],center[1],center[2]);

    // set all node indexes to -1
    currt = head;
    while (currt) {
      for (i=0; i<3; i++) currt->node[i]->index = -1;
      currt = currt->next_tri;
    }

    // find the current maximum radius
    currt = head;
    maxrad = 0.0;
    while (currt) {
      for (i=0; i<3; i++) if (currt->node[i]->index == -1) {
        d[0] = currt->node[i]->loc.x - center[0];
        d[1] = currt->node[i]->loc.y - center[1];
        d[2] = currt->node[i]->loc.z - center[2];
        rad = sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);
        if (rad > maxrad) maxrad = rad;
        currt->node[i]->index = 1;
      }
      currt = currt->next_tri;
    }
    // fprintf(stderr,"   max radius is %g\n",maxrad);

  } else {

    // if a sphere radius is defined on the command-line, use it instead
    maxrad = sphere_rad;
    center[0] = 0.;
    center[1] = 0.;
    center[2] = 0.;
    cnt=0;
  }

  // set all node indexes to -1
  currt = head;
  while (currt) {
    for (i=0; i<3; i++) currt->node[i]->index = -1;
    currt = currt->next_tri;
  }

  // for each node, reset its radius only to match the max radius
  currt = head;
  while (currt) {
    for (i=0; i<3; i++) {
      if (currt->node[i]->index == -1) {
        d[0] = currt->node[i]->loc.x - center[0];
        d[1] = currt->node[i]->loc.y - center[1];
        d[2] = currt->node[i]->loc.z - center[2];
        rad = sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);
        for (j=0; j<3; j++) d[j] *= maxrad/rad;
        currt->node[i]->loc.x = center[0] + d[0];
        currt->node[i]->loc.y = center[1] + d[1];
        currt->node[i]->loc.z = center[2] + d[2];
        currt->node[i]->index = 1;
      }
      // also reset the normal, if there is one
      if (currt->norm[i]) {
        currt->norm[i]->norm = norm(currt->node[i]->loc);
      }
    }
    currt = currt->next_tri;
  }

  // return the node count
  return(cnt);
}


/*
 * Will search through starting and the next 3 triangles looking for
 * an adjacent triangle to the base triangle, sharing two common nodes
 */
int find_adjacent_child(tri_pointer starting, tri_pointer base, node_ptr node1, node_ptr node2){

   int i, j;
   int base_side = -1;
   int adj_side = -1;
   int next_node;
   int found_match = FALSE;
   int match_tri = -1;
   tri_pointer matched_tri = NULL;

   //fprintf(stderr,"starting %ld, node2 %ld\n",starting,node2); fflush(stderr);
   if (!starting) return (FALSE);

   // identify which side of the base triangle needs matching
   for (j=0; j<3; j++) {
      if (node1 == base->node[j]) {
         base_side = j;
         break;
      }
   }

   // search only the starting tri and the next 2 (the one after that is the central tri)
   for (i=0; i<3; i++) if (!found_match) {
      for (j=0; j<3; j++) {
         if (node2 == starting->node[j]) {
            next_node = mod(j+1,3);
            if (node1 == starting->node[next_node]) {
               adj_side = j;
               found_match = TRUE;
               match_tri = i;
               matched_tri = starting;
               break;
               break;
            }
         }
      }
      // fprintf(stderr,"  checking...%d, adj side=%d\n",i,adj_side); fflush(stderr);
      // that child triangle didn't match, try the next
      starting = starting->next_tri;
   }

   // if a match is found, share adjacent triangle information
   if (found_match) {

      if (FALSE) {
      fprintf(stderr,"   found match, side %d of base and side %d of tri %d\n",base_side,adj_side,match_tri);
      fprintf(stderr,"      bewteen %g %g %g and %g %g %g\n",
                     node1->loc.x,node1->loc.y,node1->loc.z,
                     node2->loc.x,node2->loc.y,node2->loc.z);
      fprintf(stderr,"      base tri is %g %g to %g %g to %g %g\n",
                     base->node[0]->loc.x,base->node[0]->loc.y,
                     base->node[1]->loc.x,base->node[1]->loc.y,
                     base->node[2]->loc.x,base->node[2]->loc.y);
      fprintf(stderr,"      matched tri is %g %g to %g %g to %g %g\n",
                     matched_tri->node[0]->loc.x,matched_tri->node[0]->loc.y,
                     matched_tri->node[1]->loc.x,matched_tri->node[1]->loc.y,
                     matched_tri->node[2]->loc.x,matched_tri->node[2]->loc.y);
      }

      base->adjacent[base_side] = matched_tri;
      matched_tri->adjacent[adj_side] = base;
   }

   return found_match;
}


/*
 * Create a node identically between the two nodes
 * Use only the two endpoints which the new node will be between
 */
node_ptr create_midpoint(node_ptr node1, node_ptr node2) {

   node_ptr new_node = (NODE *)malloc(sizeof(NODE));

   /* new node is the midpoint of the two existing nodes */
   new_node->loc.x = (node1->loc.x + node2->loc.x)/2;
   new_node->loc.y = (node1->loc.y + node2->loc.y)/2;
   new_node->loc.z = (node1->loc.z + node2->loc.z)/2;

   /* grab 3 random numbers, the same as create_midpoint_2 and _3;
    * doing this allows the same random seed to produce similarly-
    * shaped surfaces */
   (void) rand();
   (void) rand();
   (void) rand();

#ifdef CONN
   new_node->num_conn = 0;
   new_node->max_conn = 0;
#endif

   // add it to the list!
   new_node->index = 0;		// what do we do here?
   new_node->next_node = node_head;
   node_head = new_node;

   return new_node;
}


/*
 * Create a node between the two nodes
 * Use only the two endpoints between which the new node will be
 * between
 */
node_ptr create_midpoint_2(int depth, node_ptr node1, node_ptr node2) {

   double dx,dy,dz,dr,length;
   node_ptr new_node = (NODE *)malloc(sizeof(NODE));

   // length of edge node1 - node2
   length = sqrt(pow(node1->loc.x-node2->loc.x,2)+pow(node1->loc.y-node2->loc.y,2)+pow(node1->loc.z-node2->loc.z,2));
   dr = length*base_shake*pow(base_exponent,depth);
   /* need to add some function of depth here, but fix my recursion problem first */
   dx = ((1.0+rand())/RAND_MAX-0.5)*dr;
   dy = ((1.0+rand())/RAND_MAX-0.5)*dr;
   dz = ((1.0+rand())/RAND_MAX-0.5)*dr;
   new_node->loc.x = dx + (node1->loc.x + node2->loc.x)/2;
   new_node->loc.y = dy + (node1->loc.y + node2->loc.y)/2;
   new_node->loc.z = dz + (node1->loc.z + node2->loc.z)/2;

#ifdef CONN
   new_node->num_conn = 0;
   new_node->max_conn = 0;
#endif

   // add it to the list!
   new_node->index = 0;		// what do we do here?
   new_node->next_node = node_head;
   node_head = new_node;

   return new_node;
}


/*
 * Create a node between the two nodes
 * Use the three corners of the working triangle, the new point
 * will be between nodes 1 and 2, with node 3 providing information
 * about the normal direction.
 */
node_ptr create_midpoint_3(int depth, node_ptr node1, node_ptr node2, node_ptr node3) {

   double d1,d2,d3,base_l;
   node_ptr new_node = (NODE *)malloc(sizeof(NODE));
   VEC r1, r2, r3;

   /* basis vector along side to be split */
   r1 = from(node1->loc,node2->loc);

   /* basis vector from new midpoint to 3rd node in working triangle */
   new_node->loc = midpt(node1->loc,node2->loc);
   r2 = from(new_node->loc,node3->loc);

   /* base length is the average length of the two planar vectors (r1,r2) */
   // base_l = (length(r1)+length(r2))/2.0;
   // now it's the geometric mean
   base_l = sqrt(length(r1)*length(r2));

   /* basis vector normal to the plane of the working triangle */
   r3 = norm(cross(r1,r2));

   d1 = base_shake*pow(base_exponent,depth)*((1.0+rand())/RAND_MAX-0.5);
   d2 = base_shake*pow(base_exponent,depth)*((1.0+rand())/RAND_MAX-0.5);
   /* make this code use normal_bias in the same way as create_midpoint_4 */
   d3 = base_l*normal_shake*pow(normal_exponent,depth)*((1.0+rand())/RAND_MAX-0.5+normal_bias);
   /* d3 = base_l*(normal_bias + normal_shake*pow(normal_exponent,depth)*((1.0+rand())/RAND_MAX-0.5)); */

   /* perturb the initial node placement in 3 basis directions */
   new_node->loc.x += d1*r1.x + d2*r2.x + d3*r3.x;
   new_node->loc.y += d1*r1.y + d2*r2.y + d3*r3.y;
   new_node->loc.z += d1*r1.z + d2*r2.z + d3*r3.z;

#ifdef CONN
   new_node->num_conn = 0;
   new_node->max_conn = 0;
#endif

   // add it to the list!
   new_node->index = 0;		// what do we do here?
   new_node->next_node = node_head;
   node_head = new_node;

   return new_node;
}


/*
 * Create a node between the two nodes
 * Use the three corners of the working triangle and the far corner
 * of the adjacent triangle, the new point will be between nodes 1
 * and 2, with node 3 and the adjacent tri's far node providing
 * information about the normal direction.
 *
 * note that basis vectors are NOT unit-length
 */
node_ptr create_midpoint_4(int depth, node_ptr node1, node_ptr node2, node_ptr node3, node_ptr node4) {

   double d1,d2,d3,base_l;
   node_ptr new_node = (NODE *)malloc(sizeof(NODE));
   VEC r1, r2, r3;//, base_pert;

   /* basis vector along side to be split */
   r1 = from(node1->loc,node2->loc);

   /* basis vector between 3rd node in working triangle and adj. tri's far corner */
   r2 = from(node4->loc,node3->loc);

   /* base length is the average length of the two planar vectors (r1,r2) */
   // base_l = (length(r1)+0.8*length(r2))/2.0;
   // now use geometric mean
   base_l = sqrt(length(r1)*0.64*length(r2));

   /* basis vector roughly normal to surface */
   r3 = norm(cross(r1,r2));

   d1 = base_shake*pow(base_exponent,depth)*((1.0+rand())/RAND_MAX-0.5);
   d2 = 0.8*base_shake*pow(base_exponent,depth)*((1.0+rand())/RAND_MAX-0.5);
   d3 = base_l*normal_shake*pow(normal_exponent,depth)*((1.0+rand())/RAND_MAX-0.5+normal_bias);

   /* Use this next format for d3 if you want super lumpy shapes */
   /* d3 = base_l*(normal_bias + normal_shake*pow(normal_exponent,depth)*((1.0+rand())/RAND_MAX-0.5)); */

   /* base location of new midpoint is average of 4 surrounding nodes
    * weighting the two closest nodes double ***NO*** */
   /* base_pert = from(midpt(node3->loc,node4->loc),midpt(node1->loc,node2->loc)), */
   new_node->loc = midpt(node1->loc,node2->loc);

   /* perturb the initial node placement in 3 basis directions */
   new_node->loc.x += d1*r1.x + d2*r2.x + d3*r3.x;
   new_node->loc.y += d1*r1.y + d2*r2.y + d3*r3.y;
   new_node->loc.z += d1*r1.z + d2*r2.z + d3*r3.z;

#ifdef CONN
   new_node->num_conn = 0;
   new_node->max_conn = 0;
#endif

   // add it to the list!
   new_node->index = 0;		// what do we do here?
   new_node->next_node = node_head;
   node_head = new_node;

   return new_node;
}


/*
 * Create a node between the two nodes - Take 5!
 *
 * Use the three corners of the working triangle and the far corner
 * of the adjacent triangle, the new point will be between nodes 1
 * and 2, with node 3 and the adjacent tri's far node providing
 * information about the normal direction.
 */
node_ptr create_midpoint_5( int depth, node_ptr node1, node_ptr node2,
                            node_ptr node3, node_ptr node4) {

   VEC r1,r2,r3;
   node_ptr new_node = (NODE *)malloc(sizeof(NODE));

   // basis vector along side to be split
   r1 = from(node1->loc,node2->loc);

   // basis vector between 3rd node in working triangle and adj. tri's far corner
   r2 = from(node4->loc,node3->loc);

   // basis vector roughly normal to surface
   r3 = norm(cross(r1,r2));

   // initial location of node
   new_node->loc = midpt(node1->loc,node2->loc);

   // perturb according to normal and depth, only
   (void) perturb_node_5 (&(new_node->loc), depth, r3);

#ifdef CONN
   new_node->num_conn = 0;
   new_node->max_conn = 0;
#endif

   // add it to the list!
   new_node->index = 0;		// what do we do here?
   new_node->next_node = node_head;
   node_head = new_node;

   return new_node;
}


/*
 * Move an existing node
 */
void move_existing_node_5 (int depth, node_ptr this) {

   int really_perturb = TRUE;
   int i,corner;
   //VEC normal;
   double len,temp[2];
   tri_pointer test_tri;

   // is this node on a clamped edge?
   if (clamp_edges) {
      // check all connected triangles to see if we're on a border
      //fprintf(stderr,"testing %d conn tris\n",this->num_conn);
      for (i=0; i<this->num_conn; i++) {
         test_tri = this->conn_tri[i];
         // check the two edges shared by this node
         corner = this->conn_tri_node[i];
         //fprintf(stderr,"  testing sides %d %d\n",corner,(corner+2)%3);
         if (test_tri->adjacent[corner] == NULL) really_perturb = FALSE;
         corner = (corner+2)%3;
         if (test_tri->adjacent[corner] == NULL) really_perturb = FALSE;
      }
   }

   if (really_perturb) {

      // find the normal vector
      test_tri = this->conn_tri[0];
      corner = this->conn_tri_node[0];

      // scale it by some function of its area?

      // perturb according to normal and depth, only
      // DANGER - what if a norm doesn't exist?
      (void) perturb_node_5 (&(this->loc), depth, test_tri->norm[corner]->norm);

   } else {

      // call the same number of random numbers!
      // this doesn't seem to work.
      len = 1.;
      while (len > 0.25) {
         temp[0] = rand()/(RAND_MAX+1.0) - 0.5;
         temp[1] = rand()/(RAND_MAX+1.0) - 0.5;
         len = temp[0]*temp[0]+temp[1]*temp[1];
      }
      len = rand()/(RAND_MAX+1.0);

   }

   return;
}


/*
 * Using simply a scaled normal vector (r3) and a depth, perturb a location;
 * Supports Gaussian random numbers (future)
 *
 * Meant to be called by create_midpoint_5
 */
void perturb_node_5 (VEC *loc, int depth, VEC r3) {

   double l01,len,temp[2],d1,d2,d3;
   double t1,t2;
   VEC r1,r2;

   // from the normal vector, create two more basis vectors
   l01 = sqrt(r3.x*r3.x+r3.y*r3.y+r3.z*r3.z);
   // first guess
   r1.x = 0.0;
   r1.y = 1.0;
   r1.z = 0.0;
   // find the component of this that is || to r3
   // make sure r1 and r3 aren't even close to parallel
   len = r1.x*r3.x+r1.y*r3.y+r1.z*r3.z;
   if (l01-fabs(len) < 0.01) {
      r1.x = 1.0;
      r1.y = 0.0;
      r1.z = 0.0;
      len = r1.x*r3.x+r1.y*r3.y+r1.z*r3.z;
   }
   // and subtract it off
   r1.x -= len*r3.x;
   r1.y -= len*r3.y;
   r1.z -= len*r3.z;
   // finally, normalize it
   len = sqrt(r1.x*r1.x+r1.y*r1.y+r1.z*r1.z) / l01;
   r1.x /= len;
   r1.y /= len;
   r1.z /= len;
   // and determine the third from a cross product
   r2.x = r3.y*r1.z - r1.y*r3.z;
   r2.y = r3.z*r1.x - r1.z*r3.x;
   r2.z = r3.x*r1.y - r1.x*r3.y;
   len = sqrt(r2.x*r2.x+r2.y*r2.y+r2.z*r2.z) / l01;
   r2.x /= len;
   r2.y /= len;
   r2.z /= len;
   // now, we have a basis of vectors, each of length l01!

   //fprintf(stderr,"length of r1 is %g\n",sqrt(r1.x*r1.x+r1.y*r1.y+r1.z*r1.z));
   //fprintf(stderr,"length of r2 is %g\n",sqrt(r2.x*r2.x+r2.y*r2.y+r2.z*r2.z));
   //fprintf(stderr,"length of r3 is %g\n",sqrt(r3.x*r3.x+r3.y*r3.y+r3.z*r3.z));

   if (use_gaussian_random) {

      t1 = (double)(rand())/(double)(RAND_MAX);
      t2 = (double)(rand())/(double)(RAND_MAX);
      d1 = (double)(sqrt(-2.*log(t1))*cos(2.*M_PI*t2));
      d2 = (double)(sqrt(-2.*log(t1))*sin(2.*M_PI*t2));
      // after this, both d1 and d2 are Gaussian random numbers with
      //   mean 0 and std deviation of 1
      d1 *= 0.5*base_shake*pow(base_exponent,depth);
      d2 *= 0.5*base_shake*pow(base_exponent,depth);

      t1 = (double)(rand())/(double)(RAND_MAX);
      t2 = (double)(rand())/(double)(RAND_MAX);
      d3 = (double)(sqrt(-2.*log(t1))*cos(2.*M_PI*t2));
      d3 = 0.5*(d3+normal_bias)*normal_shake*pow(normal_exponent,depth);

   } else {

   // determine the perturbation distances, making sure that the
   //  planar distribution is circular (_4 was rectangular)
   len = 1.;
   while (len > 0.25) {
      temp[0] = rand()/(RAND_MAX+1.0) - 0.5;
      temp[1] = rand()/(RAND_MAX+1.0) - 0.5;
      len = temp[0]*temp[0]+temp[1]*temp[1];
   }
   d1 = base_shake*pow(base_exponent,depth);
   d2 = d1*temp[1];
   d1 *= temp[0];

   d3 = normal_shake*pow(normal_exponent,depth)*((1.0+rand())/RAND_MAX-0.5+normal_bias);

   }

   // finally, perturb the node locations
   loc->x += d1*r1.x + d2*r2.x + d3*r3.x;
   loc->y += d1*r1.y + d2*r2.y + d3*r3.y;
   loc->z += d1*r1.z + d2*r2.z + d3*r3.z;

   return;
}


/*
 * Create a node between the two nodes - Take 6!
 *
 * For this one, just use the two node locations and their normals and
 * fit a cubic spline between them! This is so cool. No perturbations
 * away from the spline fit are done.
 *
 * Routine taken from vort3d/split.c:find_node_position_using_spline()
 */
node_ptr create_midpoint_spline (node_ptr n1, node_ptr n2, int *node_cnt) {

   int corner;
   double dl[3],norm1[3],norm2[3],fp[2][3][3],p1[3],p2[3],a[4];//,len;
   node_ptr new_node = (NODE *)malloc(sizeof(NODE));
   tri_pointer this_tri = NULL;

   // compute the length of the edge
   dl[0] = n2->loc.x - n1->loc.x;
   dl[1] = n2->loc.y - n1->loc.y;
   dl[2] = n2->loc.z - n1->loc.z;
   //len = sqrt(dl[0]*dl[0] + dl[1]*dl[1] + dl[2]*dl[2]);

   // find the normals

   // node 1
   for (int i=0; i<3; i++) norm1[i] = 0.0;
   for (int i=0; i<n1->num_conn; i++) {
      this_tri = n1->conn_tri[i];
      corner = n1->conn_tri_node[0];
      if (this_tri->norm[corner] == NULL) {
         // just use flat normal vector
         VEC e1 = from(this_tri->node[0]->loc,this_tri->node[1]->loc);
         VEC e2 = from(this_tri->node[0]->loc,this_tri->node[2]->loc);
         VEC tri_norm = cross(e1,e2);
         tri_norm = norm(tri_norm);
         norm1[0] += tri_norm.x;
         norm1[1] += tri_norm.y;
         norm1[2] += tri_norm.z;
      } else {
         // weigh in this normal
         norm1[0] += this_tri->norm[corner]->norm.x;
         norm1[1] += this_tri->norm[corner]->norm.y;
         norm1[2] += this_tri->norm[corner]->norm.z;
      }
   }
   (void) norm3(norm1);

   // node 2
   for (int i=0; i<3; i++) norm2[i] = 0.0;
   for (int i=0; i<n2->num_conn; i++) {
      this_tri = n2->conn_tri[i];
      corner = n2->conn_tri_node[0];
      if (this_tri->norm[corner] == NULL) {
         // just use flat normal vector
         VEC e1 = from(this_tri->node[0]->loc,this_tri->node[1]->loc);
         VEC e2 = from(this_tri->node[0]->loc,this_tri->node[2]->loc);
         VEC tri_norm = cross(e1,e2);
         tri_norm = norm(tri_norm);
         norm2[0] += tri_norm.x;
         norm2[1] += tri_norm.y;
         norm2[2] += tri_norm.z;
      } else {
         // weigh in this normal
         norm2[0] += this_tri->norm[corner]->norm.x;
         norm2[1] += this_tri->norm[corner]->norm.y;
         norm2[2] += this_tri->norm[corner]->norm.z;
      }
   }
   (void) norm3(norm2);


   // compute the tangential operator for each node (P = I - nn^T)
   for (int i=0; i<3; i++) {
      for (int j=0; j<3; j++) {
         fp[0][i][j] = 0.;
         fp[1][i][j] = 0.;
      }
      fp[0][i][i] = 1.;
      fp[1][i][i] = 1.;
      for (int j=0; j<3; j++) {
         fp[0][i][j] -= norm1[i]*norm1[j];
         fp[1][i][j] -= norm2[i]*norm2[j];
      }
   }
   // find the vector product of each of these with dl
   for (int i=0; i<3; i++) {
      p1[i] = 0.;
      p2[i] = 0.;
   }
   for (int i=0; i<3; i++) {
      for (int j=0; j<3; j++) {
         p1[i] += dl[j]*fp[0][i][j];
         p2[i] += dl[j]*fp[1][i][j];
      }
   }

   // compute a spline for each coordinate axis
   a[0] = n1->loc.x;
   a[1] = p1[0];
   a[2] = 3. * (n2->loc.x - n1->loc.x) - (p2[0] + 2. * p1[0]);
   a[3] = 2. * (n1->loc.x - n2->loc.x) + (p1[0] + p2[0]);
   new_node->loc.x = a[0] + 0.5*a[1] + 0.25*a[2] + 0.125*a[3];

   a[0] = n1->loc.y;
   a[1] = p1[1];
   a[2] = 3. * (n2->loc.y - n1->loc.y) - (p2[1] + 2. * p1[1]);
   a[3] = 2. * (n1->loc.y - n2->loc.y) + (p1[1] + p2[1]);
   new_node->loc.y = a[0] + 0.5*a[1] + 0.25*a[2] + 0.125*a[3];

   a[0] = n1->loc.z;
   a[1] = p1[2];
   a[2] = 3. * (n2->loc.z - n1->loc.z) - (p2[2] + 2. * p1[2]);
   a[3] = 2. * (n1->loc.z - n2->loc.z) + (p1[2] + p2[2]);
   new_node->loc.z = a[0] + 0.5*a[1] + 0.25*a[2] + 0.125*a[3];

#ifdef CONN
   new_node->num_conn = 0;
   new_node->max_conn = 0;
#endif

   // add it to the list!
   new_node->index = (*node_cnt)++;
   new_node->next_node = node_head;
   node_head = new_node;

   return new_node;
}


/*
 * Create a node at the center of the three nodes provided.
 */
node_ptr create_center_point(int depth, node_ptr node1, node_ptr node2, node_ptr node3) {

   double d1,d2,d3,base_l;
   node_ptr new_node = (NODE *)malloc(sizeof(NODE));
   VEC r1, r2, r3;

   /* basis vector along one side */
   r1 = from(node1->loc,node2->loc);

   // base location for new node
   // new_node->loc = midpt(node1->loc,node2->loc);
   new_node->loc.x = (node1->loc.x+node2->loc.x+node3->loc.x)/3.0;
   new_node->loc.y = (node1->loc.y+node2->loc.y+node3->loc.y)/3.0;
   new_node->loc.z = (node1->loc.z+node2->loc.z+node3->loc.z)/3.0;

   /* basis vector from first node to 3rd node */
   r2 = from(node1->loc,node3->loc);

   // and do this one just for the hell of it
   r3 = from(node2->loc,node3->loc);

   /* base length is the average length of the two planar vectors (r1,r2) */
   base_l = (length(r1)+length(r2)+length(r3))/3.0;

   /* basis vector normal to the plane of the working triangle */
   r3 = norm(cross(r1,r2));

   // reorient r2 to be perpendicular to r1 and r3
   r2 = norm(cross(r1,r3));

   // shouldn't we scale d1 and d2 with the element size? Oh, depth does that.
   d1 = base_shake*pow(base_exponent,depth)*((1.0+rand())/RAND_MAX-0.5);
   d2 = base_l*base_shake*pow(base_exponent,depth)*((1.0+rand())/RAND_MAX-0.5);

   /* make this code use normal_bias in the same way as create_midpoint_4 */
   // now, aren't we taking size into account twice? Once with base_l and
   //   once with pow(normal_exponent,depth) ?
   d3 = base_l*normal_shake*pow(normal_exponent,depth)*((1.0+rand())/RAND_MAX-0.5+normal_bias);

   /* perturb the initial node placement in 3 basis directions */
   new_node->loc.x += d1*r1.x + d2*r2.x + d3*r3.x;
   new_node->loc.y += d1*r1.y + d2*r2.y + d3*r3.y;
   new_node->loc.z += d1*r1.z + d2*r2.z + d3*r3.z;

#ifdef CONN
   new_node->num_conn = 0;
   new_node->max_conn = 0;
#endif

   // add it to the list!
   new_node->index = 0;		// what do we do here?
   new_node->next_node = node_head;
   node_head = new_node;

   return new_node;
}
