/* **********************************************************
 *
 *  createutil.c - Subroutines used in rockcreate
 *
 *  Mark J. Stock, mstock@umich.edu
 *
 *
 * rocktools - Tools for creating and manipulating triangular meshes
 * Copyright (C) 1999,2009,14  Mark J. Stock
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

int read_files_for_nodes (int,char**);
VEC get_spherical_random_vector (double);
VEC get_random_gaussian_vector (double);
int create_cubic_nodes (int);
int create_gaussian_nodes (int);
int create_random_walk_nodes (int);
int sphericalize_nodes (int);
int pack_into_cube ();
tri_pointer create_convex_hull ();
int fix_vertex_order (tri_pointer);


/*
 * read_files_for_nodes - parse through command-line looking for files,
 * read those files in, but only keep the raw nodes that they contain
 */
int read_files_for_nodes (int argc, char **argv) {

   int i,num_nodes;
   char infile[256];
   tri_pointer tri_head,temp;
   //node_ptr this;

   // we begin with zero nodes
   num_nodes = 0;
   tri_head = NULL;

   // Parse command-line args
   for (i=1; i<argc; i++) {
      //fprintf(stderr,"option string (%s)\n",argv[i]);
      if (argv[i][0] != '-' && !isdigit(argv[i][0])) {

         // we might have a filename here
         (void) strcpy (infile,argv[i]);
         //fprintf(stderr,"  potential file name (%s)\n",infile);

         // call the appropriate file reader
         tri_head = read_input (infile,FALSE,NULL);

         // count the nodes
         num_nodes = count_nodes();
         //fprintf(stderr,"After reading %s, have %d nodes\n",infile,num_nodes); fflush(stderr);
      }
   }

   // free all memory associated with the triangles, but keep the nodes
   while (tri_head) {
      temp = tri_head->next_tri;
      free(tri_head);
      tri_head = temp;
   }
   //fprintf(stderr,"done removing triangles\n"); fflush(stderr);

   /*
   this = node_head;
   while (this) {
      fprintf(stderr,"  node at %g %g %g\n",this->loc.x,this->loc.y,this->loc.z);
      this = this->next_node;
   }
   */

   //exit(0);
   return (num_nodes);
}


/*
 * create_cubic_nodes - subroutine to create a number of random 3D points
 * in space inside of a unit cube
 */
int create_cubic_nodes (int num_nodes) {

   int i;
   int actual_nodes = 0;
   VEC loc;

   // create a field of randomly-spaced points in [-0.5:0.5]
   for (i=0; i<num_nodes; i++) {
      loc.x = rand()/(RAND_MAX+1.0) - 0.5;
      loc.y = rand()/(RAND_MAX+1.0) - 0.5;
      loc.z = rand()/(RAND_MAX+1.0) - 0.5;
      // add a new node, give NULL for tri pointer and bin pointer
      (void) add_to_nodes_list(NULL,&actual_nodes,-1,&loc,NULL);
   }

   return actual_nodes;
}


/*
 * create_gaussian_nodes - subroutine to create a number of random 3D points
 * in space inside of a unit cube
 */
int create_gaussian_nodes (int num_nodes) {

   int i;
   int actual_nodes = 0;
   VEC loc;

   // create a field of gaussian-randomly-spaced points around the origin
   for (i=0; i<num_nodes; i++) {
      loc = get_random_gaussian_vector(1.0);
      // add a new node, give NULL for tri pointer and bin pointer
      (void) add_to_nodes_list(NULL,&actual_nodes,-1,&loc,NULL);
   }

   return actual_nodes;
}


/*
 * create_random_walk_nodes - subroutine to create a number of random 3D points
 * in space inside of a unit cube
 */
int create_random_walk_nodes (int num_nodes) {

   int i;
   int actual_nodes = 0;
   double lengthper;
   VEC loc;

   // set length jump (diffusion dictates a length scale of time^0.5)
   lengthper = 1./sqrt((float)num_nodes);

   // place the first node
   if (num_nodes > 0) {
      loc.x = 0.0;
      loc.y = 0.0;
      loc.z = 0.0;
      // add a new node, give NULL for tri pointer and bin pointer
      (void) add_to_nodes_list(NULL,&actual_nodes,-1,&loc,NULL);
   }

   /* create a field of randomly-spaced points in [0:1][0:1][0:1] */
   for (i=1; i<num_nodes; i++) {
      loc = plus(loc,get_spherical_random_vector(lengthper));
      // add a new node, give NULL for tri pointer and bin pointer
      (void) add_to_nodes_list(NULL,&actual_nodes,-1,&loc,NULL);
   }

   return actual_nodes;
}


/*
 * return a random vector, evenly distributed over a sphere of radius thislen
 */
VEC get_spherical_random_vector (double thislen) {
   double len;
   VEC test;

   len = 1.;
   while (len > 0.25) {
      test.x = rand()/(RAND_MAX+1.0) - 0.5;
      test.y = rand()/(RAND_MAX+1.0) - 0.5;
      test.z = rand()/(RAND_MAX+1.0) - 0.5;
      len = test.x*test.x + test.y*test.y + test.z*test.z;
   }
   len = thislen/sqrt(len);
   test.x *= len;
   test.y *= len;
   test.z *= len;

   return test;
}


/*
 * return a random gaussian vector with given sigma
 */
VEC get_random_gaussian_vector (double sigma) {
   VEC test,out;

   test.x = (double)(rand())/(double)(RAND_MAX);
   test.y = (double)(rand())/(double)(RAND_MAX);
   out.x = (double)(sqrt(-2.*log(test.x))*cos(2.*M_PI*test.y));
   out.y = (double)(sqrt(-2.*log(test.x))*sin(2.*M_PI*test.y));
   test.x = (double)(rand())/(double)(RAND_MAX);
   test.y = (double)(rand())/(double)(RAND_MAX);
   out.z = (double)(sqrt(-2.*log(test.x))*cos(2.*M_PI*test.y));
   out.x *= sigma;
   out.y *= sigma;
   out.z *= sigma;

   return out;
}


/*
 * pack_into_cube - compresses the particle distribution into the unit cube
 */
int pack_into_cube () {

   node_ptr this;
   double scale;
   VEC min,max;

   /* set up the max and min vectors */
   min.x = 9.9e+9;
   min.y = 9.9e+9;
   min.z = 9.9e+9;
   max.x = -9.9e+9;
   max.y = -9.9e+9;
   max.z = -9.9e+9;

   /* find the min/max in all directions */
   this = node_head;
   while (this) {
      if (this->loc.x < min.x) min.x = this->loc.x;
      if (this->loc.x > max.x) max.x = this->loc.x;
      if (this->loc.y < min.y) min.y = this->loc.y;
      if (this->loc.y > max.y) max.y = this->loc.y;
      if (this->loc.z < min.z) min.z = this->loc.z;
      if (this->loc.z > max.z) max.z = this->loc.z;
      this = this->next_node;
   }

   /* find the scaling factor to use for all 3 directions */
   if ((max.x-min.x > max.y-min.y) && (max.x-min.x > max.z-min.z))
      scale = 1.0/(max.x-min.x);
   else if ((max.y-min.y > max.x-min.x) && (max.y-min.y > max.z-min.z))
      scale = 1.0/(max.y-min.y);
   else
      scale = 1.0/(max.z-min.z);

   /* scale all nodes to fill the unit box */
   this = node_head;
   while (this) {
      this->loc.x = (this->loc.x - min.x)*scale;
      this->loc.y = (this->loc.y - min.y)*scale;
      this->loc.z = (this->loc.z - min.z)*scale;
      //fprintf(stderr,"  node %d at %g %g %g\n",this->index,this->loc.x,this->loc.y,this->loc.z);
      this = this->next_node;
   }

   return (0);
}


/*
 * sphericalize_nodes - subroutine to perturb the nodes using gravitation-like laws to
 * reduce the number of small triangles in the convex hull
 */
int sphericalize_nodes (int num_cycles) {

   int i,num_nodes;
   double diff, dist;
   node_ptr this,target;
   VEC dx;

   /* now, run a simple gravitation-like repulsion routine on the points
    * add repulsion-attraction balance code */

   /* run through it num_cycles times */
   for (i=0; i<num_cycles; i++) {

      /* echo them */
      /* fprintf(stdout,"cycle %d\n",i);
      for (j=0; j<num_nodes; j++) {
         fprintf(stdout,"%lf %lf %lf\n",node[j].x,node[j].y,node[j].z);
      } */

      this = node_head;
      while (this) {
         this->temp_loc.x = 0.0;
         this->temp_loc.y = 0.0;
         this->temp_loc.z = 0.0;
         this = this->next_node;
      }

      /* find the repulsion influence between points target and this */
      num_nodes = 0;
      target = node_head;
      while (target) {
         this = node_head;
         while (this) {
            if (this != target) {
               /* first, find the distance */
               dx.x = target->loc.x - this->loc.x;
               dx.y = target->loc.y - this->loc.y;
               dx.z = target->loc.z - this->loc.z;
               dist = sqrt(pow(dx.x,2)+pow(dx.y,2)+pow(dx.z,2));
               /* then, compute the total influence */
               diff = 0.0075/dist - dist*0.01333;	/* repulsion-attraction, balance at 0.75 */
               // diff = 0.0001/dist;	/* repulsion only */
               target->temp_loc.x += diff*dx.x;
               target->temp_loc.y += diff*dx.y;
               target->temp_loc.z += diff*dx.z;
            }
            this = this->next_node;
         }
         num_nodes++;
         target = target->next_node;
      }

      /* now, apply those changes */
      diff = 1./(float)num_nodes;
      this = node_head;
      while (this) {
         this->loc.x += this->temp_loc.x*diff;
         this->loc.y += this->temp_loc.y*diff;
         this->loc.z += this->temp_loc.z*diff;
         this = this->next_node;
      }

   }

   return(0);
}


/*
 * create_convex_hull does just that, it uses the QuickHull routine
 * http://www.cse.unsw.edu.au/~lambert/java/3d/quickhull.html
 *
 * final_head is the final trimesh head, test_head will contain the list of
 * triangles that have not passed the convexity test
 */
tri_pointer create_convex_hull () {

   int num_tri = 0;
   int num_norms = 0;
   int i,num_nodes;
   double maxdist,thisdist;
   VEC testnorm;
   node_ptr thisn,farnode;
   tri_pointer test_head,new_tri,test,this,hole_head,prev,temp,new_head,final_head;
   temp = NULL;

   // the test list is the non-finalized list of triangles
   // nodes can be pointed to from tris in either final_head or test_head
   test_head = NULL;
   final_head = NULL;

   // prepare for binning normals
   NBIN normbin;
   (void) prepare_norm_bin (&normbin);

   // make sure that the incoming triangle list is empty (otherwise we'll get confused)
   if (final_head) {
      fprintf(stderr,"ERROR (create_convex_hull): final_head is not NULL\n");
      exit(1);
   }

   // make sure we have 4 or more nodes
   num_nodes = count_nodes();
   if (num_nodes < 4) {
      fprintf(stderr,"ERROR (create_convex_hull): num_nodes (%d) is < 4\n",num_nodes);
      exit(1);
   }

   // create a first guess at a convex hull (trivial)
   //   (pick 3 nodes and make two coincident triangles pointing opposite directions)

   new_tri = (TRI*)malloc(sizeof(TRI));
   new_tri->index = num_tri;
   num_tri++;
   new_tri->node[0] = node_head;
   new_tri->node[1] = new_tri->node[0]->next_node;
   new_tri->node[2] = new_tri->node[1]->next_node;
   // set the normal
   testnorm = find_tri_normal(new_tri);
   new_tri->norm[0] = add_to_norms_list(&num_norms,&testnorm,&normbin);
   // make it the new head
   new_tri->next_tri = NULL;
   test_head = new_tri;

   new_tri = (TRI*)malloc(sizeof(TRI));
   new_tri->index = num_tri;
   num_tri++;
   new_tri->node[0] = test_head->node[0];
   new_tri->node[1] = test_head->node[2];
   new_tri->node[2] = test_head->node[1];
   // set the normal
   testnorm = find_tri_normal(new_tri);
   new_tri->norm[0] = add_to_norms_list(&num_norms,&testnorm,&normbin);
   // set all adjacent pointers to the other tri
   new_tri->adjacent[0] = test_head;
   new_tri->adjacent[1] = test_head;
   new_tri->adjacent[2] = test_head;
   test_head->adjacent[0] = new_tri;
   test_head->adjacent[1] = new_tri;
   test_head->adjacent[2] = new_tri;
   // make it the new head
   new_tri->next_tri = test_head;
   test_head = new_tri;

   // iterate through test_head list
   test = test_head;
   while (test) {

      //printf("testing tri %d joining nodes %d %d %d\n",test->index,test->node[0]->index,test->node[1]->index,test->node[2]->index);

      // are there any nodes in the hemisphere along this triangle's normal?
      // loop through all nodes looking for ones far from the plane of the triangle
      maxdist = 1.e-12;
      farnode = NULL;
      thisn = node_head;
      while (thisn) {
         thisdist = dot(test->norm[0]->norm,from(test->node[0]->loc,thisn->loc));
         //printf("  dist %d is %g\n",thisn->index,thisdist);
         if (thisdist > maxdist) {
            maxdist = thisdist;
            farnode = thisn;
         }
         thisn = thisn->next_node;
      }

      // if the above filter catches one of the tri's nodes, call it a null
      if (farnode == test->node[0]) farnode = NULL;
      if (farnode == test->node[1]) farnode = NULL;
      if (farnode == test->node[2]) farnode = NULL;

      // is the triangle in the convex hull, or does the new node grow the test hull?
      if (farnode) {

         // a farnode exists, we need to replace triangles in our hull
         //printf("  node %d is the farnode\n",farnode->index);

         // make a list of triangles in the test list that this node can see
         // (naive test: just test all tris)
         hole_head = NULL;
         this = test_head;
         prev = NULL;
         while (this) {
            // is the farnode in the "outside" hemisphere of the this triangle?
            //if (dot(this->norm[0],from(test->node[0]->loc,farnode->loc)) > 0.0) {
            if (dot(this->norm[0]->norm,from(this->node[0]->loc,farnode->loc)) > 0.0) {
               //printf("    visible tri %d has nodes %d %d %d\n",this->index,this->node[0]->index,this->node[1]->index,this->node[2]->index);
               // save the next vis test
               temp = this->next_tri;
               // take this triangle out of the test list
               if (prev) prev->next_tri = temp;
               else test_head = temp;
               // make it the new head of the hole list
               this->next_tri = hole_head;
               hole_head = this;
               // and make sure that this points to the proper next tri
               this = temp;
               // don't update previous tri because it doesn't change
            } else {
               // update previous and advance to the next tri
               prev = this;
               this = this->next_tri;
            }
         }

         // test print
         /*
         printf("  test list now looks like\n");
         this = test_head;
         while (this) {
            printf("    tri %d with nodes %d %d %d adjacents %d %d %d\n",this->index,this->node[0]->index,this->node[1]->index,this->node[2]->index,this->adjacent[0]->index,this->adjacent[1]->index,this->adjacent[2]->index);
            this = this->next_tri;
         }
         printf("  hole list now looks like\n");
         this = hole_head;
         while (this) {
            printf("    tri %d with nodes %d %d %d adjacents %d %d %d\n",this->index,this->node[0]->index,this->node[1]->index,this->node[2]->index,this->adjacent[0]->index,this->adjacent[1]->index,this->adjacent[2]->index);
            this = this->next_tri;
         } */

         // find edges of hole list
         // first, set all tri indexes in the hole list to -1
         this = hole_head;
         while (this) {
            this->index = -1;
            this = this->next_tri;
         }
         // now, loop through all edges of all tris in hole list and look at the adjacent index
         this = hole_head;
         new_head = NULL;
         while (this) {
            for (i=0; i<3; i++) {
               // is this edge on the boundary of the hole?
               if (this->adjacent[i]->index != -1) {
                  // this edge and farnode node create a new triangle!
                  //printf("  edge %d of tri %d %d %d is a hole edge\n",i,this->node[0]->index,this->node[1]->index,this->node[2]->index);
                  new_tri = (TRI*)malloc(sizeof(TRI));
                  new_tri->index = num_tri;
                  num_tri++;
                  new_tri->node[0] = farnode;
                  new_tri->node[1] = this->node[i];
                  new_tri->node[2] = this->node[(i+1)%3];
                  //printf("    new tri has nodes %d %d %d\n",new_tri->node[0]->index,new_tri->node[1]->index,new_tri->node[2]->index);
                  // set the normal
                  testnorm = find_tri_normal(new_tri);
                  new_tri->norm[0] = add_to_norms_list(&num_norms,&testnorm,&normbin);
                  // set all adjacent pointers to the other tri
                  new_tri->adjacent[0] = NULL;
                  new_tri->adjacent[1] = this->adjacent[i];
                  new_tri->adjacent[2] = NULL;
                  // reciprocate the adjancency
                  if (this->adjacent[i]->node[0] == new_tri->node[2]) {
                     this->adjacent[i]->adjacent[0] = new_tri;
                  } else if (this->adjacent[i]->node[1] == new_tri->node[2]) {
                     this->adjacent[i]->adjacent[1] = new_tri;
                  } else {
                     this->adjacent[i]->adjacent[2] = new_tri;
                  }
                  // add it to the list of new tris
                  new_tri->next_tri = new_head;
                  new_head = new_tri;
               }
            }
            this = this->next_tri;
         }

         /*
         printf("  new list now looks like\n");
         this = new_head;
         while (this) {
            printf("    tri %d with nodes %d %d %d\n",this->index,this->node[0]->index,
                   this->node[1]->index,this->node[2]->index);
            this = this->next_tri;
         } */


         // remove all triangles in "hole" list
         this = hole_head;
         while (this) {
            temp = this->next_tri;
            free(this);
            this = temp;
         }

         // and connect the triangles in the "new" list (these are a fan, node 
         // "farnode" is the middle, and all tris radiate from it)
         this = new_head;
         while (this) {
            // is the first adjacent set?
            if (!this->adjacent[0]) {
               //printf("  looking for adjacent[0] of tri %d\n",this->index);
               // thisn is the lead node of the current triangle
               thisn = this->node[1];
               //printf("    want node %d\n",thisn->index);
               // look through all other tris for that node in position 2
               temp = this->next_tri;
               while (temp) {
                  //printf("    checking tri %d\n",temp->index);
                  if (temp->node[2] == thisn) {
                     this->adjacent[0] = temp;
                     // and reciprocate
                     temp->adjacent[2] = this;
                     // break out of the loop
                     temp = NULL;
                  } else {
                     temp = temp->next_tri;
                  }
               }
               // check results
               if (!this->adjacent[0]) {
                  fprintf(stderr,"ERROR (create_convex_hull): no adjacent 0 found!!\n");
                  exit(1);
               }
            }

            // recall the middle adjacent[1] is the hole triangle's adjacent,
            //   so it should be set already
            if (!this->adjacent[1]) {
               fprintf(stderr,"ERROR (create_convex_hull): no adjacent 1 found!!\n");
               exit(1);
            }

            // is the last adjacent set?
            if (!this->adjacent[2]) {
               //printf("  looking for adjacent[2] of tri %d\n",this->index);
               // thisn is the trailing node of the current triangle
               thisn = this->node[2];
               //printf("    want node %d\n",thisn->index);
               // look through all other tris for that node in position 1
               temp = this->next_tri;
               while (temp) {
                  //printf("    checking tri %d\n",temp->index);
                  if (temp->node[1] == thisn) {
                     this->adjacent[2] = temp;
                     // and reciprocate
                     temp->adjacent[0] = this;
                     // break out of the loop
                     temp = NULL;
                  } else {
                     temp = temp->next_tri;
                  }
               }
               // check results
               if (!this->adjacent[2]) {
                  fprintf(stderr,"ERROR (create_convex_hull): no adjacent 2 found!!\n");
                  exit(1);
               }
            }

            // go to the next tri in "new" list
            this = this->next_tri;
         }

         // finally, move all tris from "new" list to head of "test" list
         // find the last entry
         this = new_head;
         while (this) {
            // temp will be the last saved non-NULL value, it's the last tri in the list
            temp = this;
            this = this->next_tri;
         }
         // temp will now point to the next tri in the test list (we're done with "test" now)
         temp->next_tri = test_head;
         // and the head of the new list will become the new "test"
         test = new_head;
         test_head = test;

         /*
         printf("  final list now looks like\n");
         this = final_head;
         while (this) {
            printf("    tri %d with nodes %d %d %d adjacents %d %d %d\n",this->index,this->node[0]->index,this->node[1]->index,this->node[2]->index,this->adjacent[0]->index,this->adjacent[1]->index,this->adjacent[2]->index);
            this = this->next_tri;
         }

         printf("  test list now looks like\n");
         this = test_head;
         while (this) {
            printf("    tri %d with nodes %d %d %d adjacents %d %d %d\n",this->index,this->node[0]->index,this->node[1]->index,this->node[2]->index,this->adjacent[0]->index,this->adjacent[1]->index,this->adjacent[2]->index);
            this = this->next_tri;
         } */


      } else {

         // a farnode does not exist, this tri must be a part of the actual convex hull!
         //printf("  no farnode exists\n");
         // take it out of the test list
         temp = test->next_tri;
         // make it the new head of the final list
         test->next_tri = final_head;
         final_head = test;
         // and make sure that test points to the new head of the test list
         test = temp;
         // test and test_head should be the same
         if (final_head == test_head) test_head = temp;

         /*
         printf("  final list now looks like\n");
         this = final_head;
         while (this) {
            printf("    tri %d with nodes %d %d %d adjacents %d %d %d\n",this->index,this->node[0]->index,this->node[1]->index,this->node[2]->index,this->adjacent[0]->index,this->adjacent[1]->index,this->adjacent[2]->index);
            this = this->next_tri;
         }

         printf("  test list now looks like\n");
         this = test_head;
         while (this) {
            printf("    tri %d with nodes %d %d %d adjacents %d %d %d\n",this->index,this->node[0]->index,this->node[1]->index,this->node[2]->index,this->adjacent[0]->index,this->adjacent[1]->index,this->adjacent[2]->index);
            this = this->next_tri;
         } */

      }
   }

   // final test
   /*
   printf("  final list is\n");
   this = final_head;
   while (this) {
      printf("    tri %d with nodes %d %d %d\n",this->index,this->node[0]->index,
             this->node[1]->index,this->node[2]->index);
      this = this->next_tri;
   } */

   return (final_head);
}


/*
 * fix_vertex_order takes the tri_head list of triangles and makes
 * sure all of them have their vertexes ordered in the correct direction
 *
 * A triangle is considered correctly oriented if its midpoint[0]->num_conn
 * equals 1, and not if the value is 0
 */
int fix_vertex_order(tri_pointer tri_head) {

   //int curr_base;	/* the index of the adjacent triangle that has been checked/fixed */
   //int adj_base;	/* the node index of the shared vertex between current and adjacent tri */
   //int curr_test;	/* the index of the next (test) node of the current tri's shared side */
   //int adj_test;	/* the index of the test node from the adjacent triangle's list */
   int num_fixed = 0;
   VEC normal, center, test;
   node_ptr test_node, dummy_node;
   tri_pointer curr_tri, dummy_tri;

   /* check each triangle for vertex order */
   curr_tri = tri_head;
   while (curr_tri) {

      /* fprintf(stderr,"checking triangle with first node at (%lf)\n",curr_tri->node[0]->loc.x); */

      /* Before running set_adjacent, check all triangles's current normal; then
         check a point one triangle away from the main triangle. If the projection
         of the vector from that point to one of the points in the current triangle
         is negative, then the vertexes are oriented properly. This takes advantage 
         of the fact that the shape is convex. This would not work for non-convex hulls */

      normal = find_normal(curr_tri->node[0]->loc,curr_tri->node[1]->loc,curr_tri->node[2]->loc);
      center = midpt(curr_tri->node[0]->loc,midpt(curr_tri->node[1]->loc,curr_tri->node[2]->loc));

      /* fprintf(stderr,"normal is (%lf %lf %lf)\n",normal.x,normal.y,normal.z); */

      /* choose a test node, any one not in the curr_tri's plane */
      test_node = node_head;
      while (test_node == curr_tri->node[0] || test_node == curr_tri->node[1] || test_node == curr_tri->node[2]) {
         test_node = test_node->next_node;
      }
      test = test_node->loc;
      /* fprintf(stderr,"test is (%lf %lf %lf)\n",test.x,test.y,test.z);
      fprintf(stderr,"dot product is (%lf)\n",dot(from(center,test),normal)); */

      /* compute dot product of center-to-test onto normal */
      if ( dot(from(center,test),normal) > 0.0 ) {

         /* if this is positive, the normal is pointing toward the rest of
          * the rock, and the vertex order must be reversed */

         /* leave node[0] and adjacent[1], and swap the others */
         dummy_tri = curr_tri->adjacent[0];
         curr_tri->adjacent[0] = curr_tri->adjacent[2];
         curr_tri->adjacent[2] = dummy_tri;
         dummy_node = curr_tri->node[1];
         curr_tri->node[1] = curr_tri->node[2];
         curr_tri->node[2] = dummy_node;
         /* don't fiddle with node_record->conn_tri array, it doesn't matter here */
         /* yes, we do have to fiddle with it, because set_adjacent hasn't been run */
         num_fixed++;
      }
         /* otherwise, the normal is in the correct direction, and the
          * vertex order is correct */

      curr_tri = curr_tri->next_tri;
   }

   return num_fixed;
}
