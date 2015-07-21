/* **********************************************************
 *
 *  convexhull.c - Subroutine for finding the convex hull of points
 *
 *  Mark J. Stock, mstock@umich.edu
 *
 *
 * rocktools - Tools for creating and manipulating triangular meshes
 * Copyright (C) 1999,2009,14-15  Mark J. Stock
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
//#include <string.h>
//#include <math.h>
//#include <ctype.h>
#include "structs.h"

tri_pointer create_convex_hull ();


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

   new_tri = alloc_new_tri();
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

   new_tri = alloc_new_tri();
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
                  new_tri = alloc_new_tri();
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

