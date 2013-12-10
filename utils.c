/*************************************************************
 *
 *  utils.c - Useful utility subroutines for rocktools
 *
 *  Mark J. Stock, mstock@umich.edu
 *
 *
 * rocktools - Tools for creating and manipulating triangular meshes
 * Copyright (C) 1999-2000,2002-2004,8  Mark J. Stock
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
#include <string.h>
#include <ctype.h>
//#include <malloc.h>
#include "structs.h"

int count_nodes();
int set_node_connectivity();
node_ptr add_to_nodes_list(tri_pointer,int*,int,VEC*,bin_ptr);
VEC vscale(double,VEC);
VEC from(VEC,VEC);
VEC plus(VEC,VEC);
double length(VEC);
double dot(VEC,VEC);
double find_area(tri_pointer);
VEC find_center(tri_pointer);
VEC find_tri_normal(tri_pointer);
VEC find_normal(VEC,VEC,VEC);
double find_tri_dist(tri_pointer,VEC);
int inside_bounds(double,double,double);


/*
 * Add a node to the list of nodes - and search for close nodes
 */
node_ptr add_to_nodes_list (tri_pointer the_tri, int* num_nodes, int index, VEC* location, bin_ptr thebin) {

   node_ptr curr_node = NULL;
   int ibin = 0;
   int found_match = FALSE;
   double thisx;
   double match_thresh = 1.e-5;     /* threshhold to match node locations */

   // first, search the list for a node close to this
   // new way to search
   if (thebin) {
      // start searching only in bin ibin
      if (thebin->axis == 0) thisx = (*location).x;
      else if (thebin->axis == 1) thisx = (*location).y;
      else thisx = (*location).z;
      ibin = (int)((thisx-thebin->start)/thebin->dx);

      // fprintf(stderr,"  search in bin %d\n",ibin); fflush(stderr);
      curr_node = thebin->b[ibin];
   } else {
      // search through all nodes, starting with the head
      curr_node = node_head;
   }
   while (curr_node) {
      //fprintf(stderr,"  does location (%g %g) match node at (%g %g)? \n",(*location).x,(*location).y,curr_node->loc.x,curr_node->loc.y); fflush(stderr);
      if (fabs(curr_node->loc.x - (*location).x) < match_thresh) {
         if (fabs(curr_node->loc.y - (*location).y) < match_thresh) {
            if (fabs(curr_node->loc.z - (*location).z) < match_thresh) {
               found_match = TRUE;
               // fprintf(stderr,"yes!\n");
               break;
            }
         } else {
            // fprintf(stderr,"no.\n");
         }
      } else {
         // fprintf(stderr,"no.\n");
      }
      if (thebin) curr_node = curr_node->next_bnode;
      else curr_node = curr_node->next_node;
   }


   /* did we find a match in the list of existing nodes? */
   if (found_match) {

      /* add some data to the specific node entry */
#ifdef CONN
      if (the_tri) add_conn_tri (curr_node, the_tri, index);
#endif

   } else {

      /* if not, create one and add it to the list */
      curr_node = (NODE *)malloc(sizeof(NODE));
      curr_node->index = (*num_nodes)++;
      curr_node->loc.x = (*location).x;
      curr_node->loc.y = (*location).y;
      curr_node->loc.z = (*location).z;
#ifdef CONN
      curr_node->num_conn = 0;
      curr_node->max_conn = 0;
      if (the_tri) add_conn_tri (curr_node, the_tri, index);
      //curr_node->conn_tri = (tri_pointer*)malloc(curr_node->max_conn*sizeof(tri_pointer));
      //curr_node->conn_tri_node = (int*)malloc(curr_node->max_conn*sizeof(int));
      //curr_node->conn_tri[0] = the_tri;
      //curr_node->conn_tri_node[0] = index;
      //curr_node->num_conn = 1;
#endif
      // add it to the head of the full list
      curr_node->next_node = node_head;
      node_head = curr_node;
      // add it to the head of the bin's list
      if (thebin) {
         curr_node->next_bnode = thebin->b[ibin];
         thebin->b[ibin] = curr_node;
      }
      // fprintf(stderr,"  adding new node at %g %g %g, num_conn= 1\n",location->x,location->y,location->z); fflush(stderr);
   }

   return curr_node;
}


#ifdef CONN
/*
 * Add a triangle to a node's connectivity lists
 */
int add_conn_tri (node_ptr curr_node, tri_pointer the_tri, int index) {

   int i;
   tri_pointer *new_conn_tri;
   int *new_conn_tri_node;

   if (curr_node->num_conn == curr_node->max_conn || curr_node->max_conn==0) {
      // malloc more room in the arrays
      if (curr_node->max_conn == 0) curr_node->max_conn = 1;
      else curr_node->max_conn *= 2;

      // extend the conn_tri array
      new_conn_tri = (tri_pointer*)malloc(curr_node->max_conn*sizeof(tri_pointer));
      for (i=0; i<curr_node->num_conn; i++)
         new_conn_tri[i] = curr_node->conn_tri[i];
      free(curr_node->conn_tri);
      curr_node->conn_tri = new_conn_tri;

      // extend the conn_tri_node array
      new_conn_tri_node = (int*)malloc(curr_node->max_conn*sizeof(int));
      for (i=0; i<curr_node->num_conn; i++)
         new_conn_tri_node[i] = curr_node->conn_tri_node[i];
      free(curr_node->conn_tri_node);
      curr_node->conn_tri_node = new_conn_tri_node;

      //fprintf(stderr,"node %d now has array length %d\n",curr_node->index,curr_node->max_conn);
   }

   // actually link the tri to the node
   curr_node->conn_tri[curr_node->num_conn] = the_tri;
   curr_node->conn_tri_node[curr_node->num_conn] = index;
   curr_node->num_conn++;

   //fprintf(stderr,"  add tri to existing node (%g %g %g), num_conn = %d\n",
   //   (*location).x,(*location).y,(*location).z,curr_node->num_conn);
   //fprintf(stderr,"node %d has conn_tri %d and c_t_node %d\n",
   //   curr_node->index,curr_node->conn_tri[curr_node->num_conn-1]->index,
   //   curr_node->conn_tri_node[curr_node->num_conn-1]);
   //fflush(stderr);

   return (curr_node->max_conn);
}
#endif


#ifdef CONN
/*
 * Determine the adjacent triangles for each triangle
 */
int set_adjacent_tris (tri_pointer tri_head) {

   int i,j;
   int num_set = 0;
   int test_node_index;
   tri_pointer this_tri = tri_head;
   tri_pointer test_tri;	// the triangle we're testing for test_node
   node_ptr test_node;		// the node pointer we're looking for
   node_ptr test_node2;

   while (this_tri) {

      //fprintf(stderr,"set adjacent to tri %d\n",this_tri->index);

      // loop through all triangles with test_node2, look for one with test_node also
      for (j=0; j<3; j++) {
         if (!this_tri->adjacent[j]) {
            test_node = this_tri->node[j];
            test_node2 = this_tri->node[(j+1)%3];
            //fprintf(stderr," find adjacent on side %d\n",j);
            //fprintf(stderr,"  loop thru tris with node %d looking for node %d\n",
            //   test_node2->index,test_node->index);
            for (i=0; i<test_node2->num_conn; i++) {
               test_tri = test_node2->conn_tri[i];
               if (test_tri) {
                  //fprintf(stderr,"  looking at tri %d, node %d\n",test_tri->index,test_node2->conn_tri_node[i]);
                  // if test_tri and this_tri are oriented the same way, this will work
                  // just test conn_tri_node[i]+1
                  test_node_index = mod(test_node2->conn_tri_node[i]+1,3);
                  //fprintf(stderr,"    node %d\n",test_tri->node[test_node_index]->index);
                  if (test_tri->node[test_node_index] == test_node) {
                     this_tri->adjacent[j] = test_tri;
                     test_tri->adjacent[test_node2->conn_tri_node[i]] = this_tri;
                     //fprintf(stderr,"  two tris (%d,%d) share side (%d) %g %g %g to %g %g %g\n",
                     //        this_tri->index,test_tri->index,(j+1)%3,
                     //        test_node2->loc.x,test_node2->loc.y,test_node2->loc.z,
                     //        test_tri->node[test_node_index]->loc.x,test_tri->node[test_node_index]->loc.y,test_tri->node[test_node_index]->loc.z);
                     break;
                  }
                  // if test_tri and this_tri are oriented oppositely, this will work
                  // also test conn_tri_node[i]+2
                  //test_node_index = test_node2->conn_tri_node[i];
                  //test_node_index = mod(test_node2->conn_tri_node[i]+2,3);
                  //fprintf(stderr,"    node %d\n",test_tri->node[test_node_index]->index);
                  //if (test_tri->node[test_node_index] == test_node) {
                  //   this_tri->adjacent[j] = test_tri;
                  //   test_tri->adjacent[test_node2->conn_tri_node[i]] = this_tri;
                  //   fprintf(stderr,"  two tris (%d,%d) share side (%d) %g %g %g to %g %g %g\n",
                  //           this_tri->index,test_tri->index,(j+1)%3,
                  //           test_node2->loc.x,test_node2->loc.y,test_node2->loc.z,
                  //           test_tri->node[test_node_index]->loc.x,test_tri->node[test_node_index]->loc.y,test_tri->node[test_node_index]->loc.z);
                  //   //i = 20;		// found it, skip out of loop
                  //   break;
                  //}
               } else {
                  break;
               }
            }
         }
      }

      this_tri = this_tri->next_tri;
   }

   return num_set;
}


/*
 * Make sure all triangles are oriented similarly, must be run before
 * set_adjacent_tris, as it doesn't fix that connectivity information.
 */
int fix_orientation (tri_pointer tri_head) {

   int i,j,k,sum;
   int flipdir,num_flipped = 0;
   int ip1,ip2,numinarray;
   int inarray[20],dir[20];
   tri_pointer this_tri;
   tri_pointer adjtri[20];	// the triangle we're testing for test_node
   node_ptr this_node;		// the node pointer we're looking for
   node_ptr adjnode[20];


   this_node = node_head;
   while (this_node) {

      // fprintf(stderr,"Looking at node %d\n",this_node->index);
      for (j=0; j<20; j++) inarray[j] = FALSE;

      // set tri 0, as nodes 0 and 1
      inarray[0] = TRUE;
      adjtri[0] = this_node->conn_tri[0];
      i = this_node->conn_tri_node[0];
      adjnode[0] = adjtri[0]->node[(i+1)%3];
      adjnode[1] = adjtri[0]->node[(i+2)%3];
      dir[0] = 0;		// define direction as OK==0
      numinarray = 1;

      // add on the rest of the triangles, one at a time
      for (k=1; k<this_node->num_conn; k++) {
      // loop through all the other triangles, looking for the one next to this
      for (j=0; j<this_node->num_conn; j++) if (!inarray[j]) {
         i = this_node->conn_tri_node[j];
         this_tri = this_node->conn_tri[j];
         ip1 = (i+1)%3;
         ip2 = (i+2)%3;
         if (this_tri->node[ip1] == adjnode[numinarray]) {
            // this tri is oriented the same as 0
            inarray[j] = TRUE;
            adjtri[numinarray] = this_tri;
            adjnode[numinarray+1] = this_tri->node[ip2];
            dir[numinarray] = 0;
            numinarray++;
         } else if (this_tri->node[ip2] == adjnode[numinarray]) {
            // this tri is oriented opposite 0
            inarray[j] = TRUE;
            adjtri[numinarray] = this_tri;
            adjnode[numinarray+1] = this_tri->node[ip1];
            dir[numinarray] = 1;
            numinarray++;
         }
      }
      }

//    fprintf(stderr,"out of %d connecting tris:\n",this_node->num_conn);
//    for (k=0; k<this_node->num_conn; k++) {
//       fprintf(stderr,"elem %d, nodes %d %d, dir %d\n",adjtri[k]->index,adjnode[k]->index,adjnode[k+1]->index,dir[k]);
//       fprintf(stderr,"   nodes in order: %d %d %d\n",adjtri[k]->node[0]->index,adjtri[k]->node[1]->index,adjtri[k]->node[2]->index);
//    }

      // now, see which elems are oriented strangely
      sum = 0;
      for (k=0; k<this_node->num_conn; k++) sum += dir[k];
      if ((double)(sum)/(double)(this_node->num_conn) > 0.5) {
         // flip all tris whose direction is 0
         flipdir = 0;
      } else {
         flipdir = 1;
      }

      // flip the tris in question
      for (k=0; k<this_node->num_conn; k++) if (dir[k]==flipdir) {
         flip_tri(adjtri[k]);
         num_flipped++;
      }

      this_node = this_node->next_node;
   }

   return(num_flipped);
}


/*
 * Flip the normal of the element|triangle
 */
int flip_tri(tri_pointer this) {

   int j;
   node_ptr temp;

   // fprintf(stderr,"flipping elem %d\n",this->index);

   // leave local node 0 the same, flip nodes 1 and 2
   for (j=0; j<this->node[1]->num_conn; j++) {
      if (this->node[1]->conn_tri[j] == this) {
         this->node[1]->conn_tri_node[j] = 2;
         break;
      }
   }
   for (j=0; j<this->node[2]->num_conn; j++) {
      if (this->node[2]->conn_tri[j] == this) {
         this->node[2]->conn_tri_node[j] = 1;
         break;
      }
   }

   temp = this->node[1];
   this->node[1] = this->node[2];
   this->node[2] = temp;

   return(0);
}
#endif


/*
 * Write out the list of nodes
 */
int write_node_list() {

   int num_nodes = 0;
   node_ptr curr_node = node_head;

   while (curr_node) {
      fprintf(stderr,"Node %d is at %lf %lf %lf\n",num_nodes,curr_node->loc.x,curr_node->loc.y,curr_node->loc.z);
      num_nodes++;
      curr_node = curr_node->next_node;
   }

   return num_nodes;
}


/*
 * Count the number of nodes
 */
int count_nodes() {

   int num_nodes = 0;
   node_ptr curr_node = node_head;

   while (curr_node) {
      num_nodes++;
      curr_node = curr_node->next_node;
   }

   return num_nodes;
}


/*
 * Write out the list of tris
 */
int write_tri_list (tri_pointer tri_head) {

   int num_tris = 0;
   tri_pointer this = tri_head;

   while (this) {
      fprintf(stderr,"Tri %d with node %d %d %d\n",this->index,
              this->node[0]->index,this->node[1]->index,this->node[2]->index);
      num_tris++;
      this = this->next_tri;
   }

   return num_tris;
}


/*
 * Return the area of the triangle in question
 */
double find_area(tri_pointer thetri) {

   double a,b,c,s,area;

   a = length(from(thetri->node[0]->loc,thetri->node[1]->loc));
   b = length(from(thetri->node[2]->loc,thetri->node[1]->loc));
   c = length(from(thetri->node[0]->loc,thetri->node[2]->loc));
   s = 0.5*(a+b+c);
   area = sqrt(s*(s-a)*(s-b)*(s-c));

   return area;
}


/*
 * Return the center of the triangle in question
 */
VEC find_center(tri_pointer thetri) {

   VEC center;

   center.x = (thetri->node[0]->loc.x+thetri->node[1]->loc.x+thetri->node[2]->loc.x) / 3.;
   center.y = (thetri->node[0]->loc.y+thetri->node[1]->loc.y+thetri->node[2]->loc.y) / 3.;
   center.z = (thetri->node[0]->loc.z+thetri->node[1]->loc.z+thetri->node[2]->loc.z) / 3.;

   return center;
}


/*
 * Return the distance from a triangle's center to a point
 */
double find_tri_dist(tri_pointer thetri,VEC point) {

   //int i;
   double dist;
   VEC center = find_center(thetri);

   //center.x = thetri->node[0]->loc.x+thetri->node[1]->loc.x+thetri->node[2]->loc.x;
   //center.y = thetri->node[0]->loc.y+thetri->node[1]->loc.y+thetri->node[2]->loc.y;
   //center.z = thetri->node[0]->loc.z+thetri->node[1]->loc.z+thetri->node[2]->loc.z;
   dist = length(from(center,point));

   return dist;
}


/*
 * Update the min/max from an array of doubles
 */
int update_minmax(double *this, VEC *min, VEC *max) {

   if (this[0] > max->x) max->x = this[0];
   if (this[1] > max->y) max->y = this[1];
   if (this[2] > max->z) max->z = this[2];
   if (this[0] < min->x) min->x = this[0];
   if (this[1] < min->y) min->y = this[1];
   if (this[2] < min->z) min->z = this[2];

   return 0;
}


/*
 * Return the cross product, v1xv2
 */
VEC cross(VEC v1, VEC v2) {

   VEC cp;

   cp.x = v1.y*v2.z - v1.z*v2.y;
   cp.y = v1.z*v2.x - v1.x*v2.z;
   cp.z = v1.x*v2.y - v1.y*v2.x;

   return cp;
}


/*
 * Return the scaled vector
 */
VEC vscale(double a, VEC p1) {

   VEC p2;

   p2.x = a*p1.x;
   p2.y = a*p1.y;
   p2.z = a*p1.z;

   return p2;
}


/*
 * Return the vector pointing from p1 to p2
 */
VEC from(VEC p1, VEC p2) {

   p2.x -= p1.x;
   p2.y -= p1.y;
   p2.z -= p1.z;

   return p2;
}


/*
 * Return the vector sum of p1 and p2
 */
VEC plus(VEC p1, VEC p2) {

   p2.x += p1.x;
   p2.y += p1.y;
   p2.z += p1.z;

   return p2;
}


/*
 * Return the midpoint of the two points
 */
VEC midpt(VEC p1, VEC p2) {

   p1.x = (p1.x+p2.x)/2.;
   p1.y = (p1.y+p2.y)/2.;
   p1.z = (p1.z+p2.z)/2.;

   return p1;
}


/*
 * Return the normalization of the vector to unit length
 */
VEC norm(VEC v1) {

   double dr = 1./sqrt(pow(v1.x,2)+pow(v1.y,2)+pow(v1.z,2));

   v1.x *= dr;
   v1.y *= dr;
   v1.z *= dr;

   return v1;
}


/*
 * Returns the length of the vector
 */
double length(VEC v1) {

   double length = sqrt(pow(v1.x,2)+pow(v1.y,2)+pow(v1.z,2));

   return length;
}


/*
 * Returns the length squared of the vector
 */
double lengthsq(VEC v1) {

   double lengthsq = pow(v1.x,2)+pow(v1.y,2)+pow(v1.z,2);

   return lengthsq;
}


/*
 * Return the normal of a defined triangle
 */
VEC find_tri_normal(tri_pointer the_tri) {

   return (find_normal(the_tri->node[0]->loc,the_tri->node[1]->loc,the_tri->node[2]->loc));
}


/*
 * Return the normal from three points
 */
VEC find_normal(VEC pt1,VEC pt2,VEC pt3) {

   double length;
   /* fprintf(stderr,"finding normal, first point is %g %g %g\n",pt1.x,pt1.y,pt1.z); */
   /* fprintf(stderr,"finding normal, first point is %g %g %g\n",pt2.x,pt2.y,pt2.z); */
   /* fprintf(stderr,"finding normal, first point is %g %g %g\n",pt3.x,pt3.y,pt3.z); */

   /* pt2 = subtract(pt2,pt1); */
   pt2.x -= pt1.x;
   pt2.y -= pt1.y;
   pt2.z -= pt1.z;
   /* pt3 = subtract(pt3,pt1); */
   pt3.x -= pt1.x;
   pt3.y -= pt1.y;
   pt3.z -= pt1.z;
   /* pt1 = cross(pt2,pt3); */
   pt1.x = pt2.y*pt3.z - pt2.z*pt3.y;
   pt1.y = pt2.z*pt3.x - pt2.x*pt3.z;
   pt1.z = pt2.x*pt3.y - pt2.y*pt3.x;
   /* fprintf(stderr,"                cross product is %g %g %g\n",pt1.x,pt1.y,pt1.z); */
   /* pt1 = normalize(pt1); */
   length = 1./sqrt((pt1.x*pt1.x) + (pt1.y*pt1.y) + (pt1.z*pt1.z));
   /* fprintf(stderr,"                length is %lf\n",length); */
   pt1.x *= length;
   pt1.y *= length;
   pt1.z *= length;

   return pt1;
}


/*
 * Returns the dot product of the two vectors
 */
double dot(VEC v1,VEC v2) {

   double dot_prod = v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;

   return dot_prod;
}


/*
 * Returns the vector projection of v1 onto v2
 *
 * Incomplete ****
 */
VEC projection(VEC v1,VEC v2) {

   /* return vector must be in direction of v2 */
   fprintf(stderr,"ERROR (projection): this routine should never be called!\n");
   exit(1);

   return v1;
}


/*
 * Returns z-angle of line from pt1 to pt2 in 2 dimensions?
 * angle is in radians, and 0=2pi=positive x axis
 */
double theta(VEC pt1,VEC pt2) {

   double angle;

   pt2.x -= pt1.x;
   pt2.y -= pt1.y;

   if (pt2.x > 0) {
      angle = atan2(pt2.y,pt2.x);
      /* fprintf(stderr,"atan2(%g,%g) = %g\n",pt2.y,pt2.x,angle); */
      if (angle < 0) angle += 2*M_PI;		/* M_PI is pi */
   } else {
      angle = M_PI - atan2(pt2.y,-1.0*pt2.x);
      /* fprintf(stderr,"atan2(%g,%g) = %g\n",pt2.y,-1.0*pt2.x,atan2(pt2.y,-1.0*pt2.x)); */
   }

   return angle;
}


/*
 * prepare the node_bin structure, given bounds
 */
void prepare_node_bin (bin_ptr bin, VEC nmin, VEC nmax) {
   int i;

   // split on longest axis
   if (nmax.x-nmin.x >= nmax.y-nmin.y && nmax.x-nmin.x >= nmax.z-nmin.z) {
      bin->axis = 0;
      bin->dx = (nmax.x-nmin.x)/(BIN_COUNT-1);
      bin->start = nmin.x - 0.5*bin->dx;
   } else if (nmax.y-nmin.y >= nmax.x-nmin.x && nmax.y-nmin.y >= nmax.z-nmin.z) {
      bin->axis = 1;
      bin->dx = (nmax.y-nmin.y)/(BIN_COUNT-1);
      bin->start = nmin.y - 0.5*bin->dx;
   } else {
      bin->axis = 2;
      bin->dx = (nmax.z-nmin.z)/(BIN_COUNT-1);
      bin->start = nmin.z - 0.5*bin->dx;
   }
   for (i=0;i<BIN_COUNT;i++) bin->b[i] = NULL;
   return;
}


/*
 * allocate memory for a two-dimensional array of floats
 */
float** allocate_2d_array_f(int nx,int ny) {

   int i;
   float **array = (float **)malloc(nx * sizeof(float *));

   array[0] = (float *)malloc(nx * ny * sizeof(float));
   for (i=1; i<nx; i++)
      array[i] = array[0] + i * ny;

   return(array);
}

int free_2d_array_f(float** array){
   free(array[0]);
   free(array);
   return(0);
}


/*
 * allocate memory for a two-dimensional array of bytes
 */
unsigned char** allocate_2d_array_uc(int nx,int ny) {

   int i;
   unsigned char **array = (unsigned char **)malloc(nx * sizeof(unsigned char *));

   array[0] = (unsigned char *)malloc(nx * ny * sizeof(unsigned char));
   for (i=1; i<nx; i++)
      array[i] = array[0] + i * ny;

   return(array);
}

int free_2d_array_uc(unsigned char** array){
   free(array[0]);
   free(array);
   return(0);
}


#ifdef ADJ_NODE
#ifdef CONN
/*
 * Determine the node connectivity, store it in a linked list
 * of node_record objects
 */
int set_node_connectivity() {

   int i,j,index;
   node_ptr curr_node = node_head;
   node_ptr test_node;

   fprintf(stderr,"Setting up node connectivity data...\n");
   fflush(stderr);

   /* Create the temp_node_record, a temporary storage list */
   while (curr_node) {
      curr_node->num_adj_nodes = 0;

      /* fprintf(stderr,"looking at node at x=%g, has adjacent nodes:\n",curr_node->loc.x); */

      /* loop over all triangles sharing this node */
      for (i=0; i<curr_node->num_conn; i++) {
         index = mod(curr_node->conn_tri_node[i]+1,3);
         test_node = curr_node->conn_tri[i]->node[index];

         /* check the node's current list of adjacent nodes for this new node */
         /* fprintf(stderr,"  looking through %d nodes\n",curr_node->num_adj_nodes); */
         for (j=0; j<curr_node->num_adj_nodes; j++)
            if (test_node == curr_node->adj_node[j]) {
               /* fprintf(stderr,"   matched at j=%d of %d\n",j,curr_node->num_adj_nodes); */
               break;
            }

         /* if it wasn't in the list, add it */
         if (j == curr_node->num_adj_nodes) {
            curr_node->adj_node[curr_node->num_adj_nodes] = test_node;
            curr_node->num_adj_nodes++;
            /* fprintf(stderr,"   added node at x=%g\n",curr_node->adj_node[j]->loc.x); */
         }

         /* and do it all over again for the other node in this tri */
         index = mod(curr_node->conn_tri_node[i]+2,3);
         test_node = curr_node->conn_tri[i]->node[index];
         /* fprintf(stderr,"  looking through %d nodes\n",curr_node->num_adj_nodes); */
         for (j=0; j<curr_node->num_adj_nodes; j++)
            if (test_node == curr_node->adj_node[j]) {
               /* fprintf(stderr,"   matched at j=%d of %d\n",j,curr_node->num_adj_nodes); */
               break;
            }
         if (j == curr_node->num_adj_nodes) {
            curr_node->adj_node[curr_node->num_adj_nodes] = test_node;
            curr_node->num_adj_nodes++;
            /* fprintf(stderr,"   added node at x=%g\n",curr_node->adj_node[j]->loc.x); */
         }
      }

      /* now, just write some test data */
      for (j=0; j<curr_node->num_adj_nodes; j++) {
         /* fprintf(stderr,"   at x=%g\n",curr_node->adj_node[j]->loc.x); */
      }

      curr_node = curr_node->next_node;
   }

   /* if we got this far, it worked, return TRUE */
   return (TRUE);
}
#endif
#endif


/*
 * inside_bounds tells whether the specific value is within the
 * two bounds specified
 */
int inside_bounds(double test_value, double low_bound, double high_bound) {

   if ((test_value > high_bound) || (test_value < low_bound)) return (FALSE);
   return (TRUE);
}

