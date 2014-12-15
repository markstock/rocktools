/*************************************************************
 *
 *  erodeutil.c - Utility subroutines for rocksmooth
 *
 *  Mark J. Stock, mstock@umich.edu
 *
 *
 * rocktools - Tools for creating and manipulating triangular meshes
 * Copyright (C) 1999,2006 Mark J. Stock
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
#include "structs.h"


/*
 * find_flow will determine the flow over all edges in the surface,
 * approximating rainfall as constant and ground losses as zero.
 */
int find_flow(tri_pointer tri_head) {

   int debug = 0;
   int i;
   //int last_solved_level;
   //int current_level;
   double min_z;
   double tri_area;
   double carried_flow;
   VEC rain_from;
   VEC tri_normal;
   node_ptr curr_node;
   node_ptr low_node;
   tri_pointer curr_tri;


   /* Set the rain incident vector */
   rain_from.x = 0.0;
   rain_from.y = 0.0;
   rain_from.z = 1.0;
   rain_from = norm(rain_from);


   /* Set all old flow rates to zero! */
   curr_node = node_head;
   while (curr_node) {
      curr_node->flow_rate = 0.0;
      curr_node = curr_node->next_node;
   }


   /* first determine the direct (first-order) flow, for each triangle
    * onto a specific node */
   curr_tri = tri_head;
   while (curr_tri) {

      /* find the lowest node, save it in curr_node */
      curr_node = curr_tri->node[0];
      min_z = curr_node->loc.z;
      for (i=1; i<3; i++) {
         if (curr_tri->node[i]->loc.z < min_z) {
            curr_node = curr_tri->node[i];
            min_z = curr_node->loc.z;
         }
      }
      if (debug > 0) fprintf(stderr,"lowest node in tri is at z= %lf",min_z);

      /* determine the area and the normal */
      tri_area = find_area(curr_tri);
      tri_normal = find_tri_normal(curr_tri);

      /* modify the total area based on the normal, to find the projected area */
      tri_area *= dot(tri_normal,rain_from);

      /* if the triangle faces down, no rain will land on it */
      if (tri_area < 0.0) tri_area = 0.0;

      /* add the first-order flow to that node */
      if (debug > 0) fprintf(stderr,", its flow is %lf",curr_node->flow_rate);
      curr_node->flow_rate += tri_area;
      if (debug > 0) fprintf(stderr," + %lf = %lf\n",tri_area,curr_node->flow_rate);

      curr_tri = curr_tri->next_tri;
   }

   /* determine the downhill node for each node, NULL if none */
   curr_node = node_head;
   while (curr_node) {

      if (curr_node->num_adj_nodes > 0) {
         if (debug > 0) fprintf(stderr,"looking for downstream node from %d adjacent nodes\n",curr_node->num_adj_nodes);
         min_z = curr_node->adj_node[0]->loc.z;
         low_node = (node_ptr) curr_node->adj_node[0];
         for (i=1; i<curr_node->num_adj_nodes; i++)
            if (curr_node->adj_node[i]->loc.z < min_z) {
               min_z = curr_node->adj_node[i]->loc.z;
               low_node = (node_ptr) curr_node->adj_node[i];
            }

         /* if the lowest adjacent node is higher than the current node,
          * then there is no downstream node, flow vanishes */
         if (min_z < curr_node->loc.z) {
            if (debug > 0) fprintf(stderr,"  curr_node is at z= %lf, lowest adj. node is at %lf\n",curr_node->loc.z,min_z);
            curr_node->downstream = (node_ptr) low_node;
         } else {
            if (debug > 0) fprintf(stderr,"  curr_node is at z= %lf, no downstream node\n",curr_node->loc.z);
            curr_node->downstream = (node_ptr) NULL;
         }

      } else {
         /* no adjacent nodes? something is wrong. */
         fprintf(stderr,"Node at x=%lf has no adjacent nodes set.\n",curr_node->loc.x);
         fprintf(stderr,"This is a problem. Exiting.\n");
         exit(0);
      }

      curr_node = curr_node->next_node;
   }


   /* now, the hard part, determine the flow at all nodes, assuming no
    * losses to groundwater */


   /* I am going to solve this the dumb, inefficient way for now. */

   /* for each node with direct flow, follow the list of downstream nodes
    * and add that amount of direct flow to each one. It's dumb, but it works. */
   /* but, first, save the amount of direct (level 0) flow into each node */
   curr_node = node_head;
   while (curr_node) {
      curr_node->temp_loc.x = curr_node->flow_rate;
      curr_node = curr_node->next_node;
   }

   curr_node = node_head;
   while (curr_node) {
      carried_flow = curr_node->temp_loc.x;	/* the original direct flow only */
      if (debug > 0) fprintf(stderr,"node has direct flow of %lf",carried_flow);
      if (carried_flow > 0.0) {
         if (debug > 0) fprintf(stderr,", flowing");
         low_node = (node_ptr) curr_node->downstream;
         while (low_node) {
            if (debug > 0) fprintf(stderr,".");
            low_node->flow_rate += carried_flow;
            low_node = (node_ptr) low_node->downstream;
         }
      }
      if (debug > 0) fprintf(stderr,"\n");
      curr_node = curr_node->next_node;
   }


/* *******************   CRAP   *********************************************
   last_solved_level = 0;
   current_level = 1;
   while () {

      // first, identify the nodes that will be solved for on this level
      curr_node = node_head;
      while (curr_node) {
         if (curr_node->index == last_solved_level) {
            curr_node->downstream->index = current_level;
         }
         curr_node = curr_node->next_node;
      }

      // then, go through all nodes at this level and solve for their flow
      curr_node = node_head;
      while (curr_node) {
         if (curr_node->index == current_level) {
            curr_node->downstream->index = current_level;
         }
         curr_node = curr_node->next_node;
      }


      // make sure last_solved_level doesn't get too high
      if (last_solved_level > 100) break;

      break;
   }
*****************************************   END CRAP   ******************* */


   // (void) write_flow_data();

   return 0;
}


/*
 * fill_basins tries to raise the level of land where there is
 * no downstream node
 */
int fill_basins(tri_pointer tri_head) {

   int i;
   double min_z;
   node_ptr curr_node;

   curr_node = node_head;
   while (curr_node) {
      if (!curr_node->downstream) {

         /* find the lowest neighbor node */
         min_z = curr_node->adj_node[0]->loc.z;
         for (i=1; i<curr_node->num_adj_nodes; i++)
            if (curr_node->adj_node[i]->loc.z < min_z)
               min_z = curr_node->adj_node[i]->loc.z;

         /* raise the current node to just above that lowest neighbor node */
         curr_node->loc.z = min_z+0.0005;

      }
      curr_node = curr_node->next_node;
   }

   return(1);
}


/*
 * Ideally, erode_surface takes the flow at all nodes and computes
 * any land levelling (+-z) that should take place. Fluvial erosion
 * depends not only on flow volume, but speed.
 */
int erode_surface(tri_pointer tri_head, double erosion_rate) {

   node_ptr curr_node;

   curr_node = node_head;
   while (curr_node) {

      // here's the erosion
      curr_node->loc.z -= erosion_rate*curr_node->flow_rate;

      curr_node = curr_node->next_node;
   }

   return(1);
}


/*
 * write_flow_data will write a debugging-like text description of the
 * solved flow
 */
int write_flow_data() {

   int num = 0;
   node_ptr curr_node = node_head;

   fprintf(stderr,"\nFlow data at time= 0.0\n");

   while (curr_node) {

      if (curr_node->downstream) {
         fprintf(stderr,"node at z= %g, downstream z= %g, total flow of %lf, and direct inflow of %lf\n",
            curr_node->loc.z,curr_node->downstream->loc.z,curr_node->flow_rate,curr_node->temp_loc.x);
      } else {
         fprintf(stderr,"node at z= %g, no downstream node, total flow of %lf, and direct inflow of %lf\n",
            curr_node->loc.z,curr_node->flow_rate,curr_node->temp_loc.x);
      }

      num++;
      curr_node = curr_node->next_node;
   }

   fprintf(stderr,"Wrote flow for %d nodes.\n\n",num);

   return num;
}

