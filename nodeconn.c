/*************************************************************
 *
 *  nodeconn.c - Useful utility subroutines for rocktools
 *
 *  Mark J. Stock, mstock@umich.edu
 *
 *
 * rocktools - Tools for creating and manipulating triangular meshes
 * Copyright (C) 1999-2004  Mark J. Stock
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
#include <malloc.h>
#define ADJ_NODE
#include "structs.h"

int set_node_connectivity();


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

