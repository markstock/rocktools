/***********************************************************
 *
 *  structs.h - data structures supporting rocktools
 *
 *  Mark J. Stock, mstock@umich.edu
 *
 *
 * rocktools - Unix tools for creating and manipulating triangular meshes
 * Copyright (C) 1999,2003,2004,2006  Mark J. Stock
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
 ***********************************************************/


/*
 * Here are some structs and pointers to hold the triangle data
 */

#define TRUE 1
#define FALSE 0
#define mod(a,b) ((a)-(b)*floor((a)/(b)))

#define DOTPER 5000
#define DPMO (DOTPER-1)

//#define MAX_NODES 200000
#define MAX_IMAGE 32768
#define BIN_COUNT 10000
/*#define MAX_CONN 20*/
#define MAX_ADJ 10

#define M_PI           3.14159265358979323846

/*
 * A basic 3D vector
 */
typedef struct three_d_vector {
   double x;
   double y;
   double z;
} VEC;

/*
 * Pointers to the two types of data structures below
 */
typedef struct node_record *node_ptr;
typedef struct tri_record *tri_pointer;


/*
 * A data structure definition of a node, shared by an
 * arbitrary number of triangles
 */
typedef struct node_record {
   int index;			/* index of the node */
   VEC loc;			/* node's location in space */

#ifdef CONN
   int num_conn;		/* number of triangles using this node */
   int max_conn;		/* length of conn_* arrays */
   tri_pointer *conn_tri;	/* pointers to triangles using it */
   int *conn_tri_node;		/* index (0,1,2) of this node in the respective connecting tri */
#endif

#ifdef ADJ_NODE
   /* rocksmooth needs some more data for each node, this includes: */
   VEC temp_loc;		/* temporary location, need to keep separate */
   int num_adj_nodes;		/* number of nodes adjacent to this node */
   node_ptr adj_node[MAX_ADJ];	/* pointers to all adjacent nodes */
#endif

#ifdef ERODE
   /* rockerode needs some more data for each node, this includes: */
   node_ptr downstream;		/* pointer to the one downstream node */
   double flow_rate;		/* mass flow rate, used for erosion */
#endif

   node_ptr next_bnode;		// next node in the bin
   node_ptr next_node;
} NODE;


/*
 * A data structure definition of a single triangle
 */
typedef struct tri_record {
   int index;			/* index of the node */
   node_ptr node[3];		/* pointer to the three nodes, CCW */
   int use_norm;		/* use the normals? 0=no, 1=yes, alt usage in rockdetail */
   VEC norm[3];			/* surface normal vectors at the nodes */
   tri_pointer adjacent[3];	/* pointers to the three adjacent triangles */
 				/* adjacent[0] refers to edge between node[0] and node[1] */

   // midpoint used in detail, inout
   node_ptr midpoint[3];	/* midpoints inherited from adj tris' splitting */
 				/* midpoint[0] refers to the midpoint location */
#ifdef DETAIL
   int splittable;		// is this tri splittable?
#endif

   tri_pointer next_tri;	/* pointer to the next element in the list */
} TRI;


/*
 * structure for a binning system for nodes
 */
typedef struct bin_record *bin_ptr;
typedef struct bin_record {
   int axis;			// along which axis does the binning run? 0=x,1=y,2=z
   double start,dx;		// starting point and bin size
   node_ptr b[BIN_COUNT];	// the bins
} BIN;


/*
 * structure for marker characteristics (rockmarker)
 */
typedef struct marker_record {

   enum marker_type_record {
      sphere,
      rectangle,
      dualcone
   } marker_type;

   double marker_size;
   int randomize_size;
   double min_size;
   double max_size;
   double size_range;

   enum density_type_record {
      by_area,
      by_elem
   } density_type;
   double density;

   int randomize_radius;
   double min_radius;
   double max_radius;
   double radius_range;
   int randomize_height;
   double min_height;
   double max_height;
   double height_range;
   int randomize_rotation;
   int randomize_normal;
   double normal_pert;
} MARKER;

/*
 * All extern declarations used by utils.c, inout.c
 */
extern node_ptr node_head;

extern node_ptr add_to_nodes_list(tri_pointer,int*,int,VEC*,bin_ptr);
extern int add_conn_tri (node_ptr, tri_pointer, int);
extern int set_adjacent_tris(tri_pointer);
extern int fix_orientation(tri_pointer);
extern int flip_tri(tri_pointer);
extern int set_node_connectivity();
extern int write_node_list();
extern int count_nodes();
extern int write_tri_list (tri_pointer);
extern double find_area(tri_pointer);
extern VEC find_center(tri_pointer);
extern int update_minmax(double*,VEC*,VEC*);
extern VEC cross(VEC,VEC);
extern VEC vscale(double,VEC);
extern VEC from(VEC,VEC);
extern VEC plus(VEC,VEC);
extern VEC midpt(VEC,VEC);
extern VEC norm(VEC);
extern double length(VEC);
extern double lengthsq(VEC);
extern double find_tri_dist(tri_pointer,VEC);
extern VEC find_tri_normal(tri_pointer);
extern VEC find_normal(VEC,VEC,VEC);
extern double dot(VEC,VEC);
extern double theta(VEC,VEC);
extern void prepare_node_bin(bin_ptr,VEC,VEC);
extern float** allocate_2d_array_f(int,int);
extern int free_2d_array_f(float**);
extern tri_pointer read_input(char*,int,tri_pointer);
extern int write_output(tri_pointer,char*,int,int,char**);
extern int find_mesh_stats(char *,VEC*,VEC*,int*,int*);
extern int get_tri(FILE*,int,tri_pointer,char*);
extern int write_tri(FILE*,int,tri_pointer,char*);
extern int inside_bounds(double,double,double);

/* end */
