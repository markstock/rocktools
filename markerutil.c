/*************************************************************
 *
 *  markerutil.c - Input and output routines for triangle meshes
 *
 *  Mark J. Stock, mstock@umich.edu
 *
 *
 * rocktools - Tools for creating and manipulating triangular meshes
 * Copyright (C) 2006  Mark J. Stock
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

int write_markers(tri_pointer, MARKER, char[4]);
int write_sphere(VEC, VEC, double, MARKER, char[4]);
int write_rectangle(VEC, VEC, VEC, double, MARKER, char[4]);
int write_cone(VEC, VEC, double, MARKER, char[4]);
int write_rad_rectangle(int, int, VEC, VEC, VEC, VEC);

/*
 * Determine the appropriate output file format from
 * the command-line and write the triangles to stdout
 */
int write_markers(tri_pointer head, MARKER marker, char format[4]) {

   int num_wrote = 0;
   double area,len,remaining_probability;
   VEC loc,norm,lx;
   tri_pointer curr_tri = head;

   // set all node indexes to -1
   curr_tri = head;
   while (curr_tri) {

      // determine the chance that there is a marker for this element
      if (marker.density_type == by_area) {
         area = find_area(curr_tri);
         remaining_probability = marker.density * area;
      } else {
         remaining_probability = marker.density;
      }

      // keep placing markers until there's no chance left
      while (remaining_probability > 0.) {

         // should we place another marker on this triangle?
         if ((1.+rand())/RAND_MAX < remaining_probability) {

            // find the area of the triangle
            area = find_area(curr_tri);

            // find the center of the triangle
            loc = find_center(curr_tri);

            // find the normal
            norm = find_tri_normal(curr_tri);

            // find the vector from node 0 to node 1 (useful for rectangle)
            lx = from(curr_tri->node[0]->loc, curr_tri->node[1]->loc);

            // in order to find the mean length
            len = marker.marker_size * sqrt(area);

            // select on marker type
            if (marker.marker_type == sphere) {
               write_sphere(loc,norm,len,marker,format);
            } else if (marker.marker_type == rectangle) {
               write_rectangle(loc,norm,lx,len,marker,format);
            } else if (marker.marker_type == dualcone) {
               write_cone(loc,norm,len,marker,format);
            } else {
               fprintf(stderr,"Improper marker type (shouldn't get here)\n");
               fprintf(stderr,"Quitting.\n");
               exit(0);
            }

            num_wrote++;
         }

         // reduce the chances of there being another one
         remaining_probability -= 1.;
      }

      // done with this element now
      curr_tri = curr_tri->next_tri;
   }

   return num_wrote;
}


/*
 * Write a sphere at the location given
 */
int write_sphere(VEC loc, VEC norm, double size, MARKER marker, char format[4]) {

   static int count = 0;
   double thissize,dist;
   VEC delta;

   // first, modify size
   if (marker.randomize_size) {
      thissize = size * (marker.min_size +
                         marker.size_range * (1.+rand())/RAND_MAX);
   } else {
      thissize = size;
   }

   // then, modify placement
   if (marker.randomize_radius) {
      // find a random position inside of a unit sphere
      dist = 2.;
      while (dist > 1.) {
         delta.x = -1. + 2.*(1.+rand())/RAND_MAX;
         delta.y = -1. + 2.*(1.+rand())/RAND_MAX;
         delta.z = -1. + 2.*(1.+rand())/RAND_MAX;
         dist = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
      }
      dist = sqrt(dist);

      // project it to the plane
      delta = from(delta, vscale(dot(delta, norm), norm));

      // and move the sphere accordingly
      loc.x += size * delta.x;
      loc.y += size * delta.y;
      loc.z += size * delta.z;
   }

   if (strncmp(format, "rad", 3) == 0) {
      fprintf(stdout,"def sphere s%d\n0 0 4 %g %g %g %g\n",count,loc.x,loc.y,loc.z,0.5*thissize);
   } else {
      fprintf(stderr,"No output file format or unsupported file format given\n");
      fprintf(stderr,"Quitting.\n");
      exit(0);
   }

   count++;

   return count;
}


/*
 * Write a rectangle at the location given
 */
int write_rectangle(VEC loc, VEC normal, VEC lx, double size,
                    MARKER marker, char format[4]) {

   static int count = 0;
   //int num = 0;
   //int i;
   double dist;
   VEC thissize,delta,ly,up;
   VEC node[8];

   // normalize the basis vectors
   lx = norm(lx);
   normal = norm(normal);

   // first, modify size
   if (marker.randomize_size) {
      thissize.x = size * (marker.min_size +
                           marker.size_range * (1.+rand())/RAND_MAX);
      thissize.y = size * (marker.min_size +
                           marker.size_range * (1.+rand())/RAND_MAX);
      thissize.z = size * (marker.min_size +
                           marker.size_range * (1.+rand())/RAND_MAX);
   } else {
      thissize.x = size;
      thissize.y = size;
      thissize.z = size;
   }

   // then, modify placement
   if (marker.randomize_radius) {

      // find a random position inside of a rectangle
      dist = 2.;
      while (dist > 1.) {
         delta.x = -1. + 2.*(1.+rand())/RAND_MAX;
         delta.y = -1. + 2.*(1.+rand())/RAND_MAX;
         delta.z = -1. + 2.*(1.+rand())/RAND_MAX;
         dist = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
      }
      dist = sqrt(dist);

      // project it to the plane
      delta = from(delta, vscale(dot(delta, normal), normal));

      // and move the rectangle accordingly
      loc.x += size * delta.x;
      loc.y += size * delta.y;
      loc.z += size * delta.z;
   }

   // need to write code for randomization of height and rotation!

   // before fiddling with the normal, set the other local axis
   ly = cross(normal, lx);

   // rotate the normal vector, if called for
   if (marker.randomize_normal) {

      // in this incarnation, randomize it more if the normal is nearly +-z

      // set the up vector
      up.x = 0.;
      up.y = 0.;
      up.z = 1.;

      // find the magnitude of "z"-ness
      //dotp = dot(up, normal);

      // find a random direction inside the sphere
      dist = 2.;
      while (dist > 1.) {
         delta.x = -1. + 2.*(1.+rand())/RAND_MAX;
         delta.y = -1. + 2.*(1.+rand())/RAND_MAX;
         delta.z = -1. + 2.*(1.+rand())/RAND_MAX;
         dist = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
      }
      //dist = sqrt(dist);

      // modify the normal direction by this amount
      normal = plus(vscale(marker.normal_pert,delta),normal);

      // and put it back to a norm
      normal = norm(normal);

      // then, reset the lx and ly vectors (before rotating it)
      lx = cross(ly, normal);
      ly = cross(normal, lx);
   }

   // lastly, randomize the local x,y axes (normal is local z)
   if (marker.randomize_rotation) {

      // pick a random rotation, 0 to 2pi
      dist = (2.*M_PI*rand())/RAND_MAX;

      // now, we take cos(dist) of lx and sin(dist) of ly to become the new lx
      lx = plus(vscale(cos(dist),lx),vscale(sin(dist),ly));

      // and recross to get ly
      ly = cross(normal, lx);

      // normalize everything again
      lx = norm(lx);
      ly = norm(ly);
   }

   //fprintf(stderr,"\nn %g %g %g\n",normal.x,normal.y,normal.z);
   //fprintf(stderr,"x %g %g %g\n",lx.x,lx.y,lx.z);
   //fprintf(stderr,"y %g %g %g\n",ly.x,ly.y,ly.z);

   // ------------- done with the randomizing, write it! --------------

   // scale the axes
   ly = norm(ly);
   lx = vscale(thissize.x, lx);
   ly = vscale(thissize.y, ly);
   normal = vscale(thissize.z, normal);

   // Now, write the faces of the rectangle, one at a time
   node[0] = plus(loc, plus(vscale(-0.5, lx), plus(vscale(-0.5, ly), vscale(-0.5, normal))));
   node[1] = plus(loc, plus(vscale(0.5, lx), plus(vscale(-0.5, ly), vscale(-0.5, normal))));
   node[2] = plus(loc, plus(vscale(-0.5, lx), plus(vscale(0.5, ly), vscale(-0.5, normal))));
   node[3] = plus(loc, plus(vscale(0.5, lx), plus(vscale(0.5, ly), vscale(-0.5, normal))));
   node[4] = plus(loc, plus(vscale(-0.5, lx), plus(vscale(-0.5, ly), vscale(0.5, normal))));
   node[5] = plus(loc, plus(vscale(0.5, lx), plus(vscale(-0.5, ly), vscale(0.5, normal))));
   node[6] = plus(loc, plus(vscale(-0.5, lx), plus(vscale(0.5, ly), vscale(0.5, normal))));
   node[7] = plus(loc, plus(vscale(0.5, lx), plus(vscale(0.5, ly), vscale(0.5, normal))));

   if (strncmp(format, "rad", 3) == 0) {
      // maybe don't write the bottom?
      // bottom
      //write_rad_rectangle(count, 0, node[0], node[2], node[3], node[1]);
      // top
      //write_rad_rectangle(count, 1, node[4], node[5], node[7], node[6]);
      write_rad_rectangle(count, 2, node[0], node[1], node[5], node[4]);
      write_rad_rectangle(count, 3, node[1], node[3], node[7], node[5]);
      write_rad_rectangle(count, 4, node[3], node[2], node[6], node[7]);
      write_rad_rectangle(count, 5, node[2], node[0], node[4], node[6]);
   //} else if (strncmp(format, "obj", 1) == 0) {
      //num_wrote = write_obj(head);
   } else {
      fprintf(stderr,"No output file format or unsupported file format given\n");
      fprintf(stderr,"Quitting.\n");
      exit(0);
   }

   count++;

   return count;
}


/*
 * write one rectangle
 */
int write_rad_rectangle(int count, int num, VEC x1, VEC x2, VEC x3, VEC x4) {
   fprintf(stdout,"def polygon p%d_%d 0 0 12\n",count,num);
   fprintf(stdout," %g %g %g\n",x1.x,x1.y,x1.z);
   fprintf(stdout," %g %g %g\n",x2.x,x2.y,x2.z);
   fprintf(stdout," %g %g %g\n",x3.x,x3.y,x3.z);
   fprintf(stdout," %g %g %g\n",x4.x,x4.y,x4.z);
   return(0);
}


/*
 * Write a cone at the location given
 */
int write_cone(VEC loc, VEC normal, double size, MARKER marker, char format[4]) {

   static int count = 0;
   double thissize,dist;
   VEC delta,p1,p2,p3,p4;

   // first, modify size
   if (marker.randomize_size) {
      thissize = size * (marker.min_size +
                         marker.size_range * (1.+rand())/RAND_MAX);
   } else {
      thissize = size;
   }

   // then, modify placement
   if (marker.randomize_radius) {
      // find a random position inside of a unit sphere
      dist = 2.;
      while (dist > 1.) {
         delta.x = -1. + 2.*(1.+rand())/RAND_MAX;
         delta.y = -1. + 2.*(1.+rand())/RAND_MAX;
         delta.z = -1. + 2.*(1.+rand())/RAND_MAX;
         dist = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
      }
      dist = sqrt(dist);

      // project it to the plane
      delta = from(delta, vscale(dot(delta, normal), normal));

      // and move the sphere accordingly
      loc.x += size * delta.x;
      loc.y += size * delta.y;
      loc.z += size * delta.z;
   }

   // find the end point of the cone
   normal = norm(normal);
   p1 = plus(loc, vscale(-size, normal));
   p2 = plus(loc, vscale(-0.5*size, normal));
   p3 = plus(loc, vscale(0.5*size, normal));
   p4 = plus(loc, vscale(size, normal));

   if (strncmp(format, "rad", 3) == 0) {
      fprintf(stdout,"outcolor cone o%d\n0 0 8 %g %g %g %g %g %g %g %g\n",
              count,loc.x,loc.y,loc.z,p4.x,p4.y,p4.z,
              0.6*thissize,0.5*thissize);
      fprintf(stdout,"incolor cup u%d\n0 0 8 %g %g %g %g %g %g %g %g\n",
              count,p3.x,p3.y,p3.z,p4.x,p4.y,p4.z,
              0.4*thissize,0.5*thissize);
      fprintf(stdout,"incolor ring i%d\n0 0 8 %g %g %g %g %g %g %g %g\n",
              count,p3.x,p3.y,p3.z,normal.x,normal.y,normal.z,
              0.,0.4*thissize);
      fprintf(stdout,"outcolor cone o%d\n0 0 8 %g %g %g %g %g %g %g %g\n",
              count,loc.x,loc.y,loc.z,p1.x,p1.y,p1.z,
              0.6*thissize,0.5*thissize);
      fprintf(stdout,"incolor cup u%d\n0 0 8 %g %g %g %g %g %g %g %g\n",
              count,p2.x,p2.y,p2.z,p1.x,p1.y,p1.z,
              0.4*thissize,0.5*thissize);
      fprintf(stdout,"incolor ring i%d\n0 0 8 %g %g %g %g %g %g %g %g\n",
              count,p2.x,p2.y,p2.z,-normal.x,-normal.y,-normal.z,
              0.,0.4*thissize);
   } else {
      fprintf(stderr,"No output file format or unsupported file format given\n");
      fprintf(stderr,"Quitting.\n");
      exit(0);
   }

   count++;

   return count;
}


