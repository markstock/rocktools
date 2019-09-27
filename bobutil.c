/*************************************************************
 *
 *  bobutil.c - Utility subroutines for use with rockbob
 *
 *  Mark J. Stock, mstock@umich.edu
 *
 * rocktools - Tools for creating and manipulating triangular meshes
 * Copyright (C) 2004-14  Mark J. Stock
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
#include <limits.h>


/* define a C preprocessor variable so that when structs.h is included,
 * it will contain extra information used only by this program */
#define MODULE_ROCKBOB
#include "structs.h"

unsigned char*** allocate_3d_array_b(int nx, int ny, int nz) {

   // allocate an array of nx pointers, one for each plane
   unsigned char ***array = (unsigned char ***)malloc(nx * sizeof(char **));

   // are we dealing with more than 2B bytes?
   // On x86-64 Linux, INT_MAX is 2147483647, LONG_MAX is 9223372036854775807
   long int totBytes = (long int)nx * (long int)ny * (long int)nz * (long int)sizeof(char);

   if (totBytes > (long int)INT_MAX/2) {
      // all data is larger than 2GB, allocate it in planes
      fprintf(stderr,"  voxel array is plane-by-plane\n");

      for (int i=0; i<nx; i++) {
         array[i] = (unsigned char **)malloc(ny * sizeof(char *));
         array[i][0] = (unsigned char *)malloc(ny * nz * sizeof(char));
         for (int j=1; j<ny; j++)
            array[i][j] = array[i][0] + j * nz;
      }

   } else {
      // we can fit it all in one malloc (it's under 2GB)
      fprintf(stderr,"  voxel array is monolithic\n");

      array[0] = (unsigned char **)malloc(nx * ny * sizeof(char *));
      for (int i=1; i<nx; i++)
         array[i] = array[0] + i * ny;

      array[0][0] = (unsigned char *)malloc(nx * ny * nz * sizeof(char));
      for (int i=0; i<nx; i++) {
         if (i!=0)
            array[i][0] = array[0][0] + i * ny * nz;
         for (int j=1; j<ny; j++)
            array[i][j] = array[i][0] + j * nz;
      }
   }

   return(array);
}

int free_3d_array_b (unsigned char*** array, int nx, int ny, int nz){
   long int totBytes = (long int)nx * (long int)ny * (long int)nz * (long int)sizeof(char);
   if (totBytes > (long int)INT_MAX/2) {
      for (int i=0; i<nx; i++) {
         free(array[i][0]);
         free(array[i]);
      }
   } else {
      free(array[0][0]);
      free(array[0]);
   }
   free(array);
   return(0);
}

unsigned char** allocate_2d_array_b(int nx, int ny) {

   unsigned char **array = (unsigned char **)malloc(nx * sizeof(unsigned char *));
   array[0] = (unsigned char *)malloc(nx * ny * sizeof(unsigned char));
   for (int i=1; i<nx; i++)
      array[i] = array[0] + i * ny;

   return(array);
}

int free_2d_array_b(unsigned char** array){
   free(array[0]);
   free(array);
   return(0);
}

// define the possible output file types
typedef enum output_format_type {
   noout,       // default is no output
   bob,         // brick of bytes (0..255)
   bos,         // brock of shorts (0..65535)
   bof          // brock of floats (IEEE floats)
} OUT_FORMAT;


/* Function to find minimum of x and y */
int min(int x, int y)
{
  return y ^ ((x ^ y) & -(x < y));
}
 
/* Function to find maximum of x and y */
int max(int x, int y)
{
  return x ^ ((x ^ y) & -(x < y));
}

/*
 * Write a 3D brick of bytes
 */
int write_bob_file_from_uchar(FILE* ofp, unsigned char*** z, int nx, int ny, int nz) {

   /* write header */
   fwrite(&nx,sizeof(int),1,ofp);
   fwrite(&ny,sizeof(int),1,ofp);
   fwrite(&nz,sizeof(int),1,ofp);

   /* write the data */
   for (int i=0; i<nx; i++) {
      fwrite(&z[i][0][0],sizeof(unsigned char),ny*nz,ofp);
   }

   /* return 0 if all went well */
   return(0);
}

//
// minimum distance from point to arbitrary triangle
// code from Omegaflow v2, Surface.F90, function pointElemDistance3d
// which is from fmmbem/src/elemnode.f:1563
//
double mdtri(double x1x, double x1y, double x1z,
             double x2x, double x2y, double x2z,
             double x3x, double x3y, double x3z,
             double px, double py, double pz) {

   // initialize with impossibly large value
   double minDist = 9.9e+9;

   // check against three corners and edges, distance first, then parameter
   double ax,ay,az;
   double bx,by,bz;
   double rx,ry,rz,rn;

   // 1(2)
   ax = px-x1x;
   ay = py-x1y;
   az = pz-x1z;
   const double d1 = pow(ax,2) + pow(ay,2) + pow(az,2);
   if (d1 < minDist) minDist = d1;
   bx = x2x-x1x;
   by = x2y-x1y;
   bz = x2z-x1z;
   const double l1 = pow(bx,2) + pow(by,2) + pow(bz,2);
   rx = by*az - ay*bz;
   ry = bz*ax - az*bx;
   rz = bx*ay - ax*by;
   rn = (pow(rx,2) + pow(ry,2) + pow(rz,2)) / l1;
   if (rn < minDist) {
      const double t = ( ax*bx + ay*by + az*bz ) / l1;
      if (t > 0.0 && t < 1.0) {
         minDist = rn;
      }
   }

   // 2(3)
   ax = px-x2x;
   ay = py-x2y;
   az = pz-x2z;
   const double d2 = pow(ax,2) + pow(ay,2) + pow(az,2);
   if (d2 < minDist) minDist = d2;
   bx = x3x-x2x;
   by = x3y-x2y;
   bz = x3z-x2z;
   const double l2 = pow(bx,2) + pow(by,2) + pow(bz,2);
   rx = by*az - ay*bz;
   ry = bz*ax - az*bx;
   rz = bx*ay - ax*by;
   rn = (pow(rx,2) + pow(ry,2) + pow(rz,2)) / l2;
   if (rn < minDist) {
      const double t = ( ax*bx + ay*by + az*bz ) / l2;
      if (t > 0.0 && t < 1.0) {
         minDist = rn;
      }
   }

   // 3(1)
   ax = px-x3x;
   ay = py-x3y;
   az = pz-x3z;
   const double d3 = pow(ax,2) + pow(ay,2) + pow(az,2);
   if (d3 < minDist) minDist = d3;
   bx = x1x-x3x;
   by = x1y-x3y;
   bz = x1z-x3z;
   const double l3 = pow(bx,2) + pow(by,2) + pow(bz,2);
   rx = by*az - ay*bz;
   ry = bz*ax - az*bx;
   rz = bx*ay - ax*by;
   rn = (pow(rx,2) + pow(ry,2) + pow(rz,2)) / l3;
   if (rn < minDist) {
      const double t = ( ax*bx + ay*by + az*bz ) / l3;
      if (t > 0.0 && t < 1.0) {
         minDist = rn;
      }
   }

   // finally, check against the prism extending from the tri face

   // first, find the triangle normal
   ax = x2x-x3x;
   ay = x2y-x3y;
   az = x2z-x3z;
   bx = x1x-x3x;
   by = x1y-x3y;
   bz = x1z-x3z;
   rx = by*az - ay*bz;
   ry = bz*ax - az*bx;
   rz = bx*ay - ax*by;
   rn = 1.0/sqrt(pow(rx,2) + pow(ry,2) + pow(rz,2));
   const double nx = rx*rn;
   const double ny = ry*rn;
   const double nz = rz*rn;

   // then find raw distance of the point to that plane
   const double cx = (x1x + x2x + x3x) / 3.0;
   const double cy = (x1y + x2y + x3y) / 3.0;
   const double cz = (x1z + x2z + x3z) / 3.0;
   rx = px-cx;
   ry = py-cy;
   rz = pz-cz;
   rn = pow( nx*rx + ny*ry + nz*rz, 2 );

   if (rn < minDist) {

      // now, is the point in the prism of the triangle?
      ax = x2x-x1x;
      ay = x2y-x1y;
      az = x2z-x1z;
      const double an = 1.0/sqrt(pow(ax,2) + pow(ay,2) + pow(az,2));
      ax *= an;
      ay *= an;
      az *= an;

      bx = ny*az - ay*nz;
      by = nz*ax - az*nx;
      bz = nx*ay - ax*ny;
      const double bn = 1.0/sqrt(pow(bx,2) + pow(by,2) + pow(bz,2));
      bx *= bn;
      by *= bn;
      bz *= bn;

      const double xx1 = (x1x-cx)*ax + (x1y-cy)*ay + (x1z-cz)*az;
      const double yy1 = (x1x-cx)*bx + (x1y-cy)*by + (x1z-cz)*bz;
      const double xx2 = (x2x-cx)*ax + (x2y-cy)*ay + (x2z-cz)*az;
      const double yy2 = (x2x-cx)*bx + (x2y-cy)*by + (x2z-cz)*bz;
      const double xx3 = (x3x-cx)*ax + (x3y-cy)*ay + (x3z-cz)*az;
      const double yy3 = (x3x-cx)*bx + (x3y-cy)*by + (x3z-cz)*bz;

      const double den = (xx1-xx3) * (yy2-yy3) - (xx2-xx3) * (yy1-yy3);

      // projection of coordinates in the plane
      const double x = rx*ax + ry*ay + rz*az;
      const double y = rx*bx + ry*by + rz*bz;

      // area coordinates of the projected point
      const double p1 = ((yy2-yy3)*(x-xx3)-(xx2-xx3)*(y-yy3)) / den;
      const double p2 = ((xx1-xx3)*(y-yy3)-(yy1-yy3)*(x-xx3)) / den;
      const double p3 = 1.0-p1-p2;

      if (p1 >= 0.0 && p2 >= 0.0 && p3 >= 0.0) {
         minDist = rn;
      }
   }

   return sqrt(minDist);
}

//
// min dist code from
// http://stackoverflow.com/questions/849211/shortest-distance-between-a-point-and-a-line-segment
//
double minimum_distance(double vx, double vy, double vz,
                        double wx, double wy, double wz,
                        double px, double py, double pz) {

  //fprintf(stderr,"v %g %g %g\n",vx,vy,vz);
  //fprintf(stderr,"w %g %g %g\n",wx,wy,wz);
  //fprintf(stderr,"p %g %g %g\n",px,py,pz);

  // Return minimum distance between line segment vw and point p
  // i.e. |w-v|^2 -  avoid a sqrt
  const double l2 = pow(vx-wx,2) + pow(vy-wy,2) + pow(vz-wz,2);

  // v == w case (sphere)
  if (l2 == 0.0) return sqrt( pow(vx-px,2) + pow(vy-py,2) + pow(vz-pz,2) );

  // Consider the line extending the segment, parameterized as v + t (w - v).
  // We find projection of point p onto the line. 
  // It falls where t = [(p-v) . (w-v)] / |w-v|^2
  const double t = ( (px-vx)*(wx-vx) + (py-vy)*(wy-vy) + (pz-vz)*(wz-vz) ) / l2;

  // Beyond the 'v' end of the segment
  if (t < 0.0) return sqrt( pow(vx-px,2) + pow(vy-py,2) + pow(vz-pz,2) );

  // Beyond the 'w' end of the segment
  else if (t > 1.0) return sqrt( pow(wx-px,2) + pow(wy-py,2) + pow(wz-pz,2) );

  // Projection falls on the segment
  const double jx = vx + t * (wx - vx);
  const double jy = vy + t * (wy - vy);
  const double jz = vz + t * (wz - vz);
  return sqrt( pow(jx-px,2) + pow(jy-py,2) + pow(jz-pz,2) );
}

/*
 * Write a voxel of the shell of a mesh
 *
 * "dx" is the voxel size
 * "thick" is the thickness of the mesh, in world units
 */
int write_bob (tri_pointer tri_head, double *xb, double *yb, double *zb,
      double dx, double thick, int diffuseSteps, char* output_format) {

   int nx, ny, nz;
   double start[3];
   double size[3];
   unsigned char*** dat = NULL;

   double xmin,xmax,ymin,ymax;		// bounds of the image
   double zmin,zmax;			// bounds in the image direction
   tri_pointer this_tri;
   OUT_FORMAT outType = noout;

   int debug_write = FALSE;
   FILE *debug_out;

   if (debug_write)
     debug_out = fopen("temp", "w");

   // set the desired output format
   if (strncmp(output_format, "bob", 3) == 0) {
      outType = bob;
   } else if (strncmp(output_format, "bos", 3) == 0) {
      outType = bos;
   } else if (strncmp(output_format, "bof", 3) == 0) {
      outType = bof;
   } else {
      fprintf(stderr,"WARNING (write_bob): output file format (%s)\n",output_format);
      fprintf(stderr,"  unrecognized. Writing bob by default.\n");
      outType = bob;
   }

   // now, actually create the data //

   // cycle through all elements, determining the aspect ratio needed
   fprintf(stderr,"Determining bounds"); fflush(stderr);
   xmin = 9.9e+9;
   xmax = -9.9e+9;
   ymin = 9.9e+9;
   ymax = -9.9e+9;
   zmin = 9.9e+9;
   zmax = -9.9e+9;
   node_ptr this_node = node_head;
   int cnt = 0;
   while (this_node) {
      double dtemp = this_node->loc.x;
      if (dtemp<xmin) xmin = dtemp;
      if (dtemp>xmax) xmax = dtemp;

      dtemp = this_node->loc.y;
      if (dtemp<ymin) ymin = dtemp;
      if (dtemp>ymax) ymax = dtemp;

      dtemp = this_node->loc.z;
      if (dtemp<zmin) zmin = dtemp;
      if (dtemp>zmax) zmax = dtemp;

      if (++cnt%DOTPER == 1) {
         fprintf(stderr,".");
         fflush(stderr);
      }

      this_node = this_node->next_node;
   }
   fprintf(stderr,"\n");

   start[0] = xmin;
   start[1] = ymin;
   start[2] = zmin;
   size[0] = xmax - xmin;
   size[1] = ymax - ymin;
   size[2] = zmax - zmin;

   fprintf(stderr,"  min/max object bounds are %g/%g %g/%g %g/%g\n",
           start[0],start[0]+size[0],
           start[1],start[1]+size[1],
           start[2],start[2]+size[2]);

   if (dx < 0.0) {
      // dx wasn't set, set it here
      if (size[0] > size[1] && size[0] > size[2]) {
         dx = size[0] / 100.;
      } else if (size[1] > size[0] && size[1] > size[2]) {
         dx = size[1] / 100.;
      } else {
         dx = size[2] / 100.;
      }
   }

   if (thick < 0.0) {
      // thickness wasn't set, set it here
      thick = 2. * dx;
   }

   //const double bufferSize = thick + 2.0 * dx;
   const double bufferSize = thick + dx*(double)(2+diffuseSteps);

   // determine volume bounds and resolution
   if (xb[0] > 0.0) {
      start[0] = xb[1];
      size[0] = xb[2];
   } else {
      start[0] -= bufferSize;
      size[0] += 2.0 * bufferSize;
   }
   nx = size[0] / dx;

   if (yb[0] > 0.0) {
      start[1] = yb[1];
      size[1] = yb[2];
   } else {
      start[1] -= bufferSize;
      size[1] += 2.0 * bufferSize;
   }
   ny = size[1] / dx;

   if (zb[0] > 0.0) {
      start[2] = zb[1];
      size[2] = zb[2];
   } else {
      start[2] -= bufferSize;
      size[2] += 2.0 * bufferSize;
   }
   nz = size[2] / dx;

   fprintf(stderr,"  min/max volume bounds are %g/%g %g/%g %g/%g\n",
           start[0],start[0]+size[0],
           start[1],start[1]+size[1],
           start[2],start[2]+size[2]);

   fprintf(stderr,"  brick will be %d x %d x %d\n",nx,ny,nz);

   // sanity check on bob size
   if (nx*(float)ny*nz > 1.0e+11 || nx > 100000 || ny > 100000 || nz > 100000) {
      fprintf(stderr,"Will not write brick-of-bytes file that large.\n");
      fflush(stderr);
      return(1);
   }

   // allocate the array(s)
   dat = allocate_3d_array_b (nx, ny, nz);

   // zero out the array
   for (int i=0; i<nx; i++) for (int j=0; j<ny; j++) for (int k=0; k<nz; k++) dat[i][j][k] = 0.0;


   // then, loop through all elements, writing to the image
   fprintf(stderr,"Writing data to voxels"); fflush(stderr);

   cnt = 0;
   this_tri = tri_head;
   while (this_tri) {

      const node_ptr n0 = this_tri->node[0];
      const node_ptr n1 = this_tri->node[1];
      const node_ptr n2 = this_tri->node[2];

      const double rad = thick/dx;

      // scale the tri into grid coords
      const double x1 = (n0->loc.x - start[0]) / dx;
      const double y1 = (n0->loc.y - start[1]) / dx;
      const double z1 = (n0->loc.z - start[2]) / dx;
      const double x2 = (n1->loc.x - start[0]) / dx;
      const double y2 = (n1->loc.y - start[1]) / dx;
      const double z2 = (n1->loc.z - start[2]) / dx;
      const double x3 = (n2->loc.x - start[0]) / dx;
      const double y3 = (n2->loc.y - start[1]) / dx;
      const double z3 = (n2->loc.z - start[2]) / dx;

      // find x,y,z range affected by this segment
      const int imin = max((int)floor(fmin(x1-rad, fmin(x2-rad, x3-rad))) - 1, 0);
      const int imax = min((int)ceil(fmax(x1+rad, fmax(x2+rad, x3+rad))) + 1, nx);
      const int jmin = max((int)floor(fmin(y1-rad, fmin(y2-rad, y3-rad))) - 1, 0);
      const int jmax = min((int)ceil(fmax(y1+rad, fmax(y2+rad, y3+rad))) + 1, ny);
      const int kmin = max((int)floor(fmin(z1-rad, fmin(z2-rad, z3-rad))) - 1, 0);
      const int kmax = min((int)ceil(fmax(z1+rad, fmax(z2+rad, z3+rad))) + 1, nz);

      // loop over that subblock
      for (int i=imin; i<imax; i++) {
      for (int j=jmin; j<jmax; j++) {
      for (int k=kmin; k<kmax; k++) {
         // how far is this node from the segment, in voxels?
         double thisDist = 2.0;

         // find distance from triangle to node center, in cell units
         //thisDist = minimum_distance(x1,y1,z1, x2,y2,z2, (double)i+0.5,(double)j+0.5,(double)k+0.5) - rad;
         thisDist = mdtri(x1,y1,z1, x2,y2,z2, x3,y3,z3, (double)i+0.5, (double)j+0.5, (double)k+0.5);
         thisDist -= rad;

         // convert that distance to an unsigned char
         // 255=inside
         // 127=right on the boundary
         // 0=more than 1 voxel away from surface
         unsigned char thisChar = 0;
         if (thisDist < -1.0) thisChar = 255;
         else if (thisDist < 1.0) thisChar = 255 - (unsigned char)(255.0*0.5*(1.0 + thisDist));

         // only update the array if this voxel is nearer to this segment
         if (thisChar > dat[i][j][k]) dat[i][j][k] = thisChar;
      }
      }
      }

      // advance
      this_tri = this_tri->next_tri;
 
      if (++cnt%DOTPER == 1) {
         fprintf(stderr,".");
         fflush(stderr);
      }
   }
   fprintf(stderr,"\n");
   fflush(stderr);

   if (debug_write)
     fclose(debug_out);


   // optionally diffuse the brick-of-whatevers
  // smooth the bof in-place (eventually put this in bobtools)
  if (diffuseSteps > 0) {

    // allocate temporary bof array - but only one plane!
    unsigned char** temp1 = allocate_2d_array_b(ny,nz);
    unsigned char** temp2 = allocate_2d_array_b(ny,nz);

    fprintf(stderr,"diffusing");
    fflush(stderr);
    for (int iter=0; iter<diffuseSteps; iter++) {
      fprintf(stderr,".");
      fflush(stderr);

      // copy first plane into temp1
      for (int j=0; j<ny; j++)
        for (int k=0; k<nz; k++)
          temp1[j][k] = dat[0][j][k];

      // iterate through planes
      for (int ix=1; ix<nx-1; ix++) {

        // do the diffusion, put it in temp2
        for (int iy=1; iy<ny-1; iy++) {
        for (int iz=1; iz<nz-1; iz++) {
          unsigned int neibsum = (unsigned int)dat[ix][iy][iz+1]
                         +(unsigned int)dat[ix][iy][iz-1]
                         +(unsigned int)dat[ix][iy+1][iz]
                         +(unsigned int)dat[ix][iy-1][iz]
                         +(unsigned int)dat[ix+1][iy][iz]
                         +(unsigned int)temp1[iy][iz];
          temp2[iy][iz] = (unsigned char)((neibsum + 6*(unsigned int)dat[ix][iy][iz] + 6) / 12);
        }
        }

        // we can overwrite plane ix-1 now
        for (int j=0; j<ny; j++)
          for (int k=0; k<nz; k++)
            dat[ix-1][j][k] = temp1[j][k];

        // and swap planes
        for (int j=0; j<ny; j++)
          for (int k=0; k<nz; k++)
            temp1[j][k] = temp2[j][k];
      }

      // copy temp1 into last plane
      for (int j=0; j<ny; j++)
        for (int k=0; k<nz; k++)
          dat[nx-2][j][k] = temp1[j][k];
    }

    free_2d_array_b(temp1);
    free_2d_array_b(temp2);
    fprintf(stderr,"\n");
    fflush(stderr);
  }


  // then, expand to allow 3D printing
  if (FALSE) {
    fprintf(stderr,"Growing base to avoid overhangs\n"); fflush(stderr);

    // depending on the angle (this should work for anything steeper than 45 degrees (45-90)
    //const float angle = 45;

    // compute adjacent and diagonal weights
    //const float tana = tan((90.0-angle)*M_PI/180.0);
    //const float aw1 = tana;
    //const float aw2 = 1.0-tana;
    //const float tast = tana/sqrt(2.0);
    //const float dw1 = tast*tast;
    //const float dw2 = tast*(1.0-tast);
    //const float dw3 = (1.0-tast)*(1.0-tast);

    // iterate through z-planes, from top to bottom
    for (int iz=nz-2; iz>0; --iz) {

      // for this layer, look up (+z) for data
      for (int ix=1; ix<nx-1; ++ix) {
      for (int iy=1; iy<ny-1; ++iy) {
        // enforce 45 degree angle (255 = inside object)
        // linearly interpolate to find diagonal values (as they are farther than 1 dx away)
        const int ne = (int)(0.5000*(float)dat[ix+1][iy+1][iz+1] + 0.0858*(float)dat[ix][iy][iz+1] +
                             0.2071*(float)dat[ix+1][iy][iz+1] + 0.2071*(float)dat[ix][iy+1][iz+1]);
        const int nw = (int)(0.5000*(float)dat[ix-1][iy+1][iz+1] + 0.0858*(float)dat[ix][iy][iz+1] +
                             0.2071*(float)dat[ix-1][iy][iz+1] + 0.2071*(float)dat[ix][iy+1][iz+1]);
        const int sw = (int)(0.5000*(float)dat[ix-1][iy-1][iz+1] + 0.0858*(float)dat[ix][iy][iz+1] +
                             0.2071*(float)dat[ix-1][iy][iz+1] + 0.2071*(float)dat[ix][iy-1][iz+1]);
        const int se = (int)(0.5000*(float)dat[ix+1][iy-1][iz+1] + 0.0858*(float)dat[ix][iy][iz+1] +
                             0.2071*(float)dat[ix+1][iy][iz+1] + 0.2071*(float)dat[ix][iy-1][iz+1]);
        // the adjacent columns are easy
        const int xneib = min((int)dat[ix-1][iy][iz+1], (int)dat[ix+1][iy][iz+1]);
        const int yneib = min((int)dat[ix][iy-1][iz+1], (int)dat[ix][iy+1][iz+1]);
        const int aneib = min(sw, ne);
        const int bneib = min(se, nw);
        const int hneib = min(xneib, yneib);
        const int dneib = min(aneib, bneib);
        const int allnb = min(hneib, dneib);
        const int currv = (int)dat[ix][iy][iz];
        dat[ix][iy][iz] = (unsigned char)max(currv, allnb);
      }
      }
    }
  }


  // finally, print the image
  if (outType == bob) {
    fprintf(stderr,"Writing BOB file"); fflush(stderr);

    // finally, write the file
    FILE *ofp = stdout;
    (void) write_bob_file_from_uchar(ofp, dat, nx, ny, nz);

    // free the memory and return
    free_3d_array_b(dat, nx, ny, nz);

    fprintf(stderr,"\n");
    fflush(stderr);
  } else {
    fprintf(stderr,"Output file format unsupported.\n"); fflush(stderr);
  }

  // replace the old list with the new list
  return(0);
}

