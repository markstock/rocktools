/*************************************************************
 *
 *  xrayutil.c - Utility subroutines for use primarily 
 *	with rockdetail
 *
 *  Mark J. Stock, mstock@umich.edu
 *
 * rocktools - Tools for creating and manipulating triangular meshes
 * Copyright (C) 2004-15  Mark J. Stock
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
#include "png.h"

#ifdef _OPENMP
  #include <omp.h>
#endif

/* define a C preprocessor variable so that when structs.h is included,
 * it will contain extra information used only by this program */
#define MODULE_ROCKXRAY
#include "structs.h"

/* generate depth-of-field using a slower Gaussian splatting scheme */
//#define USE_GAUSSIAN

png_byte** allocate_2d_array_pb(int,int,int);
int write_png_image(png_byte**,int,int,int,char*,int,int);

/*
 * Write a PGM image of the xray of the shell of a mesh
 *
 * "vz" is the view vector
 * "size" is the final image pixel resolution desired
 * "thick" is the thickness of the mesh, in world units
 * "square" forces a square image, and centers the object (TRUE|FALSE)
 * "thisq" sets quality (0=low, 1=med, 2=high, 3=very high)
 */
int write_xray (tri_pointer tri_head, VEC vz, double *xb, double *yb, int size,
      double thick, int square, double border, int thisq, double peak_crop, double gamma,
      int write_hibit, int is_solid, int num_images, char* prefix, char* output_format,
      int force_num_threads) {

   int write_pgm;			// write a PGM file
   int write_png;			// write a PNG file
   int cnt;
   int xres,yres;			// the actual image size
   float ***a;				// the array to print
   png_byte **img = NULL;		// the png array
   double dtemp,xsize,ysize,zsize,dd,ddz;
   double xmin,xmax,ymin,ymax;		// bounds of the image
   double zmin,zmax;			// bounds in the image direction
   VEC vx,vy;				// image basis vectors
   int num_norm_layers;			// number of subdivisions in normal direction
   tri_pointer this_tri;
   node_ptr this_node;

   int debug_write = FALSE;
   FILE *debug_out;

   if (debug_write)
     debug_out = fopen("temp", "w");

   // set the desired output format
   if (strncmp(output_format, "pgm", 3) == 0) {
      write_pgm = TRUE;
      write_png = FALSE;
   } else if (strncmp(output_format, "png", 3) == 0) {
      write_png = TRUE;
      write_pgm = FALSE;
   } else {
      //fprintf(stderr,"WARNING (write_xray): output file format (%s)\n",output_format);
      //fprintf(stderr,"  unrecognized. Writing PNG by default.\n");
      write_png = TRUE;
      write_pgm = FALSE;
   }

   // now, actually create the image //

   // first, find the three basis vectors: screen-x, screen-y, z (vz)
   vz = norm(vz);
   vy.x = 0.0;
   vy.y = 0.0;
   vy.z = 1.0;
   // one way: -dz was the view "from" vector
   vx = norm(cross(vy,vz));
   vy = norm(cross(vz,vx));
   // the other way: -dz is the view "direction" vector
   // vx = norm(cross(vz,vy));
   // vy = norm(cross(vx,vz));
   //fprintf(stderr,"image axes are %g %g %g  and %g %g %g\n",vx.x,vx.y,vx.z,vy.x,vy.y,vy.z);


   // then, cycle through all elements, projecting the nodes to
   //    the image plane, and determining the aspect ratio needed
   fprintf(stderr,"Determining bounds"); fflush(stderr);
   xmin = 9.9e+9;
   xmax = -9.9e+9;
   ymin = 9.9e+9;
   ymax = -9.9e+9;
   zmin = 9.9e+9;
   zmax = -9.9e+9;
   this_node = node_head;
   cnt = 0;
   while (this_node) {
      // fprintf(stderr,"this node at %g %g %g, projection is %g %g\n",this_node->loc.x,this_node->loc.y,this_node->loc.z,dot(vx,this_node->loc),dot(vy,this_node->loc));
      // project this node to image plane
      dtemp = dot(vx,this_node->loc);
      if (dtemp<xmin) xmin = dtemp;
      if (dtemp>xmax) xmax = dtemp;

      dtemp = dot(vy,this_node->loc);
      if (dtemp<ymin) ymin = dtemp;
      if (dtemp>ymax) ymax = dtemp;

      dtemp = dot(vz,this_node->loc);
      if (dtemp<zmin) zmin = dtemp;
      if (dtemp>zmax) zmax = dtemp;

      if (++cnt%DOTPER == 1) {
         fprintf(stderr,".");
         fflush(stderr);
      }
      this_node = this_node->next_node;
   }
   fprintf(stderr,"\n");
   // and, if x- and y-bounds are used (i.e. if xb[0] is greater than 0),
   //    correct these numbers to either crop off image, or to pad the image
   if (xb[0] > 0.0) {
      xmin = xb[1];
      xmax = xb[2];
   }
   if (yb[0] > 0.0) {
      ymin = yb[1];
      ymax = yb[2];
   }
   fprintf(stderr,"min/max image bounds are %g/%g and %g/%g\n",xmin,xmax,ymin,ymax);


   // determine image bounds, etc
   // no border if image bounds are within geometry
   if (xb[0] > 0.0 || yb[0] > 0.0) border = 0.0;
   xsize = (1.0+border)*(xmax-xmin);
   ysize = (1.0+border)*(ymax-ymin);
   //fprintf(stderr,"geometry size %g x %g\n",xsize,ysize);
   if (xsize > ysize) {
      xres = size;
      if (square) {
         yres = xres;
         ymin = 0.5*(ymax+ymin) - 0.5*xsize;
         xmin = 0.5*(xmax+xmin) - 0.5*xsize;
      } else {
         yres = (int)(xres*(ymax-ymin + border*(xmax-xmin)) / xsize);
         ymin -= 0.5*border*(xmax-xmin);
         xmin -= 0.5*border*(xmax-xmin);
      }
      dd = xsize/xres;
   } else {
      yres = size;
      if (square) {
         xres = yres;
         xmin = 0.5*(xmax+xmin) - 0.5*ysize;
         ymin = 0.5*(ymax+ymin) - 0.5*ysize;
      } else {
         xres = (int)(yres*(xmax-xmin + border*(ymax-ymin)) / ysize);
         xmin -= 0.5*border*(ymax-ymin);
         ymin -= 0.5*border*(ymax-ymin);
      }
      dd = ysize/yres;
   }
   fprintf(stderr,"Final image to be %d x %d\n",xres,yres);
   //fprintf(stderr,"  new xmin,ymin %g %g, dd %g\n",xmin,ymin,dd);


   // and in the z direction?
   zsize = zmax-zmin;
   zmin -= 0.01*zsize;
   zmax += 0.01*zsize;

   // set ddz, a measure of how thick each z-layer is (each output image)
   //   (add 1e-5 so that num_images of 1 doesn't futz things up)
   ddz = (zmax-zmin)/(num_images - 1.0 + 1e-5);


   // find out how many mesh layers we need (thick==-1 if not entered)
   if (thick > 0.0) {
      float fthick;
      //if (hiq) fthick = 3.3*thick/dd;
      //else fthick = 2.0*thick/dd;
      fthick = (2.0+1.3*thisq) * thick/dd;
      num_norm_layers = (int)fthick;
      if (num_norm_layers < 1) num_norm_layers = 1;
      fprintf(stderr,"Using %d layers (%g raw thickness)\n",num_norm_layers,fthick); fflush(stderr);
   } else {
      num_norm_layers = 1;
      thick = 0.0;
   }
   if (is_solid) {
      num_norm_layers = 1;
      thick = 0.0;
   }


   // allocate the array(s)
   a = (float***)malloc(num_images*sizeof(float**));
   for (int i=0; i<num_images; i++) {
      a[i] = allocate_2d_array_f(xres,yres);
   }
   

   // zero out the arrays
   for (int inum=0; inum<num_images; inum++)
      for (int i=0; i<xres; i++)
        for (int j=0; j<yres; j++)
           a[inum][i][j] = 0.0;


   // then, loop through all elements, writing to the image
   fprintf(stderr,"Writing data to image plane"); fflush(stderr);

#ifdef _OPENMP
   // set parallelism
   // how many threads do we need? 2x as many cores, but limit by image size
   int num_threads = (1+xres/1024)*(1+yres/1024);
   if (num_threads > 2*omp_get_num_procs()) num_threads = 2*omp_get_num_procs();
   if (16*num_threads > xres) num_threads = xres/16;
   if (force_num_threads > 0) num_threads = force_num_threads;
   if (num_threads < 1) num_threads = 1;
   if (num_threads > omp_get_max_threads()) num_threads = omp_get_max_threads();
   omp_set_num_threads(num_threads);
   fprintf(stderr," using %d threads",num_threads); fflush(stderr);

   // how many memory blocks do we need?
   int num_locks = 8*num_threads;
   if (num_threads == 1) num_locks = 1;

   // initialize image memory locks
   omp_lock_t memlock[num_locks];
   int lock_min[num_locks];
   int lock_max[num_locks];
   for (int i=0; i<num_locks; i++) {
      omp_init_lock(&memlock[i]);
      lock_min[i] = (i*xres)/num_locks;
      lock_max[i] = ((i+1)*xres)/num_locks - 1;
      //fprintf(stderr,"\n lock %d from %d to %d",i,lock_min[i],lock_max[i]); fflush(stderr);
   }

   // finally, set equal-spaced tri_head for each thread
   // we will not need to do this if OpenMP could do task-parallel (it can)
   tri_pointer tri_heads[num_threads+1];
   // first loop: just count
   cnt = 0;
   this_tri = tri_head;
   while (this_tri) {
      this_tri = this_tri->next_tri;
      cnt++;
   }
   int tris_per_thread = cnt/num_threads;
   // second loop: set pointers
   cnt = 0;
   this_tri = tri_head;
   while (this_tri) {
      if (++cnt%tris_per_thread == 1) {
         //fprintf(stderr,"\n mark tri %d number %d",cnt,cnt/tris_per_thread); fflush(stderr);
         tri_heads[cnt/tris_per_thread] = this_tri;
      }
      this_tri = this_tri->next_tri;
   }
   tri_heads[num_threads] = NULL;
#endif

   // begin parallel section
   // how do we have each thread march through the same linked list?
   // i.e. each thread grabs the next triangle and processes it
#pragma omp parallel private(cnt,this_tri)
{
   cnt = 0;
#ifdef _OPENMP
   this_tri = tri_heads[omp_get_thread_num()];
   while (this_tri != tri_heads[omp_get_thread_num()+1]) {
#else
   this_tri = tri_head;
   while (this_tri) {
#endif

      const node_ptr n0 = this_tri->node[0];
      const node_ptr n1 = this_tri->node[1];
      const node_ptr n2 = this_tri->node[2];

      // first, see if the triangle is anywhere near the actual view window

      // check x-direction first
      double minpos = 9.9e+9;
      double maxpos = -9.9e+9;
      for (int i=0; i<3; i++) {
         // find location in image coordinates
         const double pos = dot(vx,this_tri->node[i]->loc) - xmin;
         if (pos > maxpos) maxpos = pos;
         if (pos < minpos) minpos = pos;
      }
      if ((int)(floor((maxpos+thick)/dd)) < -1 ||
          (int)(floor((minpos-thick)/dd)) > xres+1) {
         // skip this tri
         this_tri = this_tri->next_tri;
         continue;
      }

      // then check y-direction
      minpos = 9.9e+9;
      maxpos = -9.9e+9;
      for (int i=0; i<3; i++) {
         // find location in image coordinates
         const double pos = dot(vy,this_tri->node[i]->loc) - ymin;
         if (pos > maxpos) maxpos = pos;
         if (pos < minpos) minpos = pos;
      }
      if ((int)(floor((maxpos+thick)/dd)) < -1 ||
          (int)(floor((minpos-thick)/dd)) > yres+1) {
         // skip this tri
         this_tri = this_tri->next_tri;
         continue;
      }

      // use x-array coordinates to determine which lock(s) to get;
      //   because the x coord is the slower-varying index in memory
      // so, what?
#ifdef _OPENMP
      int lowbound = floor((maxpos+thick)/dd) - 1;
      int highbound = floor((minpos-thick)/dd) + 1;
      // loop through locks, grabbing the right ones
      for (int i=0; i<num_locks; i++ ) {
         if (lowbound > lock_max[i] || highbound < lock_min[i]) {
            // this tri does not use these pixel rows, skip
         } else {
            // this tri will write to pixels in this band
            omp_set_lock(&memlock[i]);
         }
      }
#endif


      // break the tri down into sub-triangles in the triangle plane
      const double area = find_area(this_tri);
      if (isnan(area)) {
         fprintf(stderr,"\nfound tri with nan area, skipping");
         this_tri = this_tri->next_tri;
         continue;
      }
      double sidelen = sqrt(area);
      if (is_solid) sidelen *= 3.;

      //if (area < min_area) {
         //fprintf(stderr,"\nmin area %g",area);
         //min_area = area;
      //}

      // new method:
      //int subdivisions;
      //if (hiq) subdivisions = (int)(4.0*sidelen/dd);
      //else subdivisions = (int)(2.5*sidelen/dd);
      int subdivisions = (int)((2.5+1.5*thisq)*sidelen/dd);
      if (subdivisions < 1) subdivisions = 1;
      // if (subdivisions > 400) subdivisions = 400;
      //fprintf(stderr,"sidelen/dd is %g, area is %g, sidelen is %g\n",sidelen/dd,area,sidelen);

      // fprintf(stderr,"this tri is at %g %g %g, %g %g %g, %g %g %g\n",n0->loc.x,n0->loc.y,n0->loc.z,n1->loc.x,n1->loc.y,n1->loc.z,n2->loc.x,n2->loc.y,n2->loc.z);

      // now, subdivide in the tri-normal direction to simulate the
      //    thickness of the triangular prism
      VEC trinorm = find_normal(n0->loc,n1->loc,n2->loc);
      // scale the normal vector to half the thickness
      if (!is_solid) { 
         trinorm.x *= 0.5*thick;
         trinorm.y *= 0.5*thick;
         trinorm.z *= 0.5*thick;
      }

      // reset D to be 3*subdivisions
      const int D = 3*subdivisions;

      // base density is a scaled area measure
      double factor = (1.e+5)*area/(double)(subdivisions*subdivisions)/(dd*dd);
      if (is_solid) { 
         // flip if triangle points away from viewer
         // OMG, this is wrong!!!!! It should be sgn(dot())
         // Wait, we're using the absolute area of the tri, so it's OK.
         factor *= dot(trinorm,vz);
      }

      if (debug_write) {
         //fprintf(stderr,"%g %g %g %g\n",n0->loc.x,n0->loc.y,n0->loc.z,dot(trinorm,vz));
         fprintf(debug_out,"%g %g %g %g\n",n0->loc.x,n0->loc.y,n0->loc.z,dot(trinorm,vz));
      }


      // put the value on the grid by subdivision
      for (int k=0; k<num_norm_layers; k++) {

      // displacement of given layer in normal dir.
      const double norm_disp = (double)(2*k+1)/(double)(num_norm_layers) - 1.0;

      for (int i=0; i<subdivisions; i++) {

         for (int j=0; j<(2*i+1); j++) {

            int A,B,C;
            if (j%2 == 0) {
               A = 3*(subdivisions-i)-2;
               B = (i-j/2)*3+1;
               C = D-A-B;
            } else {
               A = 3*(subdivisions-i)-1;
               B = (i-(j+1)/2)*3+2;
               C = D-A-B;
            }
            VEC ec;
            ec.x = (A*n0->loc.x+B*n1->loc.x+C*n2->loc.x)/(double)(D);
            ec.y = (A*n0->loc.y+B*n1->loc.y+C*n2->loc.y)/(double)(D);
            ec.z = (A*n0->loc.z+B*n1->loc.z+C*n2->loc.z)/(double)(D);
            // fprintf(stderr,"   point at %g %g %g\n",ec.x,ec.y,ec.z);

            // and perturb it in the triangle-normal direction
            ec.x += trinorm.x*norm_disp;
            ec.y += trinorm.y*norm_disp;
            ec.z += trinorm.z*norm_disp;

            // find location in image coordinates (vx, vy, vz are normalized)
            double xpos = dot(vx,ec) - xmin;
            double ypos = dot(vy,ec) - ymin;
            double zpos = dot(vz,ec) - zmin;

#ifdef USE_GAUSSIAN
            fprintf(stderr,"GAUSSIAN kernel unsupports with multiple layers\n");
            exit(1);

            // circle of confusion radius in image units (sigma)
            const double rad = fabs(zpos-0.45) + dd;
            const double cnst = 1./pow(rad,2);
            //fprintf(stderr,"%g %g %g  %g\n",xpos,ypos,zpos,rad);
            if (rad < 0.) {
               // change to 5., run at 10x res to make smooth dots
            } else if (rad < 0.03) {
               // one sample per pixel
               int sx = (int)((xpos-3.*rad)/dd);
               if (sx<0) sx=0;
               int ex = (int)((xpos+3.*rad)/dd);
               if (ex>=xres) ex=xres-1;
               int sy = (int)((ypos-3.*rad)/dd);
               if (sy<0) sy=0;
               int ey = (int)((ypos+3.*rad)/dd);
               if (ey>=yres) ey=yres-1;
               //fprintf(stderr,"%d:%d %d:%d  %g %g %g  %g\n",sx,ex,sy,ey,xpos,ypos,zpos,rad);
               for (int ix=sx; ix<ex; ix++) {
                  const double dx = xpos - ix*dd;	// distance in image units
                  const double drr = pow(dx,2);
                  for (int iy=sy; iy<ey; iy++) {
                     const double dy = ypos - iy*dd;
                     const double dr = drr + pow(dy,2);
                     //fprintf(stderr,"  %d %d  %g %g  %g\n",ix,iy,dx,dy,dr);
                     a[0][ix][iy] += factor*cnst*exp(-0.5*dr*cnst);
                  }
               }
                     //exit(0);
            } else {
               // nothing (too diffuse)
            }
#else  // not Gaussian

            // find lower-left pixel coordinate (be able to accept negative quantities)
            const int xloc = (int)(floor(xpos/dd));
            const int yloc = (int)(floor(ypos/dd));
            // fprintf(stderr,"   which is %g %g, cell %d %d\n",xpos,ypos,xloc,yloc);

            // reset xpos, ypos as local cell coordinates
            xpos = xpos/dd - xloc;
            ypos = ypos/dd - yloc;
            // if (xpos < 0 || ypos < 0)
               // fprintf(stderr,"   subcell coords %g %g, weight %g\n",xpos,ypos,factor);

            // write a little blob at each point, use area weighting
            if (num_images == 1) {

               // finally, set the density of this point
               double rfactor = factor;
               if (is_solid) {
                  // zpos is always positive---it's the raw distance from the zmin plane
                  rfactor *= zpos;
               }

               if (xloc > -1 && xloc < xres) {
                  double rtemp = rfactor*(1.0-xpos);
                  if (yloc > -1 && yloc < yres)
                     a[0][xloc][yloc]     += rtemp*(1.0-ypos);
                  if (yloc > -2 && yloc+1 < yres)
                     a[0][xloc][yloc+1]   += rtemp*(ypos);
               }
               if (xloc > -2 && xloc+1 < xres) {
                  double rtemp = rfactor*(xpos);
                  if (yloc > -1 && yloc < yres)
                     a[0][xloc+1][yloc]   += rtemp*(1.0-ypos);
                  if (yloc > -2 && yloc+1 < yres)
                     a[0][xloc+1][yloc+1] += rtemp*(ypos);
               }

            } else {
               // multiple layers: we need to be careful with solid images

               // finally, set the density of this point
               double rfactor = factor;

               // (re)set the z position (indicates which layers/images to which to write)
               const int zloc = (int)(floor(zpos/ddz));
               zpos = zpos/ddz - zloc;
               double rtemp = 0.0;
               double stemp = 0.0;

               if (is_solid) {

                  // NOT DONE!!!
                  // must loop over lots of layers
                  // need to calculate multipliers based on TSC interpolation!!!
                  // choice is between sharp layers (no interp) and soft layers (TSC)
                  double zsq = zpos*zpos;
                  double zinv = 2.0-pow(1.0-zpos,2);

                  if (xloc > -1 && xloc < xres) {
                     rtemp = rfactor*(1.0-xpos);
                     if (yloc > -1 && yloc < yres)
                        stemp = rtemp*(1.0-ypos);
                        a[zloc+1][xloc][yloc]     += stemp*zsq;
                        a[zloc][xloc][yloc]       += stemp*zinv;
                        for (int inum=zloc-1; inum>-1; inum--)
                           a[inum][xloc][yloc]    += stemp*2.0;
                     if (yloc > -2 && yloc+1 < yres)
                        stemp = rtemp*(ypos);
                        a[zloc+1][xloc][yloc+1]   += stemp*zsq;
                        a[zloc][xloc][yloc+1]     += stemp*zinv;
                        for (int inum=zloc-1; inum>-1; inum--)
                           a[inum][xloc][yloc+1]  += stemp*2.0;
                  }
                  if (xloc > -2 && xloc+1 < xres) {
                     rtemp = rfactor*(xpos);
                     if (yloc > -1 && yloc < yres)
                        stemp = rtemp*(1.0-ypos);
                        a[zloc+1][xloc+1][yloc]   += stemp*zsq;
                        a[zloc][xloc+1][yloc]     += stemp*zinv;
                        for (int inum=zloc-1; inum>-1; inum--)
                           a[inum][xloc+1][yloc]  += stemp*2.0;
                     if (yloc > -2 && yloc+1 < yres)
                        stemp = rtemp*(ypos);
                        a[zloc+1][xloc+1][yloc+1] += stemp*zsq;
                        a[zloc][xloc+1][yloc+1]   += stemp*zinv;
                        for (int inum=zloc-1; inum>-1; inum--)
                           a[inum][xloc+1][yloc+1]+= stemp*2.0;
                  }

               } else {

                  if (xloc > -1 && xloc < xres) {
                     rtemp = rfactor*(1.0-xpos);
                     if (yloc > -1 && yloc < yres)
                        stemp = rtemp*(1.0-ypos);
                        a[zloc][xloc][yloc]       += stemp*(1.0-zpos);
                        a[zloc+1][xloc][yloc]     += stemp*(zpos);
                     if (yloc > -2 && yloc+1 < yres)
                        stemp = rtemp*(ypos);
                        a[zloc][xloc][yloc+1]     += stemp*(1.0-zpos);
                        a[zloc+1][xloc][yloc+1]   += stemp*(zpos);
                  }
                  if (xloc > -2 && xloc+1 < xres) {
                     rtemp = rfactor*(xpos);
                     if (yloc > -1 && yloc < yres)
                        stemp = rtemp*(1.0-ypos);
                        a[zloc][xloc+1][yloc]     += stemp*(1.0-zpos);
                        a[zloc+1][xloc+1][yloc]   += stemp*(zpos);
                     if (yloc > -2 && yloc+1 < yres)
                        stemp = rtemp*(ypos);
                        a[zloc][xloc+1][yloc+1]   += stemp*(1.0-zpos);
                        a[zloc+1][xloc+1][yloc+1] += stemp*(zpos);
                  }

               } // if not solid
            } // if multiple images

#endif  // not Gaussian
         }
      }
      }

#ifdef _OPENMP
      // loop through locks, releasing them all
      for (int i=0; i<num_locks; i++ ) {
         omp_unset_lock(&memlock[i]);
      }
#endif

      if (++cnt%DOTPER == 1) {
         fprintf(stderr,".");
         fflush(stderr);
      }

      this_tri = this_tri->next_tri;
   }
} // end omp section
   fprintf(stderr,"\n");

   if (debug_write)
     fclose(debug_out);

   // finally, print the image

   float maxval = 0.;
   if (write_pgm) {

      if (num_images == 1) 
         fprintf(stderr,"Writing PGM image");
      else
         fprintf(stderr,"Writing %d PGM images", num_images);
      fflush(stderr);

      // gamma-correct and check for peak value
      for (int inum=0; inum<num_images; inum++) {
         for (int j=yres-1; j>-1; j--) {
            for (int i=0; i<xres; i++) {
               a[inum][i][j] = exp(gamma*log(a[inum][i][j]));
               if (a[inum][i][j] > maxval) maxval = a[inum][i][j];
            }
         }
      }
      fprintf(stderr,", maxval is %g\n",maxval); fflush(stderr);

      // peak cropping is now a command-line option
      fprintf(stderr,", maxval is %g",maxval);
      if (!is_solid) {
         maxval *= peak_crop;
         fprintf(stderr,", peak-cropped maxval is %g",maxval);
      }
      fprintf(stderr,"\n"); fflush(stderr);

      // scale values and write files
      for (int inum=0; inum<num_images; inum++) {

         FILE* ofp;
         if (num_images == 1) {
            ofp = stdout;
         } else {
            char file_name[255];
            sprintf(file_name, "%s_%02d.png", prefix, inum);
            ofp = fopen(file_name, "wb");
         }

         if (write_hibit) {
            // write header
            fprintf(ofp,"P2\n%d %d\n%d\n",xres,yres,65535);
            // write data
            for (int j=yres-1; j>-1; j--) {
               for (int i=0; i<xres; i++) {
                  int printval = (int)(a[inum][i][j]*65536.0/maxval);
                  if (printval > 65535) printval = 65535;
                  fprintf(ofp,"%d\n",printval);
               }
            }
         } else {
            // write header
            fprintf(ofp,"P2\n%d %d\n%d\n",xres,yres,255);
            // write data
            for (int j=yres-1; j>-1; j--) {
               for (int i=0; i<xres; i++) {
                  int printval = (int)(a[inum][i][j]*256.0/maxval);
                  if (printval > 255) printval = 255;
                  fprintf(ofp,"%d\n",printval);
               }
            }
         }

         if (num_images != 1) {
            fclose(ofp);
         }
      }

   } else if (write_png) {

      if (num_images == 1) 
         fprintf(stderr,"Writing PNG image");
      else
         fprintf(stderr,"Writing %d PNG images", num_images);
      fflush(stderr);

      // gamma-correct and check for peak value
      for (int inum=0; inum<num_images; inum++) {
         for (int i=0; i<xres; i++) {
            for (int j=0; j<yres; j++) {
               a[inum][i][j] = exp(gamma*log(a[inum][i][j]));
               if (a[inum][i][j] > maxval) maxval = a[inum][i][j];
            }
         }
      }

      fprintf(stderr,", maxval is %g",maxval);
      if (!is_solid) {
         maxval *= peak_crop;
         fprintf(stderr,", peak-cropped maxval is %g",maxval);
      }
      fprintf(stderr,"\n"); fflush(stderr);

      // scale all values
      if (write_hibit) {
         img = allocate_2d_array_pb(xres,yres,16);
         for (int inum=0; inum<num_images; inum++) {
            for (int j=yres-1; j>-1; j--) {
               for (int i=0; i<xres; i++) {
                  int printval = (int)(a[inum][i][j]*65536.0/maxval);
                  if (printval<0) printval = 0;
                  if (printval>65535) printval = 65535;
                  img[yres-1-j][2*i] = (png_byte)(printval/256);
                  img[yres-1-j][2*i+1] = (png_byte)(printval%256);
               }
            }
            write_png_image(img,yres,xres,16,prefix,inum,num_images);
         }
      } else {
         img = allocate_2d_array_pb(xres,yres,8);
         for (int inum=0; inum<num_images; inum++) {
            for (int j=yres-1; j>-1; j--) {
               for (int i=0; i<xres; i++) {
                  int printval = (int)(a[inum][i][j]*256.0/maxval);
                  if (printval<0) printval = 0;
                  if (printval>255) printval = 255;
                  img[yres-1-j][i] = (png_byte)printval;
               }
            }
            write_png_image(img,yres,xres,8,prefix,inum,num_images);
         }
      }

   }

   /* replace the old list with the new list */
   return(0);
}


/*
 * allocate memory for a two-dimensional array of png_byte
 */
png_byte** allocate_2d_array_pb(int nx, int ny, int depth) {

   int i,bytesperpixel;
   png_byte **array;

   if (depth <= 8) bytesperpixel = 1;
   else bytesperpixel = 2;
   array = (png_byte **)malloc(ny * sizeof(png_byte *));
   array[0] = (png_byte *)malloc(bytesperpixel * nx * ny * sizeof(png_byte));

   for (i=1; i<ny; i++)
      array[i] = array[0] + i * bytesperpixel * nx;

   return(array);
}


/*
 * write a png file
 */
// int write_png_image(char *file_name,png_byte** image,int xres,int yres,int depth) {
int write_png_image(png_byte** image,int xres,int yres,int depth,char *prefix,int img_num,int num_images) {

   png_uint_32 height,width;
   FILE *fp;
   png_structp png_ptr;
   png_infop info_ptr;
   //png_colorp palette;
   //png_voidp user_error_ptr;
   char file_name[255];

   // if we were given a prefix, write to a file instead of stdout
   if (num_images == 1) {
      fp = stdout;
   } else {
      sprintf(file_name, "%s_%02d.png", prefix, img_num);
      fp = fopen(file_name, "wb");
   }

   if (fp == NULL)
      return (-1);

   height=yres;
   width=xres;

   /* Create and initialize the png_struct with the desired error handler
    * functions.  If you want to use the default stderr and longjump method,
    * you can supply NULL for the last three parameters.  We also check that
    * the library version is compatible with the one used at compile time,
    * in case we are using dynamically linked libraries.  REQUIRED.
    */
   png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,
      NULL, NULL, NULL);
      // user_error_ptr, user_error_fn, user_warning_fn);

   if (png_ptr == NULL)
   {
      fclose(fp);
      return (-1);
   }

   /* Allocate/initialize the image information data.  REQUIRED */
   info_ptr = png_create_info_struct(png_ptr);
   if (info_ptr == NULL)
   {
      fclose(fp);
      png_destroy_write_struct(&png_ptr,(png_infopp)NULL);
      return (-1);
   }

   /* Set error handling.  REQUIRED if you aren't supplying your own
    * error handling functions in the png_create_write_struct() call.
    */
   if (setjmp(png_jmpbuf(png_ptr)))
   {
      /* If we get here, we had a problem reading the file */
      fclose(fp);
      png_destroy_write_struct(&png_ptr, &info_ptr);
      return (-1);
   }

   /* set up the output control if you are using standard C streams */
   png_init_io(png_ptr, fp);

   /* Set the image information here.  Width and height are up to 2^31,
    * bit_depth is one of 1, 2, 4, 8, or 16, but valid values also depend on
    * the color_type selected. color_type is one of PNG_COLOR_TYPE_GRAY,
    * PNG_COLOR_TYPE_GRAY_ALPHA, PNG_COLOR_TYPE_PALETTE, PNG_COLOR_TYPE_RGB,
    * or PNG_COLOR_TYPE_RGB_ALPHA.  interlace is either PNG_INTERLACE_NONE or
    * PNG_INTERLACE_ADAM7, and the compression_type and filter_type MUST
    * currently be PNG_COMPRESSION_TYPE_BASE and PNG_FILTER_TYPE_BASE. REQUIRED
    */
   png_set_IHDR(png_ptr, info_ptr, height, width, depth, PNG_COLOR_TYPE_GRAY,
      PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

   /* Optional gamma chunk is strongly suggested if you have any guess
    * as to the correct gamma of the image.
    */
   png_set_gAMA(png_ptr, info_ptr, 2.2);

   /* Write the file header information.  REQUIRED */
   png_write_info(png_ptr, info_ptr);

   /* One of the following output methods is REQUIRED */
   // png_write_image(png_ptr, row_pointers);
   png_write_image(png_ptr, image);

   /* It is REQUIRED to call this to finish writing the rest of the file */
   png_write_end(png_ptr, info_ptr);

   /* clean up after the write, and free any memory allocated */
   png_destroy_write_struct(&png_ptr, &info_ptr);

   /* close the file, if it was a file */
   if (num_images != 1) fclose(fp);

   /* that's it */
   return (0);
}


