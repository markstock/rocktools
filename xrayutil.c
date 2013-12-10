/*************************************************************
 *
 *  xrayutil.c - Utility subroutines for use primarily 
 *	with rockdetail
 *
 *  Mark J. Stock, mstock@umich.edu
 *
 * rocktools - Tools for creating and manipulating triangular meshes
 * Copyright (C) 2004-13  Mark J. Stock
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


/* define a C preprocessor variable so that when structs.h is included,
 * it will contain extra information used only by this program */
#define MODULE_ROCKXRAY
#include "structs.h"

png_byte** allocate_2d_array_pb(int,int,int);
int write_png_image(png_byte**,int,int,int);

/*
 * Write a PGM image of the xray of the shell of a mesh
 *
 * "vz" is the view vector
 * "size" is the final image pixel resolution desired
 * "thick" is the thickness of the mesh, in world units
 * "square" forces a square image, and centers the object (TRUE|FALSE)
 * "hiq" toggles higher-quality rendering (TRUE|FALSE)
 */
int write_xray (tri_pointer tri_head, VEC vz, double *xb, double *yb, int size,
      double thick, int square, double border, int hiq, double peak_crop,
      int write_hibit, int is_solid, char* output_format) {

   // int is_solid = TRUE;		// xray interior, not just boundary
   int write_pgm;			// write a PGM file
   int write_png;			// write a PNG file
   //int write_hibit = TRUE;		// write a 16-bit image (instead of 8)
   int i,j,k,cnt;
   int xres,yres;			// the actual image size
   int subdivisions;
   float **a;				// the array to print
   png_byte **img = NULL;		// the png array
   double dtemp,xsize,ysize,dd;
   double xmin,xmax,ymin,ymax;		// bounds of the image
   double zmin,zmax;			// bounds in the image direction
   double xpos,ypos,factor,rfactor;
   double zpos;
   double proj_area,area,sidelen;
   int A,B,C,D;
   int xloc,yloc;			// the actual image size
   int printval;
   float maxval = 0.;
   VEC vx,vy,ec;			// image basis vectors
   int num_layers;			// number of subdivisions in normal direction
   double norm_disp;			// displacement of given layer in normal dir.
   VEC trinorm;				// triangle normal vector
   tri_pointer this_tri;
   node_ptr this_node,n0,n1,n2;
//#define USE_GAUSSIAN
#ifdef USE_GAUSSIAN
   int ix,iy,sx,sy,ex,ey;
   double dx,dy,dr,drr,rad,cnst;
#endif

   double pos,maxpos,minpos;
   //double min_area = 9.9e+9;
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
   i = 0;
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
      if (i/DOTPER == (i+DPMO)/DOTPER) {
         fprintf(stderr,".");
         fflush(stderr);
      }
      i++;
      this_node = this_node->next_node;
   }
   fprintf(stderr,"\n");
   // and, if x- and y-bounds are used, correct these numbers to either crop off
   //    image, or to pad the image
   if (xb[0] > 0.0) {
      xmin = xb[1];
      xmax = xb[2];
   }
   if (yb[0] > 0.0) {
      ymin = yb[1];
      ymax = yb[2];
   }
   // fprintf(stderr,"min/max image bounds are %g/%g and %g/%g\n",xmin,xmax,ymin,ymax);


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
         yres = (int)(size*((ymax-ymin)/(xmax-xmin)) + border*(xmax-xmin));
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
         xres = (int)(size*((xmax-xmin)/(ymax-ymin)) + border*(ymax-ymin));
         xmin -= 0.5*border*(ymax-ymin);
         ymin -= 0.5*border*(ymax-ymin);
      }
      dd = ysize/yres;
   }
   //fprintf(stderr,"Final image to be %d x %d\n",xres,yres);
   //fprintf(stderr,"  new xmin,ymin %g %g, dd %g\n",xmin,ymin,dd);

   // find out how many layers we need (thick==-1 if not entered)
   if (thick > 0.0) {
      if (hiq) num_layers = (int)(3.3*thick/dd);
      else num_layers = (int)(2*thick/dd);
      if (num_layers < 1) num_layers = 1;
   } else {
      num_layers = 1;
      thick = 0.0;
   }
   if (is_solid) {
      num_layers = 1;
      thick = 0.0;
   }
   fprintf(stderr,"Using %d layers\n",num_layers); fflush(stderr);


   // allocate the array(s)
   a = allocate_2d_array_f(xres,yres);
   if (write_png) {
      if (write_hibit)
         img = allocate_2d_array_pb(xres,yres,16);
      else
         img = allocate_2d_array_pb(xres,yres,8);
   }

   // zero out the array
   for (i=0; i<xres; i++) for (j=0; j<yres; j++) a[i][j] = 0.0;


   // then, loop through all elements, writing to the image
   fprintf(stderr,"Writing data to image plane"); fflush(stderr);

   cnt = 0;
   this_tri = tri_head;
   while (this_tri) {

      n0 = this_tri->node[0];
      n1 = this_tri->node[1];
      n2 = this_tri->node[2];

      // first, see if the triangle is anywhere near the actual view window

      // check x-direction first
      minpos = 9.9e+9;
      maxpos = -9.9e+9;
      for (i=0; i<3; i++) {
         // find location in image coordinates
         pos = dot(vx,this_tri->node[i]->loc) - xmin;
         if (pos > maxpos) maxpos = pos;
         if (pos < minpos) minpos = pos;
      }
      //if (maxpos+thick < -0.1*xsize || minpos-thick > 1.1*xsize) {
      if ((int)(floor((maxpos+thick)/dd)) < -1 ||
          (int)(floor((minpos-thick)/dd)) > xres+1) {
         // skip this tri
         this_tri = this_tri->next_tri;
         continue;
      }

      // then check y-direction
      minpos = 9.9e+9;
      maxpos = -9.9e+9;
      for (i=0; i<3; i++) {
         // find location in image coordinates
         pos = dot(vy,this_tri->node[i]->loc) - ymin;
         if (pos > maxpos) maxpos = pos;
         if (pos < minpos) minpos = pos;
      }
      //if (maxpos+thick < -0.1*ysize || minpos-thick > 1.1*ysize) {
      if ((int)(floor((maxpos+thick)/dd)) < -1 ||
          (int)(floor((minpos-thick)/dd)) > yres+1) {
         // skip this tri
         this_tri = this_tri->next_tri;
         continue;
      }

      // break the tri down into sub-triangles in the triangle plane
      area = find_area(this_tri);
      if (isnan(area)) {
         fprintf(stderr,"\nfound tri with nan area, skipping");
         this_tri = this_tri->next_tri;
         continue;
      }
      sidelen = sqrt(area);

      //if (area < min_area) {
         //fprintf(stderr,"\nmin area %g",area);
         //min_area = area;
      //}

      // new method:
      if (is_solid) sidelen *= 3.;
      if (hiq) subdivisions = (int)(4.0*sidelen/dd);
      else subdivisions = (int)(2.5*sidelen/dd);
      if (subdivisions < 1) subdivisions = 1;
      // if (subdivisions > 400) subdivisions = 400;
      //fprintf(stderr,"sidelen/dd is %g, area is %g, sidelen is %g\n",sidelen/dd,area,sidelen);

      // for is_solid, find the area of the triangle in the image plane
      

      // now, subdivide in the tri-normal direction to simulate the
      //    thickness of the triangular prism
      trinorm = find_normal(n0->loc,n1->loc,n2->loc);
      // scale the normal vector to half the thickness
      if (!is_solid) { 
         trinorm.x *= 0.5*thick;
         trinorm.y *= 0.5*thick;
         trinorm.z *= 0.5*thick;
      }

      // reset D to be 3*subdivisions
      D = 3*subdivisions;
      factor = (1.0e+5)*area/(double)(subdivisions*subdivisions)/(dd*dd);
      // fprintf(stderr,"this tri is at %g %g %g, %g %g %g, %g %g %g\n",n0->loc.x,n0->loc.y,n0->loc.z,n1->loc.x,n1->loc.y,n1->loc.z,n2->loc.x,n2->loc.y,n2->loc.z);

      // for is_solid, determine which way the triangle is facing
      if (is_solid) { 
         proj_area = area*dot(trinorm,vz);
         factor = (1.e+5)*proj_area/(double)(subdivisions*subdivisions)/(dd*dd);
         // if (dot(trinorm,vz) > 0.) factor = factor*-1.0;
      }

      if (debug_write) {
         //fprintf(stderr,"%g %g %g %g\n",n0->loc.x,n0->loc.y,n0->loc.z,dot(trinorm,vz));
         fprintf(debug_out,"%g %g %g %g\n",n0->loc.x,n0->loc.y,n0->loc.z,dot(trinorm,vz));
      }


      // put the value on the grid by subdivision
      for (k=0; k<num_layers; k++) {

      norm_disp = (double)(2*k+1)/(double)(num_layers) - 1.0;

      for (i=0; i<subdivisions; i++) {

         for (j=0; j<(2*i+1); j++) {

            if (j%2 == 0) {
               A = 3*(subdivisions-i)-2;
               B = (i-j/2)*3+1;
               C = D-A-B;
            } else {
               A = 3*(subdivisions-i)-1;
               B = (i-(j+1)/2)*3+2;
               C = D-A-B;
            }
            ec.x = (A*n0->loc.x+B*n1->loc.x+C*n2->loc.x)/(double)(D);
            ec.y = (A*n0->loc.y+B*n1->loc.y+C*n2->loc.y)/(double)(D);
            ec.z = (A*n0->loc.z+B*n1->loc.z+C*n2->loc.z)/(double)(D);
            // fprintf(stderr,"   point at %g %g %g\n",ec.x,ec.y,ec.z);

            // and perturb it in the triangle-normal direction
            ec.x += trinorm.x*norm_disp;
            ec.y += trinorm.y*norm_disp;
            ec.z += trinorm.z*norm_disp;

            // find location in image coordinates
            xpos = dot(vx,ec) - xmin;
            ypos = dot(vy,ec) - ymin;

#ifdef USE_GAUSSIAN
            // circle of confusion radius in image units (sigma)
            zpos = dot(vz,ec) - zmin;
            rad = fabs(zpos-0.45) + dd;
            cnst = 1./pow(rad,2);
            //fprintf(stderr,"%g %g %g  %g\n",xpos,ypos,zpos,rad);
            if (rad < 0.) {
               // change to 5., run at 10x res to make smooth dots
            } else if (rad < 0.03) {
               // one sample per pixel
               sx = (int)((xpos-3.*rad)/dd);
               if (sx<0) sx=0;
               ex = (int)((xpos+3.*rad)/dd);
               if (ex>=xres) ex=xres-1;
               sy = (int)((ypos-3.*rad)/dd);
               if (sy<0) sy=0;
               ey = (int)((ypos+3.*rad)/dd);
               if (ey>=yres) ey=yres-1;
               //fprintf(stderr,"%d:%d %d:%d  %g %g %g  %g\n",sx,ex,sy,ey,xpos,ypos,zpos,rad);
               for (ix=sx; ix<ex; ix++) {
                  dx = xpos - ix*dd;	// distance in image units
                  drr = pow(dx,2);
                  for (iy=sy; iy<ey; iy++) {
                     dy = ypos - iy*dd;
                     dr = drr + pow(dy,2);
                     //fprintf(stderr,"  %d %d  %g %g  %g\n",ix,iy,dx,dy,dr);
                     a[ix][iy] += factor*cnst*exp(-0.5*dr*cnst);
                  }
               }
                     //exit(0);
            } else {
               // nothing (too diffuse)
            }
#else

            // find lower-left pixel coordinate (be able to accept negative quantities)
            xloc = (int)(floor(xpos/dd));
            yloc = (int)(floor(ypos/dd));
            // fprintf(stderr,"   which is %g %g, cell %d %d\n",xpos,ypos,xloc,yloc);

            // reset xpos, ypos as local cell coordinates
            xpos = xpos/dd - xloc;
            ypos = ypos/dd - yloc;
            // if (xpos < 0 || ypos < 0)
               // fprintf(stderr,"   subcell coords %g %g, weight %g\n",xpos,ypos,factor);

            // and find the elevation, in image coords
            zpos = dot(vz,ec) - zmin;

            // finally, set the density of this point
            if (is_solid) {
               //if (zpos < 0. && is_solid) fprintf(stderr,"\nERROR, zpos %g",zpos);
               rfactor = factor*zpos;
            } else if (FALSE) {
               if (fabs(10.*zpos-floor(10.*zpos)-0.1) < 0.05) {
                 rfactor = factor*(1.+cos(M_PI*20.*(10.*zpos-floor(10.*zpos)-0.5)));
               } else {
                 rfactor = 0.;
               }
            } else {
               rfactor = factor;
            }

            // write a little blob at each point, use area weighting
            if (xloc > -1 && xloc < xres) {
               if (yloc > -1 && yloc < yres)
                  a[xloc][yloc]     += rfactor*(1.0-xpos)*(1.0-ypos);
               if (yloc > -2 && yloc+1 < yres)
                  a[xloc][yloc+1]   += rfactor*(1.0-xpos)*(ypos);
            }
            if (xloc > -2 && xloc+1 < xres) {
               if (yloc > -1 && yloc < yres)
                  a[xloc+1][yloc]   += rfactor*(xpos)    *(1.0-ypos);
               if (yloc > -2 && yloc+1 < yres)
                  a[xloc+1][yloc+1] += rfactor*(xpos)    *(ypos);
            }
#endif
            /*
            if (xloc > -2 && xloc+1 < xres && yloc > -2 && yloc+1 < yres) {
            if (isnan(a[xloc+1][yloc+1])) {
               fprintf(stderr,"\na[%d][%d] was %lf, factor is %lf, pos is %g %g\n",xloc+1,yloc+1,a[xloc+1][yloc+1],factor,xpos,ypos);
               fprintf(stderr,"  ec is %g %g %g, normdisp is %g\n",ec.x,ec.y,ec.z,norm_disp);
               fprintf(stderr,"  trinorm is %g %g %g, area is %g\n",trinorm.x,trinorm.y,trinorm.z,area);
            }
            }
            */
         }
      }
      }

      if (cnt/DOTPER == (cnt+DPMO)/DOTPER) {
         fprintf(stderr,".");
         fflush(stderr);
      }
      cnt++;
      this_tri = this_tri->next_tri;
   }
   fprintf(stderr,"\n");

   if (debug_write)
     fclose(debug_out);

   // finally, print the image

   if (write_pgm) {

   fprintf(stderr,"Writing PGM image"); fflush(stderr);
   // gamma-correct and check for peak value
   for (j=yres-1; j>-1; j--) {
      for (i=0; i<xres; i++) {
         // printval = (int)(a[i][j]);
         // if (isnan(a[i][j])) fprintf(stderr,"\na[%d][%d] was %lf\n",i,j,a[i][j]);
         // if (a[i][j]<0) fprintf(stderr,"\na[%d][%d] was %lf\n",i,j,a[i][j]);
         // this is gamma correction, leave it in, it helps lower-res images
         // um, it doesn't actually do anything? Okay.
         //a[i][j] = exp(0.5*log(a[i][j]));
         if (a[i][j] > maxval) maxval = a[i][j];
      }
   }
   fprintf(stderr,", maxval is %g\n",maxval); fflush(stderr);

   // peak cropping is now a command-line option
   //if (!is_solid) maxval *= 0.8;
   if (!is_solid) maxval *= peak_crop;

   // scale all values
   if (write_hibit) {
      for (j=yres-1; j>-1; j--) {
         for (i=0; i<xres; i++) {
            a[i][j] *= 65536.0/maxval;
            // if (printval > maxval) maxval = printval;
         }
      }
      // write header
      fprintf(stdout,"P2\n%d %d\n%d\n",xres,yres,65535);
      // write data
      for (j=yres-1; j>-1; j--) {
         for (i=0; i<xres; i++) {
            // printval = (unsigned short)(a[i][j]);
            printval = (int)(a[i][j]);
            if (printval > 65535) printval = 65535;
            fprintf(stdout,"%d\n",printval);
         }
      }
   } else {
      for (j=yres-1; j>-1; j--) {
         for (i=0; i<xres; i++) {
            a[i][j] *= 256.0/maxval;
            // if (printval > maxval) maxval = printval;
         }
      }
      // write header
      fprintf(stdout,"P2\n%d %d\n%d\n",xres,yres,255);
      // write data
      for (j=yres-1; j>-1; j--) {
         for (i=0; i<xres; i++) {
            // printval = (unsigned short)(a[i][j]);
            printval = (int)(a[i][j]);
            if (printval > 255) printval = 255;
            fprintf(stdout,"%d\n",printval);
         }
      }
   }

   } else if (write_png) {

   fprintf(stderr,"Writing PNG image"); fflush(stderr);
   // gamma-correct and check for peak value
   for (j=yres-1; j>-1; j--) {
      for (i=0; i<xres; i++) {
         a[i][j] = exp(0.5*log(a[i][j]));
         if (a[i][j] > maxval) maxval = a[i][j];
      }
   }
   fprintf(stderr,", maxval is %g\n",maxval); fflush(stderr);
   if (!is_solid) maxval *= 0.8;
   // scale all values
   if (write_hibit) {
      for (j=yres-1; j>-1; j--) {
         for (i=0; i<xres; i++) {
            printval = (int)(a[i][j]*65536.0/maxval);
            if (printval<0) printval = 0;
            if (printval>65535) printval = 65535;
            img[yres-1-j][2*i] = (png_byte)(printval/256);
            img[yres-1-j][2*i+1] = (png_byte)(printval%256);
         }
      }
      write_png_image(img,yres,xres,16);
   } else {
      for (j=yres-1; j>-1; j--) {
         for (i=0; i<xres; i++) {
            printval = (int)(a[i][j]*256.0/maxval);
            if (printval<0) printval = 0;
            if (printval>255) printval = 255;
            img[yres-1-j][i] = (png_byte)printval;
         }
      }
      write_png_image(img,yres,xres,8);
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
int write_png_image(png_byte** image,int xres,int yres,int depth) {

   png_uint_32 height,width;
   FILE *fp;
   png_structp png_ptr;
   png_infop info_ptr;
   //png_colorp palette;
   //png_voidp user_error_ptr;
   // char *file_name = "out.png";

   height=yres;
   width=xres;

   /* open the file */
   // fp = fopen(file_name, "wb");
   fp = stdout;
   if (fp == NULL)
      return (-1);

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

   /* close the file */
   // fclose(fp);

   /* that's it */
   return (0);
}


