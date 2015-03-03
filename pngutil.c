/*
 *  pngutil.c - Utility subroutines for use primarily with rockdetail
 *
 *  Mark J. Stock, mstock@umich.edu
 *
 * rocktools - Tools for creating and manipulating triangular meshes
 * Copyright (C) 2004-12  Mark J. Stock
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
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#define png_infopp_NULL (png_infopp)NULL
#define int_p_NULL (int*)NULL
#include "png.h"
#include "structs.h"

float** read_png (char*, float, float, int*, int*);
png_byte** allocate_2d_array_pb (int,int,int);
int free_2d_array_pb (png_byte**);
tri_pointer generate_heightmesh (tri_pointer, float**, int, int, int, int, double, int, int, double, double, int);
double find_optimum_offset_z (int, int, float**, int, int, double, double);
tri_pointer add_tris_with_nodes (tri_pointer, int*, node_ptr, node_ptr, node_ptr);
tri_pointer add_tris_with_nodes_tcs (tri_pointer, int*, node_ptr, node_ptr, node_ptr, text_ptr, text_ptr, text_ptr);
tri_pointer add_twotris_by_fournodes (tri_pointer, int*, node_ptr, node_ptr, node_ptr, node_ptr);
tri_pointer add_manytris_by_fournodes (tri_pointer, int*, int*, BIN*, node_ptr, node_ptr, node_ptr, node_ptr);
unsigned char** prepare_flag_array (int, int, int, int, int, int);
unsigned char** allocate_2d_array_uc (int,int);
int fillet_bottom_heights (unsigned char**, float**, int, int, double, double);


/*
 * read a PNG, scale it, save it to 1 channel
 */
float** read_png (char *infile, float redmin, float redrange, int* nx, int* ny) {

   int i,j;
   int high_depth;
   FILE *fp;
   unsigned char header[8];
   float **red;
   png_uint_32 height,width;
   int bit_depth,color_type,interlace_type;
   png_structp png_ptr;
   png_infop info_ptr;
   png_byte **img;

   // check the file
   fp = fopen(infile,"rb");
   if (fp==NULL) {
      fprintf(stderr,"Could not open input file %s\n",infile);
      fflush(stderr);
      exit(0);
   }

   // check to see that it's a PNG
   fread (&header, 1, 8, fp);
   if (png_sig_cmp(header, 0, 8)) {
      fprintf(stderr,"File %s is not a PNG\n",infile);
      fflush(stderr);
      exit(0);
   }

   /* Create and initialize the png_struct with the desired error handler
    * functions.  If you want to use the default stderr and longjump method,
    * you can supply NULL for the last three parameters.  We also supply the
    * the compiler header file version, so that we know if the application
    * was compiled with a compatible version of the library.  REQUIRED
    */
   png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING,
      NULL, NULL, NULL);

   /* Allocate/initialize the memory for image information.  REQUIRED. */
   info_ptr = png_create_info_struct(png_ptr);
   if (info_ptr == NULL) {
      fclose(fp);
      png_destroy_read_struct(&png_ptr, png_infopp_NULL, png_infopp_NULL);
      exit(0);
   }

   /* Set error handling if you are using the setjmp/longjmp method (this is
    * the normal method of doing things with libpng).  REQUIRED unless you
    * set up your own error handlers in the png_create_read_struct() earlier.  */
   if (setjmp(png_jmpbuf(png_ptr))) {
      /* Free all of the memory associated with the png_ptr and info_ptr */
      png_destroy_read_struct(&png_ptr, &info_ptr, png_infopp_NULL);
      fclose(fp);
      /* If we get here, we had a problem reading the file */
      exit(0);
   }

   /* One of the following I/O initialization methods is REQUIRED */
   /* Set up the input control if you are using standard C streams */
   png_init_io(png_ptr, fp);

   /* If we have already read some of the signature */
   png_set_sig_bytes(png_ptr, 8);

   /* The call to png_read_info() gives us all of the information from the
    * PNG file before the first IDAT (image data chunk).  REQUIRED */
   png_read_info(png_ptr, info_ptr);

   png_get_IHDR(png_ptr, info_ptr, &width, &height, &bit_depth, &color_type,
       &interlace_type, int_p_NULL, int_p_NULL);

   /* Set up the data transformations you want.  Note that these are all
    * optional.  Only call them if you want/need them.  Many of the
    * transformations only work on specific types of images, and many
    * are mutually exclusive.  */

   /* tell libpng to strip 16 bit/color files down to 8 bits/color */
   //png_set_strip_16(png_ptr);

   /* Extract multiple pixels with bit depths of 1, 2, and 4 from a single
    * byte into separate bytes (useful for paletted and grayscale images).  */
   png_set_packing(png_ptr);

   /* Expand paletted colors into true RGB triplets */
   //if (color_type == PNG_COLOR_TYPE_PALETTE)
   //   png_set_palette_rgb(png_ptr);

   /* Expand grayscale images to the full 8 bits from 1, 2, or 4 bits/pixel */
   if (color_type == PNG_COLOR_TYPE_GRAY && bit_depth < 8)
      png_set_expand_gray_1_2_4_to_8(png_ptr);

   /* Optional call to gamma correct and add the background to the palette
    * and update info structure.  REQUIRED if you are expecting libpng to
    * update the palette for you (ie you selected such a transform above).
    */
   //png_read_update_info(png_ptr, info_ptr);

   // check image type for applicability
   if (bit_depth != 8 && bit_depth != 16) {
     fprintf(stderr,"INCOMPLETE: read_png expect 8-bit or 16-bit images\n");
     fprintf(stderr,"   bit_depth: %d\n",bit_depth);
     fprintf(stderr,"   file: %s\n",infile);
     exit(0);
   }
   if (color_type != PNG_COLOR_TYPE_GRAY && color_type != PNG_COLOR_TYPE_RGB) {
     fprintf(stderr,"INCOMPLETE: read_png expect grayscale or RGB images\n");
     fprintf(stderr,"   color_type: %d\n",color_type);
     fprintf(stderr,"   file: %s\n",infile);
     exit(0);
   }

   // set channels
   if (color_type != PNG_COLOR_TYPE_GRAY) {
     fprintf(stderr,"ERROR: not expecting 3-channel PNG, but input is 3-channel\n");
     fprintf(stderr,"  file (%s)",infile);
     fprintf(stderr,"  Convert file to grayscale and try again.\n");
     exit(0);
   }

   // set specific bit depth
   if (bit_depth == 16) high_depth = TRUE;
   else high_depth = FALSE;

   // allocate the space for the image array
   img = allocate_2d_array_pb(width,height,bit_depth);

   /* Now it's time to read the image.  One of these methods is REQUIRED */
   png_read_image(png_ptr, img);

   /* read rest of file, and get additional chunks in info_ptr - REQUIRED */
   png_read_end(png_ptr, info_ptr);

   /* At this point you have read the entire image */

   /* clean up after the read, and free any memory allocated - REQUIRED */
   png_destroy_read_struct(&png_ptr, &info_ptr, png_infopp_NULL);

   /* close the file */
   fclose(fp);

   // allocate space for the 2D array of floats
   red = allocate_2d_array_f(width,height);

   // monochrome image, read data from red array

   // no scaling, 16-bit per channel
   if (high_depth) {
      for (j=height-1; j>=0; j--) {
         for (i=0; i<width; i++) {
            red[i][j] = redmin+redrange*(img[height-1-j][2*i]*256+img[height-1-j][2*i+1])/65535.;
         }
      }

   // no scaling, 8-bit per channel
   } else {
      for (j=height-1; j>=0; j--) {
         for (i=0; i<width; i++) {
            red[i][j] = redmin+redrange*img[height-1-j][i]/255.;
         }
      }
   }

   // free the data array
   free_2d_array_pb(img);

   // set the sizes so that we can understand them
   (*nx) = width;
   (*ny) = height;

   return(red);
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

int free_2d_array_pb(png_byte** array){
   free(array[0]);
   free(array);
   return(0);
}


/*
 * generate_heightmesh
 *
 * do the heavy lifting to generate a trimesh of a heightfield
 *
 * inputs:
 *   tri_head - add all triangles here
 *   hf       - 2d array of floats
 *   nx,ny    - dimensions of hf
 *   thick    - thickness in world units
 *   nx,ny    - dimensions of hf
 */
tri_pointer generate_heightmesh (tri_pointer tri_head, float **hf, int nx, int ny,
                                 int do_bottom, int do_trans, double depth,
                                 int do_legs, int do_walls, double thick, double inset,
                                 int do_texture_coords) {

   int i,j,k,l;
   int num_tri = 0;
   int num_nodes = 0;
   int num_texts = 0;
   int rotate_for_curvature = TRUE;
   int need_side;
   int ithick,iinset;
   float **bf = NULL;
   double hscale = 1.0;
   double meanelev,thiselev,dx,dy,ixx,iyy,ixy,tra; //,eigv1,eigv2;
   unsigned char **legflag = NULL;
   VEC nmin,nmax,location;
   UV tc;
   node_ptr nodell,nodelr,nodeur,nodeul,temp_node;
   tri_pointer new_tri = NULL;
   BIN nodebin;
   TBIN textbin;

   // find x and y scales
   if (nx > ny) {
      hscale = 1./(double)(nx-1);
   } else {
      hscale = 1./(double)(ny-1);
   }

   // how many cells thick are the walls/legs?
   ithick = (int)(thick/hscale + 0.8);
   if (ithick < 1) ithick = 1;
   //fprintf(stderr,"thick %g %g %d\n",thick,thick/hscale,ithick);
   if (do_legs) fprintf(stderr,"  legs are %d cells thick\n",ithick);
   if (do_walls) fprintf(stderr,"  walls are %d cells thick\n",ithick);
   iinset = (int)(inset/hscale + 0.8);
   if (do_walls || do_legs) fprintf(stderr,"  and inset %d cells\n",iinset);
   //fprintf(stderr,"inset %g %g %d\n",inset,inset/hscale,iinset);

   // prepare the logical array of legs/walls
   legflag = prepare_flag_array (do_walls,do_legs,nx-1,ny-1,ithick,iinset);

   // prepare the array of bottom heights
   if (do_bottom) {
      fprintf(stderr,"computing bottom surface elevations...\n"); fflush(stderr);
      // make space for the array
      bf = allocate_2d_array_f(nx,ny);
      // set the values so that we can re-use them

      if (do_trans) {
         // simple: constant thickness, vertically
         for (i=0; i<nx; i++) {
            for (j=0; j<ny; j++) {
               bf[i][j] = -hf[i][j];
            }
         }

      } else {
         // complex: never allow thickness to be less than "depth"
         for (i=0; i<nx; i++) {
            for (j=0; j<ny; j++) {
               bf[i][j] = find_optimum_offset_z(i,j,hf,nx,ny,hscale,depth);
            }
         }
      }

      // additional: add fillets around walls and legs
      if (TRUE) {
         fprintf(stderr,"computing bottom surface fillets...\n"); fflush(stderr);
         fillet_bottom_heights (legflag,bf,nx,ny,depth,hscale);
      }
   }


   // bin data
   nmax.x = hscale*(nx-0.5);
   nmax.y = hscale*(ny-0.5);
   nmax.z = 0.;
   nmin.x = hscale*-0.5;
   nmin.y = hscale*-0.5;
   nmin.z = 0.01;

   // Initialize bin structure
   (void) prepare_node_bin (&nodebin,nmin,nmax);
   (void) prepare_texture_bin (&textbin);

   fprintf(stderr,"creating triangles...\n"); fflush(stderr);
   for (i=0; i<nx-1; i++) {
      for (j=0; j<ny-1; j++) {

         // always write the top ---------------------------------------

         // lower left corner
         location.x = hscale*(double)i;
         location.y = hscale*(double)j;
         location.z = (double)hf[i][j];
         //fprintf(stdout,"%d %d  %g %g %g\n",i,j,location.x,location.y,location.z);
         //fflush(stdout);
         // if we don't need connectivity info, this call ignores tri1 and 0
         nodell = add_to_nodes_list(new_tri,&num_nodes,0,&location,&nodebin);
         tc.x = location.x;
         tc.y = location.y;
         text_ptr textll = add_to_textures_list(&num_texts,&tc,&textbin);

         // lower right corner
         location.x = hscale*(double)(i+1);
         location.y = hscale*(double)j;
         location.z = (double)hf[i+1][j];
         nodelr = add_to_nodes_list(new_tri,&num_nodes,0,&location,&nodebin);
         tc.x = location.x;
         tc.y = location.y;
         text_ptr textlr = add_to_textures_list(&num_texts,&tc,&textbin);

         // upper right corner
         location.x = hscale*(double)(i+1);
         location.y = hscale*(double)(j+1);
         location.z = (double)hf[i+1][j+1];
         nodeur = add_to_nodes_list(new_tri,&num_nodes,0,&location,&nodebin);
         tc.x = location.x;
         tc.y = location.y;
         text_ptr textur = add_to_textures_list(&num_texts,&tc,&textbin);

         // upper left corner
         location.x = hscale*(double)i;
         location.y = hscale*(double)(j+1);
         location.z = (double)hf[i][j+1];
         nodeul = add_to_nodes_list(new_tri,&num_nodes,0,&location,&nodebin);
         //fprintf(stdout,"%d %d  %g %g %g\n",i,j,location.x,location.y,location.z);
         tc.x = location.x;
         tc.y = location.y;
         text_ptr textul = add_to_textures_list(&num_texts,&tc,&textbin);

         // curvature estimation
         //if (i==1 && j==1) {
         if (rotate_for_curvature) {
            // find mean
            meanelev = 0.25*(hf[i][j]+hf[i+1][j]+hf[i][j+1]+hf[i+1][j+1]);
            //fprintf(stdout,"mean  %g\n",meanelev);
            ixx = 0.;
            iyy = 0.;
            ixy = 0.;
            //meanelev = 0.;
            for (k=-1; k<3; k++) {
               dx = k-0.5;
               if (i+k >= 0 && i+k < nx) {
               for (l=-1; l<3; l++) {
                  dy = l-0.5;
                  if (j+l >= 0 && j+l < ny) {
                     thiselev = hf[i+k][j+l] - meanelev;
                     //thiselev = -fabs((float)(dx-dy));
                     //fprintf(stdout," %g",hf[i+k][j+l]);
                     //fprintf(stdout," %g",thiselev);
                     ixx += thiselev * dx*dx;
                     ixy += thiselev * dx*dy;
                     iyy += thiselev * dy*dy;
                  }
               }
               }
               //fprintf(stdout,"\n");
            }
            //fprintf(stdout,"%d %d  %g\n",i,j,ixy);
            //fprintf(stdout,"  %g %g\n",ixx,ixy);
            //fprintf(stdout,"  %g %g\n",ixy,iyy);
            tra = ixx+iyy;
            //eigv1 = 0.5 * (-tra - sqrt(tra*tra - 4.*(ixx*iyy-ixy*ixy)));
            //eigv2 = 0.5 * (-tra + sqrt(tra*tra - 4.*(ixx*iyy-ixy*ixy)));
            //fprintf(stdout,"invariants  %g %g\n",tra,ixx*iyy-ixy*ixy);
            //fprintf(stdout,"eigs  %g %g\n",eigv1,eigv2);
            // if both eigs are near zero, then it's near flat
            // if both eigs are positive, it's a local maximum
            // if both eigs are negative, it's a local minimum
            // if eigv1=-eigv2, it's a perfect saddle point
            // forget all that shit, if ixy is negative, cross one way, positive the other
            //exit(0);

            // flip diagonal if needed (by rotating nodes)
            if (ixy*tra > 0.0) {
               temp_node = nodell;
               nodell = nodelr;
               nodelr = nodeur;
               nodeur = nodeul;
               nodeul = temp_node;
            }
         }

         // make the two tris
         if (do_texture_coords) {
            tri_head = add_tris_with_nodes_tcs (tri_head,&num_tri,nodell,nodelr,nodeur,textll,textlr,textur);
            tri_head = add_tris_with_nodes_tcs (tri_head,&num_tri,nodell,nodeur,nodeul,textll,textur,textul);
         } else {
            tri_head = add_tris_with_nodes (tri_head,&num_tri,nodell,nodelr,nodeur);
            tri_head = add_tris_with_nodes (tri_head,&num_tri,nodell,nodeur,nodeul);
         }
         //tri_head = add_twotris_by_fournodes (tri_head,&num_tri,nodell,nodelr,nodeur,nodeul);


         // sometimes write the bottom ---------------------------------
         if (do_bottom) {

            // lower left corner
            location.x = hscale*(double)i;
            location.y = hscale*(double)j;
            location.z = (double)bf[i][j];
            if (legflag[i][j] != 0) location.z = 0.0;
            nodell = add_to_nodes_list(new_tri,&num_nodes,0,&location,&nodebin);
            // lower right corner
            location.x = hscale*(double)(i+1);
            location.y = hscale*(double)j;
            location.z = (double)bf[i+1][j];
            if (legflag[i][j] != 0) location.z = 0.0;
            nodelr = add_to_nodes_list(new_tri,&num_nodes,0,&location,&nodebin);
            // upper right corner
            location.x = hscale*(double)(i+1);
            location.y = hscale*(double)(j+1);
            location.z = (double)bf[i+1][j+1];
            if (legflag[i][j] != 0) location.z = 0.0;
            nodeur = add_to_nodes_list(new_tri,&num_nodes,0,&location,&nodebin);
            // upper left corner
            location.x = hscale*(double)i;
            location.y = hscale*(double)(j+1);
            location.z = (double)bf[i][j+1];
            if (legflag[i][j] != 0) location.z = 0.0;
            nodeul = add_to_nodes_list(new_tri,&num_nodes,0,&location,&nodebin);

            // flip diagonal if needed (by rotating nodes)
            if (rotate_for_curvature) {
               // find mean
               meanelev = 0.25*(bf[i][j]+bf[i+1][j]+bf[i][j+1]+bf[i+1][j+1]);
               ixx = 0.;
               iyy = 0.;
               ixy = 0.;
               for (k=-1; k<3; k++) {
                  dx = k-0.5;
                  if (i+k >= 0 && i+k < nx) {
                  for (l=-1; l<3; l++) {
                     dy = l-0.5;
                     if (j+l >= 0 && j+l < ny) {
                        thiselev = bf[i+k][j+l] - meanelev;
                        ixx += thiselev * dx*dx;
                        ixy += thiselev * dx*dy;
                        iyy += thiselev * dy*dy;
                        }
                     }
                  }
               }
               tra = ixx+iyy;
               if (ixy*tra > 0.0) {
                  temp_node = nodell;
                  nodell = nodelr;
                  nodelr = nodeur;
                  nodeur = nodeul;
                  nodeul = temp_node;
               }
            }

            // make the two tris
            tri_head = add_tris_with_nodes (tri_head,&num_tri,nodell,nodeul,nodeur);
            tri_head = add_tris_with_nodes (tri_head,&num_tri,nodell,nodeur,nodelr);
            //tri_head = add_twotris_by_fournodes (tri_head,&num_tri,nodell,nodeul,nodeur,nodelr);

            // does this cell (between 2x2 nodes) need side panels?

            // does it need a left panel?
            need_side = FALSE;
            if (i==0) {
               need_side = TRUE;
            } else if (legflag[i][j] != 0) {
               if (legflag[i-1][j] == 0) need_side = TRUE;
            }
            if (need_side) {
               // lower left corner
               location.x = hscale*(double)i;
               location.y = hscale*(double)(j+1);
               if (legflag[i][j] != 0) {
                  location.z = 0.0;
               } else {
                  location.z = (double)bf[i][j+1];
               }
               nodell = add_to_nodes_list(new_tri,&num_nodes,0,&location,&nodebin);
               // lower right corner
               location.x = hscale*(double)i;
               location.y = hscale*(double)j;
               if (legflag[i][j] != 0) {
                  location.z = 0.0;
               } else {
                  location.z = (double)bf[i][j];
               }
               nodelr = add_to_nodes_list(new_tri,&num_nodes,0,&location,&nodebin);
               // upper right corner
               location.x = hscale*(double)i;
               location.y = hscale*(double)j;
               if (i==0) {
                  location.z = (double)hf[i][j];
               } else {
                  location.z = (double)bf[i][j];
               }
               nodeur = add_to_nodes_list(new_tri,&num_nodes,0,&location,&nodebin);
               // upper left corner
               location.x = hscale*(double)i;
               location.y = hscale*(double)(j+1);
               if (i==0) {
                  location.z = (double)hf[i][j+1];
               } else {
                  location.z = (double)bf[i][j+1];
               }
               nodeul = add_to_nodes_list(new_tri,&num_nodes,0,&location,&nodebin);
               // make the two tris
               tri_head = add_manytris_by_fournodes (tri_head,&num_tri,&num_nodes,&nodebin,nodell,nodelr,nodeur,nodeul);
            }

            // does it need a lower panel?
            need_side = FALSE;
            if (j==0) {
               need_side = TRUE;
            } else if (legflag[i][j] != 0) {
               if (legflag[i][j-1] == 0) need_side = TRUE;
            }
            if (need_side) {
               // lower left corner
               location.x = hscale*(double)i;
               location.y = hscale*(double)j;
               if (legflag[i][j] != 0) {
                  location.z = 0.0;
               } else {
                  location.z = (double)bf[i][j];
               }
               nodell = add_to_nodes_list(new_tri,&num_nodes,0,&location,&nodebin);
               // lower right corner
               location.x = hscale*(double)(i+1);
               location.y = hscale*(double)j;
               if (legflag[i][j] != 0) {
                  location.z = 0.0;
               } else {
                  location.z = (double)bf[i+1][j];
               }
               nodelr = add_to_nodes_list(new_tri,&num_nodes,0,&location,&nodebin);
               // upper right corner
               location.x = hscale*(double)(i+1);
               location.y = hscale*(double)j;
               if (j==0) {
                  location.z = (double)hf[i+1][j];
               } else {
                  location.z = (double)bf[i+1][j];
               }
               nodeur = add_to_nodes_list(new_tri,&num_nodes,0,&location,&nodebin);
               // upper left corner
               location.x = hscale*(double)i;
               location.y = hscale*(double)j;
               if (j==0) {
                  location.z = (double)hf[i][j];
               } else {
                  location.z = (double)bf[i][j];
               }
               nodeul = add_to_nodes_list(new_tri,&num_nodes,0,&location,&nodebin);
               // make the two tris
               tri_head = add_manytris_by_fournodes (tri_head,&num_tri,&num_nodes,&nodebin,nodell,nodelr,nodeur,nodeul);
            }

            // does it need a right panel?
            need_side = FALSE;
            if (i==nx-2) {
               need_side = TRUE;
            } else if (legflag[i][j] != 0) {
               if (legflag[i+1][j] == 0) need_side = TRUE;
            }
            if (need_side) {
               // lower left corner
               location.x = hscale*(double)(i+1);
               location.y = hscale*(double)j;
               if (legflag[i][j] != 0) {
                  location.z = 0.0;
               } else {
                  location.z = (double)bf[i+1][j];
               }
               nodell = add_to_nodes_list(new_tri,&num_nodes,0,&location,&nodebin);
               // lower right corner
               location.x = hscale*(double)(i+1);
               location.y = hscale*(double)(j+1);
               if (legflag[i][j] != 0) {
                  location.z = 0.0;
               } else {
                  location.z = (double)bf[i+1][j+1];
               }
               nodelr = add_to_nodes_list(new_tri,&num_nodes,0,&location,&nodebin);
               // upper right corner
               location.x = hscale*(double)(i+1);
               location.y = hscale*(double)(j+1);
               if (i==nx-2) {
                  location.z = (double)hf[i+1][j+1];
               } else {
                  location.z = (double)bf[i+1][j+1];
               }
               nodeur = add_to_nodes_list(new_tri,&num_nodes,0,&location,&nodebin);
               // upper left corner
               location.x = hscale*(double)(i+1);
               location.y = hscale*(double)j;
               if (i==nx-2) {
                  location.z = (double)hf[i+1][j];
               } else {
                  location.z = (double)bf[i+1][j];
               }
               nodeul = add_to_nodes_list(new_tri,&num_nodes,0,&location,&nodebin);
               // make the two tris
               tri_head = add_manytris_by_fournodes (tri_head,&num_tri,&num_nodes,&nodebin,nodell,nodelr,nodeur,nodeul);
            }

            // does it need an upper panel?
            need_side = FALSE;
            if (j==ny-2) {
               need_side = TRUE;
            } else if (legflag[i][j] != 0) {
               if (legflag[i][j+1] == 0) need_side = TRUE;
            }
            if (need_side) {
               // lower left corner
               location.x = hscale*(double)(i+1);
               location.y = hscale*(double)(j+1);
               if (legflag[i][j] != 0) {
                  location.z = 0.0;
               } else {
                  location.z = (double)bf[i+1][j+1];
               }
               nodell = add_to_nodes_list(new_tri,&num_nodes,0,&location,&nodebin);
               // lower right corner
               location.x = hscale*(double)i;
               location.y = hscale*(double)(j+1);
               if (legflag[i][j] != 0) {
                  location.z = 0.0;
               } else {
                  location.z = (double)bf[i][j+1];
               }
               nodelr = add_to_nodes_list(new_tri,&num_nodes,0,&location,&nodebin);
               // upper right corner
               location.x = hscale*(double)i;
               location.y = hscale*(double)(j+1);
               if (j==ny-2) {
                  location.z = (double)hf[i][j+1];
               } else {
                  location.z = (double)bf[i][j+1];
               }
               nodeur = add_to_nodes_list(new_tri,&num_nodes,0,&location,&nodebin);
               // upper left corner
               location.x = hscale*(double)(i+1);
               location.y = hscale*(double)(j+1);
               if (j==ny-2) {
                  location.z = (double)hf[i+1][j+1];
               } else {
                  location.z = (double)bf[i+1][j+1];
               }
               nodeul = add_to_nodes_list(new_tri,&num_nodes,0,&location,&nodebin);
               // make the two tris
               tri_head = add_manytris_by_fournodes (tri_head,&num_tri,&num_nodes,&nodebin,nodell,nodelr,nodeur,nodeul);
            }

         } // end if do_bottom

      } // end for j=0...
   } // end for i=0...


   // de-allocate hf and bf
   free_2d_array_f(hf);
   if (do_bottom) free_2d_array_f(bf);

   fprintf(stderr,"Nodes: %d\n",num_nodes);
   fprintf(stderr,"Triangles: %d\n",num_tri);

   return(tri_head);
}


/*
 * find the lowest position on the vertical line through i,j, below hf[i][j]
 *   that is exactly depth distance from the heightfield
 */
double find_optimum_offset_z (int i, int j, float** hf, int nx, int ny,
                              double hscale, double doffset) {

   int ii,jj;
   int irange;	// check this many rings away from i,j
   int imin,imax;
   int jmin,jmax;
   double lowest,dx,dy,dz;
   double do2 = doffset*doffset;

   // how far out do we need to test?
   irange = (int)(doffset/hscale) + 1;
   //fprintf(stderr,"irange is %d\n",irange);
   //exit(0);

   imin = i-irange;
   if (imin < 0) imin = 0;
   imax = i+irange+1;
   if (imax > nx) imax = nx;
   jmin = j-irange;
   if (jmin < 0) jmin = 0;
   jmax = j+irange+1;
   if (jmax > ny) jmax = ny;

   // initialize the lowest point here
   lowest = hf[i][j] - doffset;

   // loop over possible ranges, looking for lower points
   for (ii=imin; ii<imax; ii++) {
      dx = pow(hscale*(ii-i),2);
      for (jj=jmin; jj<jmax; jj++) {
         dy = pow(hscale*(jj-j),2);
         //dz = pow(hf[ii][jj]);
         dz = pow(hf[ii][jj]-lowest,2);
         //dh = sqrt(dy);
         // is this node (ii,jj) too close to lowest?
         //if (dx+dy+dz < do2) {
            //fprintf(stderr,"  %d %d   %d %d  %g  %g  %g\n",i,j,ii,jj,sqrt(do2),sqrt(dz),sqrt(dx+dy+dz));
            // move the "lowest" point down some
            //lowest = hf[ii][jj] - sqrt(do2-dx-dy);
         //}
         if (dx+dy < do2) {
            dz = hf[ii][jj] - sqrt(do2-dx-dy);
            if (dz < lowest) lowest = dz;
         }
         //fprintf(stderr,"  %d %d   %d %d  %g %g %g  %g  %g\n",i,j,ii,jj,sqrt(dx),sqrt(dy),sqrt(dz),sqrt(dx+dy+dz),sqrt(do2));
      }
   }

   //return (hf[i][j] - doffset);
   return (lowest);
}


/*
 * thicken the mesh somewhat near the walls and legs
 */
int fillet_bottom_heights (unsigned char **legflag, float **bf,
                           int nx, int ny, double depth, double hscale) {

   // how many iterations of the diffusion?
   int numits = (int)(depth/hscale) + 1;
   int iter,i,j;
   float **df1 = NULL;
   float **df2 = NULL;
   float **dft = NULL;

   // make a temporary array of values to diffuse
   df1 = allocate_2d_array_f(nx-1,ny-1);
   df2 = allocate_2d_array_f(nx-1,ny-1);

   // first, set an elevation for all cells
   for (i=0; i<nx-1; i++) {
      for (j=0; j<ny-1; j++) {
         if (legflag[i][j] != 0) df1[i][j] = 1.0;
         else df1[i][j] = 0.0;
      }
   }

   // first, count the distance from a leg/wall
   // do this with a simple diffusion process
   for (iter=0; iter<numits; iter++) {
      // run diffusion over the whole grid
      for (i=1; i<nx-2; i++) {
         for (j=1; j<ny-2; j++) {
            //df2[i][j] = df1[i][j] + 0.125 * (df1[i-1][j]+df1[i-1][j]+df1[i][j-1]+df1[i][j+1]-4.*df1[i][j]);
            df2[i][j] = 0.5*df1[i][j] + 0.125*(df1[i-1][j]+df1[i+1][j]+df1[i][j-1]+df1[i][j+1]);
         }
      }
      // and copy back to first array
      dft = df2;
      df2 = df1;
      df1 = dft;
   }
   free_2d_array_f(df2);

   // then, adjust the elevation "bf" for each cell
   for (i=1; i<nx-1; i++) {
      for (j=1; j<ny-1; j++) {
         bf[i][j] -= depth * 0.25*(df1[i-1][j-1]+df1[i-1][j]+df1[i][j-1]+df1[i][j]);
      }
   }
   free_2d_array_f(df1);

   return (0);
}


/*
 * make a flag array, 0=no leg, 1=leg, 2..10 slightly shorter legs (eases simplification)
 * gah! just click on "planar simplification" in meshlab->remesh->QED
 */
unsigned char** prepare_flag_array (int do_walls, int do_legs,
                                    int nx, int ny,
                                    int ithick, int iinset) {

   int i,j;
   unsigned char** flag;

   // allocate the flag array
   flag = allocate_2d_array_uc (nx,ny);

   for (i=0; i<nx; i++) {
      for (j=0; j<ny; j++) {
         // default: no leg here
         flag[i][j] = 0;
      }
   }

   // set wall points
   if (do_walls) {
      i = iinset;
      for (; i<iinset+ithick; i++) {
         for (j=iinset; j<ny-iinset; j++) flag[i][j] = 1;
      }
      for (; i<nx-iinset-ithick; i++) {
         for (j=iinset; j<iinset+ithick; j++) flag[i][j] = 1;
         for (j=ny-iinset-ithick; j<ny-iinset; j++) flag[i][j] = 1;
      }
      for (; i<nx-iinset; i++) {
         for (j=iinset; j<ny-iinset; j++) flag[i][j] = 1;
      }
   }

   // set leg points
   if (do_legs) {
      for (i=iinset; i<iinset+ithick; i++) {
         for (j=iinset; j<iinset+ithick; j++) flag[i][j] = 1;
         for (j=ny-iinset-ithick; j<ny-iinset; j++) flag[i][j] = 1;
      }
      for (i=nx-iinset-ithick; i<nx-iinset; i++) {
         for (j=iinset; j<iinset+ithick; j++) flag[i][j] = 1;
         for (j=ny-iinset-ithick; j<ny-iinset; j++) flag[i][j] = 1;
      }
   }

   return (flag);
}


/*
 * use this often: add a tri
 *
 * always draw a pair of triangles with the diag from ll to ur
 */
tri_pointer add_tris_with_nodes (tri_pointer tri_head, int* num_tri,
                                 node_ptr n1, node_ptr n2, node_ptr n3) {
   // initialize a tri
   tri_pointer new_tri = alloc_new_tri();
   new_tri->index = (*num_tri);
   (*num_tri) = (*num_tri)+1;

   // set the node pointers
   new_tri->node[0] = n1;
   new_tri->node[1] = n2;
   new_tri->node[2] = n3;

   // add it on as the new head of the list
   if (tri_head) {
      new_tri->next_tri = tri_head;
      tri_head = new_tri;
   } else {
      tri_head = new_tri;
      tri_head->next_tri = NULL;
   }

   return (tri_head);
}

tri_pointer add_tris_with_nodes_tcs (tri_pointer tri_head, int* num_tri,
                                     node_ptr n1, node_ptr n2, node_ptr n3,
                                     text_ptr t1, text_ptr t2, text_ptr t3) {
   // initialize a tri
   tri_pointer new_tri = alloc_new_tri();
   new_tri->index = (*num_tri);
   (*num_tri) = (*num_tri)+1;

   // set the node pointers
   new_tri->node[0] = n1;
   new_tri->node[1] = n2;
   new_tri->node[2] = n3;
   new_tri->texture[0] = t1;
   new_tri->texture[1] = t2;
   new_tri->texture[2] = t3;

   // add it on as the new head of the list
   if (tri_head) {
      new_tri->next_tri = tri_head;
      tri_head = new_tri;
   } else {
      tri_head = new_tri;
      tri_head->next_tri = NULL;
   }

   return (tri_head);
}

tri_pointer add_twotris_by_fournodes (tri_pointer tri_head, int* num_tri,
                              node_ptr ll, node_ptr lr,
                              node_ptr ur, node_ptr ul) {
   // initialize a tri
   tri_pointer new_tri = alloc_new_tri();
   new_tri->index = (*num_tri);
   (*num_tri) = (*num_tri)+1;

   // set the node pointers
   new_tri->node[0] = ll;
   new_tri->node[1] = lr;
   new_tri->node[2] = ur;
   // add it on as the new head of the list
   if (tri_head) {
      new_tri->next_tri = tri_head;
      tri_head = new_tri;
   } else {
      tri_head = new_tri;
      tri_head->next_tri = NULL;
   }

   // and the other one
   new_tri = alloc_new_tri();
   new_tri->index = (*num_tri);
   (*num_tri) = (*num_tri)+1;
   // set the node pointers
   new_tri->node[0] = ll;
   new_tri->node[1] = ur;
   new_tri->node[2] = ul;
   // add it on as the new head of the list
   if (tri_head) {
      new_tri->next_tri = tri_head;
      tri_head = new_tri;
   } else {
      tri_head = new_tri;
      tri_head->next_tri = NULL;
   }

   return (tri_head);
}

/*
 * use this often: add a list of tris corresponding to a vertical leg
 */
tri_pointer add_manytris_by_fournodes (tri_pointer tri_head, int* num_tri,
                              int* num_nodes, BIN* nodebin,
                              node_ptr ll, node_ptr lr,
                              node_ptr ur, node_ptr ul) {

   static int numCells = 20;
   int i;
   VEC location;
   tri_pointer new_tri = NULL;
   node_ptr left[numCells+1];
   node_ptr right[numCells+1];

   left[0] = ll;
   left[numCells] = ul;
   right[0] = lr;
   right[numCells] = ur;

   // first, make the interim nodes
   for (i=1; i<numCells; i++) {
      double frac = (double)i / (double)(numCells);
      double omfrac = 1.0 - frac;
      location.x = omfrac*ll->loc.x + frac*ul->loc.x;
      location.y = omfrac*ll->loc.y + frac*ul->loc.y;
      location.z = omfrac*ll->loc.z + frac*ul->loc.z;
      left[i] = add_to_nodes_list(new_tri,num_nodes,0,&location,nodebin);
      location.x = omfrac*lr->loc.x + frac*ur->loc.x;
      location.y = omfrac*lr->loc.y + frac*ur->loc.y;
      location.z = omfrac*lr->loc.z + frac*ur->loc.z;
      right[i] = add_to_nodes_list(new_tri,num_nodes,0,&location,nodebin);
   }

   // then, make the pairs of triangles
   for (i=0; i<numCells; i++) {
      tri_head = add_twotris_by_fournodes(tri_head, num_tri, left[i], right[i], right[i+1], left[i+1]);
   }

   return (tri_head);
}


