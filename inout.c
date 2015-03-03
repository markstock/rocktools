/*************************************************************
 *
 *  inout.c - Input and output routines for triangle meshes
 *
 *  Mark J. Stock, mstock@umich.edu
 *
 *
 * rocktools - Tools for creating and manipulating triangular meshes
 * Copyright (C) 1999,2002-4,6,13-14  Mark J. Stock
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
#include <ctype.h>

#include "structs.h"

tri_pointer read_input(char[80],int,tri_pointer);
tri_pointer read_raw(char[80],tri_pointer);
tri_pointer read_tin(char[80],tri_pointer);
tri_pointer read_obj(char[80],int,tri_pointer);
tri_pointer read_msh(char[80],tri_pointer);

int write_output(tri_pointer, char[4], int, int, char**);
int write_raw(tri_pointer, int);
int write_tin(tri_pointer, int);
int write_obj(tri_pointer, int, int, char**);
int write_seg(tri_pointer, int, char**);
int write_rad(tri_pointer, int);
int write_pov(tri_pointer, int);
int write_rib(tri_pointer, int);

int find_mesh_stats(char *,VEC*,VEC*,int*,int*);
int get_tri(FILE*,int,tri_pointer);
int write_tri(FILE*,int,tri_pointer);


/*
 * Determine file type and run appropriate read routine
 */
tri_pointer read_input(char infile[80],int invert,tri_pointer tri_head) {

   //int dummy;
   char extension[4];		/* filename extension if infile */
   tri_pointer new_tri_head;	/* the pointer to the first triangle */

   /* Determine the input file format from the .XXX extension, and read it */
   strncpy(extension,infile+strlen(infile)-3,4);
   if (strncmp(extension, "raw", 3) == 0)
      new_tri_head = read_raw(infile,tri_head);
   else if (strncmp(extension, "obj", 3) == 0)
      new_tri_head = read_obj(infile,invert,tri_head);
   else if (strncmp(extension, "tin", 3) == 0)
      new_tri_head = read_tin(infile,tri_head);
   else if (strncmp(extension, "msh", 3) == 0)
      new_tri_head = read_msh(infile,tri_head);
   else {
      fprintf(stderr,"Input filename extension is assumed to be (%s)\n",extension);
      fprintf(stderr,"This is either not a supported input file format, or you\n");
      fprintf(stderr,"   need to use a proper extension (last 3 characters in filename).\n");
      exit(0);
   }

   // make sure all tris are oriented in the same direction
   //dummy = fix_orientation(new_tri_head);
   //fprintf(stderr,"Flipped %d normals\n",dummy);

   // determine the adjacent triangles for each triangle
#ifdef CONN
   (void) set_adjacent_tris(new_tri_head);
#endif

   return new_tri_head;
}


/*
 * Read in a Wavefront (.obj) file
 *
 * For now, assume that the Wavefront file was written in short
 * form (meaning that every node only appears once, and many
 * triangles may reference the same node
 */
//#pragma GCC push_options
//#pragma GCC optimize ("O0")
tri_pointer __attribute__((optimize("O0"))) read_obj(char filename[80],int invert,tri_pointer tri_head) {

   int inode = 0;
   int inorm = 0;
   int itext = 0;
   int num_tri = 0;
   int num_nodes = 0;
   int num_norms = 0;
   int num_texts = 0;
   char onechar[2],anotherchar[2];
   char twochar[3];
   char sbuf[128];
   char xs[20],ys[20],zs[20];
   VEC *loc = NULL;
   VEC *normal = NULL;
   UV *texture = NULL;
   VEC test,nmin,nmax;
   node_ptr new_node = NULL;
   norm_ptr new_norm = NULL;
   text_ptr new_text = NULL;
   tri_pointer new_tri = NULL;
   BIN nodebin;
   NBIN normbin;
   TBIN textbin;
   FILE *fp;

   // fprintf(stderr,"Try to read %s\n",filename); fflush(stderr);

   // open the .obj file for reading
   fp = fopen(filename,"r");
   if (fp==NULL) {
      fprintf(stderr,"Could not open input file %s\n",filename);
      fflush(stderr);
      exit(0);
   }
   fprintf(stderr,"Prescanning %s",filename);
   fflush(stderr);

   nmax.x = -9.999e+9;
   nmax.y = -9.999e+9;
   nmax.z = -9.999e+9;
   nmin.x = 9.999e+9;
   nmin.y = 9.999e+9;
   nmin.z = 9.999e+9;

   // read the input file and update statistics
   while (fread(&onechar,sizeof(char),1,fp) == 1) {

      // only read nodes, we're putting them into bins on this pass
      if (onechar[0] == 'v') {

         fread(&anotherchar,sizeof(char),1,fp);
         // fprintf(stdout,"  (%c)\n",anotherchar); fflush(stdout);

         if (isspace(anotherchar[0])) {
            // read a vertex location
            fscanf(fp,"%s %s %s",xs,ys,zs);
            //fprintf(stderr,"%d  %s %s %s\n",i,xs,ys,zs); fflush(stderr);
            test.x = atof(xs);
            test.y = atof(ys);
            test.z = atof(zs);
            if (test.x > nmax.x) nmax.x = test.x;
            if (test.y > nmax.y) nmax.y = test.y;
            if (test.z > nmax.z) nmax.z = test.z;
            if (test.x < nmin.x) nmin.x = test.x;
            if (test.y < nmin.y) nmin.y = test.y;
            if (test.z < nmin.z) nmin.z = test.z;
            //fprintf(stderr,"%d  %g %g %g\n",i,test.x,test.y,test.z); fflush(stderr);
            if (inode/DOTPER == (inode+DPMO)/DOTPER) fprintf(stderr,".");
            inode++;
         } else if (anotherchar[0] == 'n') {
            // it's a vertex normal
            inorm++;
         } else if (anotherchar[0] == 't') {
            // it's a texture coordinate
            itext++;
         } else {
            // if its not identifiable, skip it, do not scale, do not write
         }

      } else {
         // if its not identifiable, skip it, do not scale, do not write
      }
      fscanf(fp,"%[^\n]",sbuf);	// read line beyond first char
      fscanf(fp,"%[\n]",twochar);	// read newline
   }
   fprintf(stderr,"\n");
   fclose(fp);

   // now, we know how many nodes there will be, size and malloc the array
   loc = (VEC*)malloc(inode*sizeof(VEC));
   // fprintf(stderr,"Found %d vertexes\n",i);

   if (inorm > 0) normal = (VEC*)malloc(inorm*sizeof(VEC));
   if (itext > 0) texture = (UV*)malloc(itext*sizeof(UV));

   // Initialize bin structures
   // fprintf(stderr,"%g %g %g   %g %g %g\n",nmin.x,nmin.y,nmin.z,nmax.x,nmax.y,nmax.z);
   //nodebin.dx = (nmax.x-nmin.x)/(BIN_COUNT-1);
   //nodebin.start = nmin.x - 0.5*nodebin.dx;
   //for (i=0;i<BIN_COUNT;i++) nodebin.b[i] = NULL;
   (void) prepare_node_bin (&nodebin,nmin,nmax);
   (void) prepare_norm_bin (&normbin);
   (void) prepare_texture_bin (&textbin);


   // open the .obj file for reading
   fp = fopen(filename,"r");
   if (fp==NULL) {
      fprintf(stderr,"Could not open input file %s\n",filename);
      fflush(stderr);
      exit(0);
   }
   fprintf(stderr,"Opening file %s...",filename);
   fflush(stderr);

   // node counter variable!
   inode = 0;
   // node normal counter variable!
   inorm = 0;
   // node texture counter variable!
   itext = 0;
   // line count variable
   int nlines = 0;

   // read the input file and update statistics
   while (fread(&onechar,sizeof(char),1,fp) == 1) {

      nlines++;
      // fprintf(stdout,"(%c)",onechar); fflush(stdout);

      if (onechar[0] == '#') {
         // read a comment line
         // fprintf(stdout,"#%s\n",sbuf);	// write comment

      } else if (onechar[0] == 'v') {

         fread(&anotherchar,sizeof(char),1,fp);
         // fprintf(stdout,"  (%c)\n",anotherchar); fflush(stdout);

         if (isspace(anotherchar[0])) {
            // read a vertex location
            fscanf(fp,"%s %s %s",xs,ys,zs);
            // fprintf(stderr,"%d  (%s) (%s) (%s)\n",i,xs,ys,zs); fflush(stderr);
            loc[inode].x = atof(xs);
            loc[inode].y = atof(ys);
            loc[inode].z = atof(zs);
            inode++;
         } else if (anotherchar[0] == 'n') {
            // read a vertex normal
            fscanf(fp,"%s %s %s",xs,ys,zs);
            // fprintf(stderr,"n %d  %s %s %s\n",k,xs,ys,zs); fflush(stderr);
            normal[inorm].x = atof(xs);
            normal[inorm].y = atof(ys);
            normal[inorm].z = atof(zs);
            inorm++;
         } else if (anotherchar[0] == 't') {
            // read a vertex texture coordinate
            fscanf(fp,"%s %s",xs,ys);
            // fprintf(stderr,"n %d  %s %s %s\n",k,xs,ys,zs); fflush(stderr);
            texture[itext].x = atof(xs);
            texture[itext].y = atof(ys);
            itext++;
         } else {
            // if its not identifiable, skip it, do not scale, do not write
         }

      } else if (onechar[0] == 'f') {
         // read a triangle line
         new_tri = alloc_new_tri();
         new_tri->index = num_tri;
         //fprintf(stdout,"t");
         //fprintf(stderr,"\nline %d, tri %d\n",nlines,new_tri->index);
         volatile int node_index[3] = {0,0,0};
         volatile int uv_index[3] = {0,0,0};
         volatile int norm_index[3] = {0,0,0};

         // read each set of values, like "39//123" or "5" or "1192/2524"
         for (int nn=0; nn<3; nn++) {
            char newval[32];
            fscanf(fp,"%s",newval);
            //if (nlines==247973) fprintf(stderr,"\n  key is %s\n",newval);
            // can't use sscanf, n//n isn't matched
            //int i1,i2,i3;
            //sscanf(newval,"%d/%d/%d",&i1,&i2,&i3);
            //if (nlines==247973) fprintf(stderr,"  tokens are %d %d %d\n",i1,i2,i3);
            // can't use strtok either, because n//n returns val/val and not val/null/val
            //char* token = strtok(newval,"/");
            //if (nlines==247973) fprintf(stderr,"  token is %s\n",token);
            //token = strtok(NULL,"/");
            //if (nlines==247973) fprintf(stderr,"  token is %s\n",token);
            //token = strtok(NULL,"/");
            //if (nlines==247973) fprintf(stderr,"  token is %s\n",token);

            // look for a slash
            int i;
            for (i=0; i<strlen(newval); i++) {
               if (newval[i] == '/' || newval[i] == '\0') break;
            }
            //if (nlines==247973) fprintf(stderr,"    i is %d\n",i);
            char sub1[10];
            strncpy(sub1,newval,i);
            //if (nlines==247973) fprintf(stderr,"    sub is %s\n",sub1);
            sub1[i] = '\0';
            if (i>0) node_index[nn] = atoi(sub1);
            //if (nlines==247973) fprintf(stderr,"    node %d is %d\n",0,node_index[0]);
            //if (nlines==247973) fprintf(stderr,"    node %d is %d\n",nn,node_index[nn]);

            // should we continue?
            if (newval[i] == '/') {
              int j;
              for (j=i+1; j<strlen(newval); j++) {
                 if (newval[j] == '/' || newval[j] == '\0') break;
              }
              char sub2[10];
              strncpy(sub2,newval+i+1,j);
              sub2[j-i-1] = '\0';
              if (j-i > 1) uv_index[nn] = atoi(sub2);
              //if (nlines==247973) fprintf(stderr,"    node %d is %d\n",nn,node_index[nn]);
              //if (nlines==247973) fprintf(stderr,"    nodeA %d is %d\n",0,node_index[0]);

              // should we continue?
              if (newval[j] == '/') {
                 int k;
                 for (k=j+1; k<strlen(newval); k++) {
                    if (newval[k] == '/' || newval[k] == '\0') break;
                 }
                 //if (nlines==247973) fprintf(stderr,"    nodeB %d is %d\n",0,node_index[0]);
                 char sub3[10];
                 strncpy(sub3,newval+j+1,k);
                 //if (nlines==247973) fprintf(stderr,"    nodeC %d is %d\n",0,node_index[0]);
                 sub3[k-j-1] = '\0';
                 //if (nlines==247973) fprintf(stderr,"    nodeD %d is %d\n",0,node_index[0]);
                 //if (nlines==247973) fprintf(stderr,"    norm %d is %d\n",nn,norm_index[nn]);
                 if (k-j > 1) norm_index[nn] = atoi(sub3);
                 //if (nlines==247973) fprintf(stderr,"    norm %d is %d\n",nn,norm_index[nn]);
                 //if (nlines==247973) fprintf(stderr,"    nodeE %d is %d\n",0,node_index[0]);
                 //if (nlines==247973) fprintf(stderr,"    node %d is %d\n",nn,node_index[nn]);
              }
            }
         }

         for (int j=0; j<3; j++) {
            // flip the normals here, by reversing the node order
            int targetj = j;
            if (invert) targetj = 2-j;

            if (node_index[j] < 1) {
               fprintf(stderr,"ERROR (read_obj): node index in file is <1 on line %d\nQuitting.\n",nlines);
               fprintf(stderr,"  nodes %d %d %d\n",node_index[0],node_index[1],node_index[2]);
               fprintf(stderr,"  texts %d %d %d\n",uv_index[0],uv_index[1],uv_index[2]);
               fprintf(stderr,"  norms %d %d %d\n",norm_index[0],norm_index[1],norm_index[2]);
               exit(1);
            }
            new_node = add_to_nodes_list(new_tri,&num_nodes,j,&loc[node_index[j]-1],&nodebin);
            new_tri->node[targetj] = new_node;
         }

         // set the node's normals here, if they were read in
         for (int j=0; j<3; j++) {
            int targetj = j;
            if (invert) targetj = 2-j;

            if (norm_index[j] < 1) {
               new_tri->norm[targetj] = NULL;
            } else {
               new_norm = add_to_norms_list(&num_norms,&normal[norm_index[j]-1],&normbin);
               new_tri->norm[targetj] = new_norm;
            }
         }

         // set the node's texture coord here, if they were read in
         for (int j=0; j<3; j++) {
            int targetj = j;
            if (invert) targetj = 2-j;

            if (uv_index[j] < 1) {
               new_tri->texture[targetj] = NULL;
            } else {
               new_text = add_to_textures_list(&num_texts,&texture[uv_index[j]-1],&textbin);
               new_tri->texture[targetj] = new_text;
            }
         }

         // set the adjacent triangle and midpoint pointers to NULL
         for (int j=0; j<3; j++) {
            new_tri->adjacent[j] = NULL;
            new_tri->midpoint[j] = NULL;
         }
         num_tri++;

         // add it on as the new head of the list
         if (tri_head) {
            new_tri->next_tri = tri_head;
            tri_head = new_tri;
         } else {
            tri_head = new_tri;
            tri_head->next_tri = NULL;
         }
         // fprintf(stderr,"  tri head is %d \n",tri_head->index);

         //fscanf(fp,"%[\n]",twochar);	/* read newline */
         if (num_tri/DOTPER == (num_tri+DPMO)/DOTPER)
            fprintf(stderr,".");
      }

      // skip the rest of the line, do not scale, do not write
      fscanf(fp,"%[^\n]",sbuf);	// read line beyond first char
      fscanf(fp,"%[\n]",twochar);	// read newline
   }
   fclose(fp);

   fprintf(stderr,"\n  %d tris\n",num_tri);
   fprintf(stderr,"  %d nodes\n",num_nodes);
   if (num_texts > 0) fprintf(stderr,"  %d texture coords\n",num_texts);
   if (num_norms > 0) fprintf(stderr,"  %d normals\n",num_norms);

   return(tri_head);
}
//#pragma GCC pop_options


/*
 * Read in a RAW file
 *
 * This subroutine assumes that the nodes of all triangles are 
 * ordered in a counter-clockwise pattern when viewed from the
 * top (outside).
 */
tri_pointer read_raw (char filename[80],tri_pointer tri_head) {

   int i,j,icnt;
   int num_tri = 0;
   int num_nodes = 0;
   int num_norms = 0;
   char twochar[2];
   char sbuf[512];
   char d[18][32];
   VEC location,normal,nmin,nmax;
   double temp[18];
   node_ptr new_node;
   norm_ptr new_norm;
   //text_ptr new_text;
   tri_pointer new_tri = NULL;
   BIN nodebin;
   NBIN normbin;
   FILE *fp;

   // open the .raw file for reading
   fp = fopen(filename,"r");
   if (fp==NULL) {
      fprintf(stderr,"Could not open input file %s\n",filename);
      fflush(stderr);
      exit(0);
   }
   fprintf(stderr,"Prescanning %s",filename);
   fflush(stderr);

   // important, icnt is the node counter variable!
   icnt = 0;
   nmax.x = -9.e+9;
   nmax.y = -9.e+9;
   nmax.z = -9.e+9;
   nmin.x = 9.e+9;
   nmin.y = 9.e+9;
   nmin.z = 9.e+9;
 
   // read the input file and update statistics

   // read a line from the input file
   while (fscanf(fp,"%[^\n]",sbuf) != EOF) {

      sscanf(sbuf,"%s %s %s %s %s %s %s %s %s",d[0],d[1],d[2],d[3],d[4],d[5],d[6],d[7],d[8]);

      if (!(int)isdigit(d[0][0]) && d[0][0] != '+' && d[0][0] != '-') {
         // read a comment line, or some other descriptor
         fscanf(fp,"%[\n]",twochar);	// read up to newline
         // fprintf(stdout,"%s\n",sbuf);	// write comment

      } else {
         // read a triangle line
         for (i=0; i<9; i++) temp[i] = atof(d[i]);

         update_minmax (&temp[0], &nmin, &nmax);
         update_minmax (&temp[3], &nmin, &nmax);
         update_minmax (&temp[6], &nmin, &nmax);
         if (icnt/DOTPER == (icnt+DPMO)/DOTPER) fprintf(stderr,".");
         icnt++;

         fscanf(fp,"%[^\n]",sbuf);	// read up to newline
         fscanf(fp,"%[\n]",twochar);	// read newline
      }
   }
   fprintf(stderr,"\n");
   fclose(fp);

   //fprintf(stderr,"found %d tris\n",icnt);

   // Initialize bin structure
   // fprintf(stderr,"%g %g %g   %g %g %g\n",nmin.x,nmin.y,nmin.z,nmax.x,nmax.y,nmax.z);
   //nodebin.dx = (nmax.x-nmin.x)/(BIN_COUNT-1);
   //nodebin.start = nmin.x - 0.5*nodebin.dx;
   //for (i=0;i<BIN_COUNT;i++) nodebin.b[i] = NULL;
   (void) prepare_node_bin (&nodebin,nmin,nmax);
   (void) prepare_norm_bin (&normbin);
 
 
   // open the file for reading
   fp = fopen(filename,"r");
   if (fp==NULL) {
      fprintf(stderr,"Could not open input file %s\n",filename);
      fflush(stderr);
      exit(0);
   }
   fprintf(stderr,"Opening file %s...",filename);
   fflush(stderr);

   // read a line from the input file
   while (fscanf(fp,"%[^\n]",sbuf) != EOF) {

      // read a triangle line
      sscanf(sbuf,"%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s",d[0],d[1],d[2],d[3],d[4],d[5],d[6],d[7],d[8],d[9],d[10],d[11],d[12],d[13],d[14],d[15],d[16],d[17]);

      if (!(int)isdigit(d[0][0]) && d[0][0] != '+' && d[0][0] != '-') {
         // read a comment line, or some other descriptor
         fscanf(fp,"%[\n]",twochar);	// read up to newline
         // fprintf(stdout,"%s\n",sbuf);	// write comment

      } else {
         for (i=0; i<18; i++) temp[i] = atof(d[i]);

         // and make the tri
         new_tri = alloc_new_tri();
         new_tri->index = num_tri;
         num_tri++;
         //fprintf(stderr,"tri %d\n",num_tri); fflush(stderr);

         // if we have 9 or more numbers, then we have a triangle
         // fprintf(stderr,"loc %g %g %g\n",location.x,location.y,location.z); fflush(stderr);
         //new_node = add_to_nodes_list(new_tri,&num_nodes,j,&location,&nodebin);
         location.x = temp[0];
         location.y = temp[1];
         location.z = temp[2];
         new_node = add_to_nodes_list(new_tri,&num_nodes,0,&location,&nodebin);
         new_tri->node[0] = new_node;
         //fprintf(stderr,"%d  %g %g %g\n",num_tri,location.x,location.y,location.z); fflush(stderr);

         location.x = temp[3];
         location.y = temp[4];
         location.z = temp[5];
         new_node = add_to_nodes_list(new_tri,&num_nodes,1,&location,&nodebin);
         new_tri->node[1] = new_node;
         //fprintf(stderr,"%d  %g %g %g\n",num_tri,location.x,location.y,location.z); fflush(stderr);

         location.x = temp[6];
         location.y = temp[7];
         location.z = temp[8];
         new_node = add_to_nodes_list(new_tri,&num_nodes,2,&location,&nodebin);
         new_tri->node[2] = new_node;
         //fprintf(stderr,"%d  %g %g %g\n",num_tri,location.x,location.y,location.z); fflush(stderr);

         /* set the node's normals here, if they were read in */

         // currently not supported (j<=9 in here)
         if ((int)isdigit(d[9][0]) || d[9][0] == '+' || d[9][0] == '-') {
            normal.x = temp[9];
            normal.y = temp[10];
            normal.z = temp[11];
            new_norm = add_to_norms_list(&num_norms,&normal,&normbin);
            new_tri->norm[0] = new_norm;

            normal.x = temp[12];
            normal.y = temp[13];
            normal.z = temp[14];
            new_norm = add_to_norms_list(&num_norms,&normal,&normbin);
            new_tri->norm[1] = new_norm;

            normal.x = temp[15];
            normal.y = temp[16];
            normal.z = temp[17];
            new_norm = add_to_norms_list(&num_norms,&normal,&normbin);
            new_tri->norm[2] = new_norm;

         } else {
            new_tri->norm[0] = NULL;
            new_tri->norm[1] = NULL;
            new_tri->norm[2] = NULL;
         }

         //if (fabs(new_tri->norm[2].x) < 1.e-5 &&
         //    fabs(new_tri->norm[2].y) < 1.e-5 &&
         //    fabs(new_tri->norm[2].z) < 1.e-5 ) {
         //   fprintf(stderr,"Found zero norm! j=%d\n",j);
         //   for (i=0; i<18; i++) fprintf(stderr,"  %d %g\n",i,temp[i]);
         //   exit(0);
         //}

         /* set the adjacent triangle and midpoint pointers to NULL */
         for (j=0; j<3; j++) {
            new_tri->adjacent[j] = NULL;
            new_tri->midpoint[j] = NULL;
         }

         /* add it on as the new head of the list */
         if (tri_head) {
            new_tri->next_tri = tri_head;
            tri_head = new_tri;
         } else {
            tri_head = new_tri;
            tri_head->next_tri = NULL;
         }

         fscanf(fp,"%[^\n]",sbuf);      // read line beyond first char
         fscanf(fp,"%[\n]",twochar);	// read newline
         if (num_tri/DOTPER == (num_tri+DPMO)/DOTPER) fprintf(stderr,".");

      }
   }
   fclose(fp);
   fprintf(stderr,"%d tris\n",num_tri);

   return(tri_head);
}


/*
 * Read in a TIN file
 *
 * This subroutine assumes that the nodes of all triangles are 
 * ordered in a counter-clockwise pattern when viewed from the
 * top (outside).
 */
tri_pointer read_tin (char filename[80],tri_pointer tri_head) {

   int i,j,icnt,err;
   int num_tri = 0;
   int num_nodes = 0;
   int num_norms = 0;
   char onechar;
   char twochar[2];
   char sbuf[256];
   char xs[20];
   VEC location,normal,nmin,nmax;
   double temp[9];
   node_ptr new_node;
   norm_ptr new_norm;
   //tri_pointer tri_head = NULL;
   tri_pointer new_tri = NULL;
   BIN nodebin;
   NBIN normbin;
   //long int fpos;
   FILE *fp;

   // open the file for reading
   fp = fopen(filename,"r");
   if (fp==NULL) {
      fprintf(stderr,"Could not open input file %s\n",filename);
      fflush(stderr);
      exit(0);
   }
   fprintf(stderr,"Prescanning %s",filename);
   fflush(stderr);

   // important, icnt is the node counter variable!
   icnt = 0;
   nmax.x = -9.e+9;
   nmax.y = -9.e+9;
   nmax.z = -9.e+9;
   nmin.x = 9.e+9;
   nmin.y = 9.e+9;
   nmin.z = 9.e+9;
 
   // read the input file and update statistics
   while (fread(&onechar,sizeof(char),1,fp) == 1) {
 
      if (onechar == '#') {
         /* read a comment line */
         fscanf(fp,"%[^\n]",sbuf);      /* read comment beyond '#' */
         fscanf(fp,"%[\n]",twochar);    /* read newline */
         // fprintf(stdout,"#%s\n",sbuf);       /* write comment */
                                                                                                 
      // if this is a data line, read it
      } else if (onechar == 't') {
         /* read a triangle line */

         // read three nodes
         for (j=0; j<9; j++) {
            //fpos = ftell (fp);
            err = fscanf(fp,"%s",xs);
            // check this value for non-math character
            if (err < 1 || xs[0] < 43 || xs[0] == 47 || xs[0] > 57) {
               // if it is, fill the rest of the data with zeros
               for (i=j; i<9; i++) temp[i] = 0.;
               // do not rewind the file
               //fseek (fp, fpos, SEEK_SET);
               // and break out of this loop
               break;
            } else {
               // if it's a math character, save it
               temp[j] = atof(xs);
               //fprintf(stderr," %g",temp[j]);
            }
         }
         //fprintf(stderr,"\nj %d\n",j);

         // if we have a complete triangle
         if (j>=9) {
            // fprintf(stdout,"loc %g %g %g\n",location.x,location.y,location.z); exit(0);
            update_minmax (&temp[0], &nmin, &nmax);
            update_minmax (&temp[3], &nmin, &nmax);
            update_minmax (&temp[6], &nmin, &nmax);
            if (icnt/DOTPER == (icnt+DPMO)/DOTPER) fprintf(stderr,".");
            icnt++;
         }

         fscanf(fp,"%[^\n]",sbuf);      // read line beyond first char
         fscanf(fp,"%[\n]",twochar);	/* read newline */
      } else {
         fscanf(fp,"%[^\n]",sbuf);      // read line beyond first char
         fscanf(fp,"%[\n]",twochar);	/* read newline */
      }
   }
   fprintf(stderr,"\n");
   fclose(fp);

   //fprintf(stderr,"found %d tris\n",icnt);

   // Initialize bin structure
   // fprintf(stderr,"%g %g %g   %g %g %g\n",nmin.x,nmin.y,nmin.z,nmax.x,nmax.y,nmax.z);
   //nodebin.dx = (nmax.x-nmin.x)/(BIN_COUNT-1);
   //nodebin.start = nmin.x - 0.5*nodebin.dx;
   //for (i=0;i<BIN_COUNT;i++) nodebin.b[i] = NULL;
   (void) prepare_node_bin (&nodebin,nmin,nmax);
   (void) prepare_norm_bin (&normbin);
 
 
   // open the file for reading
   fp = fopen(filename,"r");
   if (fp==NULL) {
      fprintf(stderr,"Could not open input file %s\n",filename);
      fflush(stderr);
      exit(0);
   }
   fprintf(stderr,"Opening file %s...",filename);
   fflush(stderr);

   // read the input file and update statistics
   while (fread(&onechar,sizeof(char),1,fp) == 1) {

      if (onechar == '#') {
         /* read a comment line */
         fscanf(fp,"%[^\n]",sbuf);	/* read comment beyond '#' */
         fscanf(fp,"%[\n]",twochar);	/* read newline */
         // fprintf(stdout,"#%s\n",sbuf);	// write comment

      // if this is a data line, read it
      } else if (onechar == 'n' || onechar == 't') {
         /* read a triangle line */

         // read three nodes or normals or both
         for (j=0; j<9; j++) {
            //fpos = ftell (fp);
            err = fscanf(fp,"%s",xs);
            // check this value for non-math character
            if (err < 1 || xs[0] < 43 || xs[0] == 47 || xs[0] > 57) {
               // if it is, fill the rest of the data with zeros
               for (i=j; i<9; i++) temp[i] = 0.;
               // do not rewind the file
               //fseek (fp, fpos, SEEK_SET);
               // and break out of this loop
               break;
            } else {
               // if it's a math character, save it
               temp[j] = atof(xs);
            }
         }

         // if there weren't even 9 numbers, then we have no triangle or
         // normal on this line
         if (j<9) continue;

         // set the normals
         if (onechar == 'n') {

            // if there is an active tri, save the normals
            if (new_tri) {
               normal.x = temp[0];
               normal.y = temp[1];
               normal.z = temp[2];
               new_norm = add_to_norms_list(&num_norms,&normal,&normbin);
               new_tri->norm[0] = new_norm;

               normal.x = temp[3];
               normal.y = temp[4];
               normal.z = temp[5];
               new_norm = add_to_norms_list(&num_norms,&normal,&normbin);
               new_tri->norm[1] = new_norm;

               normal.x = temp[6];
               normal.y = temp[7];
               normal.z = temp[8];
               new_norm = add_to_norms_list(&num_norms,&normal,&normbin);
               new_tri->norm[2] = new_norm;
            }

         // set the triangle
         } else {

         // and make the tri
         new_tri = alloc_new_tri();
         new_tri->index = num_tri;
         num_tri++;
         //fprintf(stderr,"tri %d\n",num_tri); fflush(stderr);

         // if we have 9 or more numbers, then we have a triangle
         // fprintf(stderr,"loc %g %g %g\n",location.x,location.y,location.z); fflush(stderr);
         //new_node = add_to_nodes_list(new_tri,&num_nodes,j,&location,&nodebin);
         location.x = temp[0];
         location.y = temp[1];
         location.z = temp[2];
         new_node = add_to_nodes_list(new_tri,&num_nodes,0,&location,&nodebin);
         new_tri->node[0] = new_node;
         //fprintf(stderr,"%d  %g %g %g\n",num_tri,location.x,location.y,location.z); fflush(stderr);

         location.x = temp[3];
         location.y = temp[4];
         location.z = temp[5];
         new_node = add_to_nodes_list(new_tri,&num_nodes,1,&location,&nodebin);
         new_tri->node[1] = new_node;
         //fprintf(stderr,"%d  %g %g %g\n",num_tri,location.x,location.y,location.z); fflush(stderr);

         location.x = temp[6];
         location.y = temp[7];
         location.z = temp[8];
         new_node = add_to_nodes_list(new_tri,&num_nodes,2,&location,&nodebin);
         new_tri->node[2] = new_node;
         //fprintf(stderr,"%d  %g %g %g\n",num_tri,location.x,location.y,location.z); fflush(stderr);

         /* set the adjacent triangle and midpoint pointers to NULL */
         for (j=0; j<3; j++) {
            new_tri->adjacent[j] = NULL;
            new_tri->midpoint[j] = NULL;
         }

         /* add it on as the new head of the list */
         if (tri_head) {
            new_tri->next_tri = tri_head;
            tri_head = new_tri;
         } else {
            tri_head = new_tri;
            tri_head->next_tri = NULL;
         }

         } // end setting the triangle

         fscanf(fp,"%[\n]",twochar);	/* read newline */
         if (num_tri/DOTPER == (num_tri+DPMO)/DOTPER) fprintf(stderr,".");

      } else {
         fscanf(fp,"%[^\n]",sbuf);      // read line beyond first char
         fscanf(fp,"%[\n]",twochar);	/* read newline */
      }
   }
   fclose(fp);
   fprintf(stderr,"%d tris\n",num_tri);

   return(tri_head);
}


/*
 * Read in a GMSH Mesh (.msh) file
 */
tri_pointer read_msh (char filename[80],tri_pointer tri_head) {

   int i,j;
   int num_tri = 0;
   int num_nodes = 0;
   int nodespace,num_elems;
   int index,node_index[3];
   //char onechar,anotherchar;
   char twochar[2];
   char sbuf[128];
   //char is[20],xs[30],ys[30],zs[30];
   char t[10][24];	// string tokens
   node_ptr *node_arry;		// vertex indexes
   node_ptr new_node = NULL;
   //tri_pointer tri_head = NULL;
   tri_pointer new_tri = NULL;
   FILE *fp;


   // open the .msh file for reading
   fp = fopen(filename,"r");
   if (fp==NULL) {
      fprintf(stderr,"Could not open input file %s\n",filename);
      fflush(stderr);
      exit(0);
   }
   fprintf(stderr,"Opening file %s...",filename);
   fflush(stderr);

   // first line is "$NOD"
   fscanf(fp,"%[^\n]",sbuf);	// read comment beyond '#'
   fscanf(fp,"%[\n]",twochar);	// read newline

   // second line is the number of nodes
   fscanf(fp,"%s",t[0]);
   num_nodes = atoi(t[0]);
   fscanf(fp,"%[^\n]",sbuf);	// read comment beyond '#'
   fscanf(fp,"%[\n]",twochar);	// read newline

   // now, we know how many nodes there will be: size and malloc the array
   nodespace = num_nodes+1000;
   node_arry = (node_ptr*)malloc(nodespace*sizeof(node_ptr));

   // and read all of the nodes
   for (i=0; i<num_nodes; i++) {
      fscanf(fp,"%s %s %s %s",t[0],t[1],t[2],t[3]);
      index = atoi(t[0]);

      // actually create the node instead of calling add_to_nodes_list
      new_node = (NODE *)malloc(sizeof(NODE));
      new_node->index = i;
      new_node->loc.x = atof(t[1]);
      new_node->loc.y = atof(t[2]);
      new_node->loc.z = atof(t[3]);
#ifdef CONN
      new_node->num_conn = 0;
      new_node->max_conn = 0;
#endif
      new_node->next_node = node_head;
      node_head = new_node;

      // add it to our array for easy indexing
      node_arry[index] = new_node;

      fscanf(fp,"%[^\n]",sbuf);	// read comment beyond '#'
      fscanf(fp,"%[\n]",twochar);	// read newline
   }

   // read lines with "$ENDNOD" and "$ELM"
   fscanf(fp,"%[^\n]",sbuf);	// read comment beyond '#'
   fscanf(fp,"%[\n]",twochar);	// read newline
   fscanf(fp,"%[^\n]",sbuf);	// read comment beyond '#'
   fscanf(fp,"%[\n]",twochar);	// read newline

   // then, read the number of elements
   fscanf(fp,"%s",t[0]);
   num_elems = atoi(t[0]);
   fscanf(fp,"%[^\n]",sbuf);	// read comment beyond '#'
   fscanf(fp,"%[\n]",twochar);	// read newline


   // and read all of the elements
   for (i=0; i<num_elems; i++) {
      fscanf(fp,"%s %s %s %s %s",t[0],t[1],t[2],t[3],t[4]);

      // if 2nd number is 1, it's a line; is 2, it's a tri
      if (atoi(t[1]) == 2) {

         // read a triangle line
         new_tri = alloc_new_tri();
         new_tri->index = num_tri;
         fscanf(fp,"%s %s %s",t[0],t[1],t[2]);
         node_index[0] = atoi(t[0]);
         node_index[1] = atoi(t[1]);
         node_index[2] = atoi(t[2]);
         for (j=0; j<3; j++) {
            if (node_index[j] < 1) {
               fprintf(stderr,"ERROR (read_msh): node index in file is <1\nQuitting.");
               exit(1);
            }
            new_tri->node[j] = node_arry[node_index[j]];
#ifdef CONN
            add_conn_tri(new_tri->node[j], new_tri, j);
            //new_tri->node[j]->conn_tri[new_tri->node[j]->num_conn] = new_tri;
            //new_tri->node[j]->conn_tri_node[new_tri->node[j]->num_conn] = j;
            //new_tri->node[j]->num_conn++;
#endif
         }

         // set the adjacent triangle and midpoint pointers to NULL
         for (j=0; j<3; j++) {
            new_tri->adjacent[j] = NULL;
            new_tri->midpoint[j] = NULL;
         }
         num_tri++;

         // add it on as the new head of the list
         if (tri_head) {
            new_tri->next_tri = tri_head;
            tri_head = new_tri;
         } else {
            tri_head = new_tri;
            tri_head->next_tri = NULL;
         }
         // fprintf(stderr,"  tri head is %d \n",tri_head->index);

         fscanf(fp,"%[\n]",twochar);	/* read newline */
         if (num_tri/DOTPER == (num_tri+DPMO)/DOTPER)
            fprintf(stderr,".");

      } else {

         //read to the end of the line and discard the info
         fscanf(fp,"%[^\n]",sbuf);	// read comment beyond '#'
         fscanf(fp,"%[\n]",twochar);	// read newline

      }
   }

   fclose(fp);
   fprintf(stderr,"%d tris\n",num_tri);

   return(tri_head);
}


/*
 * Determine the appropriate output file format from
 * the command-line and write the triangles to stdout
 */
int write_output(tri_pointer head, char format[4], int keep_norms, int argc, char **argv) {

   int num_wrote;

   if (strncmp(format, "raw", 3) == 0)
      num_wrote = write_raw(head, keep_norms);
   else if (strncmp(format, "pov", 1) == 0)
      num_wrote = write_pov(head, keep_norms);
   else if (strncmp(format, "rad", 3) == 0)
      num_wrote = write_rad(head, keep_norms);
   else if (strncmp(format, "tin", 1) == 0)
      num_wrote = write_tin(head, keep_norms);
   else if (strncmp(format, "obj", 1) == 0)
      num_wrote = write_obj(head, keep_norms, argc, argv);
   else if (strncmp(format, "rib", 2) == 0)
      num_wrote = write_rib(head, keep_norms);
   else if (strncmp(format, "seg", 2) == 0)
      num_wrote = write_seg(head, argc, argv);
   else {
      fprintf(stderr,"No output file format or unsupported file format given, using raw\n");
      num_wrote = write_raw(head, keep_norms);
   }

   return num_wrote;
}


/*
 * Write out a RAW file
 */
int write_raw (tri_pointer head, int keep_norms) {

   int i;
   int num = 0;
   tri_pointer curr_tri = head;

   fprintf(stderr,"Writing triangles to stdout");
   fflush(stderr);

   // run thru the triangle list and write them out
   while (curr_tri) {

      // write the node locations
      for (i=0; i<3; i++) {
         fprintf(stdout,"%lf %lf %lf ",curr_tri->node[i]->loc.x,curr_tri->node[i]->loc.y,curr_tri->node[i]->loc.z);
      }
      // and the node normals
      for (i=0; i<3; i++) {
         if (keep_norms && curr_tri->norm[i]) {
            fprintf(stdout,"%lf %lf %lf ",curr_tri->norm[i]->norm.x,curr_tri->norm[i]->norm.y,curr_tri->norm[i]->norm.z);
         }
      }
      fprintf(stdout,"\b\n");

      num++;
      if (num/DOTPER == (num+DPMO)/DOTPER)
         fprintf(stderr,".");
      curr_tri = curr_tri->next_tri;
   }
   fprintf(stderr,"\n");

   fprintf(stderr,"Wrote RAW ASCII information, %d triangles\n",num);
   return num;
}


/*
 * Write out a TIN file
 */
int write_tin (tri_pointer head, int keep_norms) {

   int i;
   int num = 0;
   tri_pointer curr_tri = head;

   fprintf(stderr,"Writing triangles to stdout");
   fflush(stderr);

   /* run thru the triangle list and write them out */
   while (curr_tri) {

      /* write the node normals */
      if (keep_norms && curr_tri->norm[0] && curr_tri->norm[1] && curr_tri->norm[2]) {
         fprintf(stdout,"n");
         for (i=0; i<3; i++) {
            fprintf(stdout," %lf %lf %lf",curr_tri->norm[i]->norm.x,curr_tri->norm[i]->norm.y,curr_tri->norm[i]->norm.z);
         }
         fprintf(stdout,"\n");
      }

      /* write the node locations */
      fprintf(stdout,"t");
      for (i=0; i<3; i++) {
         fprintf(stdout," %lf %lf %lf",curr_tri->node[i]->loc.x,curr_tri->node[i]->loc.y,curr_tri->node[i]->loc.z);
      }
      fprintf(stdout,"\n");

      num++;
      if (num/DOTPER == (num+DPMO)/DOTPER)
         fprintf(stderr,".");
      curr_tri = curr_tri->next_tri;
   }
   fprintf(stderr,"\n");

   fprintf(stderr,"Wrote TIN ASCII information, %d triangles\n",num);
   return num;
}


/*
 * Write out a Wavefront .obj file
 *
 * compact notation is where each unique node is listed once, if this
 * is set to FALSE, preceding each triangle will be the three nodes' locations
 */
int write_obj(tri_pointer head, int keep_norms, int argc, char **argv) {

   int tri_ct = 0;
   int node_ct = 0;
   int norm_ct = 0;
   int text_ct = 0;
   // int normal_ct = 0;
   // node_ptr curr_node = node_head;
   tri_pointer curr_tri = head;

   fprintf(stderr,"Writing triangles in Wavefront .obj format to stdout\n");
   fflush(stderr);
   fprintf(stdout,"# Triangle mesh written from rocktools\n#");
   for (int i=0; i<argc; i++) {
      fprintf(stdout," %s",argv[i]);
   }
   fprintf(stdout,"\no tri_mesh\n");

   // set all node and norm indexes to -1
   curr_tri = head;
   while (curr_tri) {
      // fprintf(stderr,"tri %d\n",curr_tri->index); fflush(stderr);
      // for (i=0; i<3; i++) fprintf(stderr,"  node %d\n",curr_tri->node[i]->index); fflush(stderr);
      for (int i=0; i<3; i++) curr_tri->node[i]->index = -1;
      for (int i=0; i<3; i++) if (curr_tri->norm[i] != NULL) curr_tri->norm[i]->index = -1;
      for (int i=0; i<3; i++) if (curr_tri->texture[i] != NULL) curr_tri->texture[i]->index = -1;
      curr_tri = curr_tri->next_tri;
   }

   // march through all triangles, writing nodes as they appear, 
   //    and setting indexes as they appear
   curr_tri = head;
   while (curr_tri) {
      for (int i=0; i<3; i++) {
         if (curr_tri->node[i]->index == -1) {
            // this node has not appeared, write it!
            if (isnan(curr_tri->node[i]->loc.x) ||
                isnan(curr_tri->node[i]->loc.y) ||
                isnan(curr_tri->node[i]->loc.z)) {
              curr_tri->node[i]->loc.x = curr_tri->node[i%3]->loc.x;
              curr_tri->node[i]->loc.y = curr_tri->node[i%3]->loc.y;
              curr_tri->node[i]->loc.z = curr_tri->node[i%3]->loc.z;
            }
            fprintf(stdout,"v %12.7e %12.7e %12.7e\n",curr_tri->node[i]->loc.x,curr_tri->node[i]->loc.y,curr_tri->node[i]->loc.z);
            // and set the index
            curr_tri->node[i]->index = ++node_ct;
         }
      }

      // write any valid normals
      if (keep_norms) {
         // write the unique norms
         for (int i=0; i<3; i++) {
            if (curr_tri->norm[i] != NULL) {
               if (curr_tri->norm[i]->index == -1 &&
                  !isnan(curr_tri->norm[i]->norm.x) &&
                  !isnan(curr_tri->norm[i]->norm.y) &&
                  !isnan(curr_tri->norm[i]->norm.z)) {
                  // this node has not appeared and is not nan, write it!
                  fprintf(stdout,"vn %12.7e %12.7e %12.7e\n",curr_tri->norm[i]->norm.x,curr_tri->norm[i]->norm.y,curr_tri->norm[i]->norm.z);
                  // and set the index
                  curr_tri->norm[i]->index = ++norm_ct;
               }
            }
         }
      } else {
            for (int i=0; i<3; i++) free(curr_tri->norm[i]);
      }

      // write any valid texture coordinates
      // write the unique norms
      for (int i=0; i<3; i++) {
         if (curr_tri->texture[i] != NULL) {
            if (curr_tri->texture[i]->index == -1 &&
               !isnan(curr_tri->texture[i]->uv.x) &&
               !isnan(curr_tri->texture[i]->uv.y)) {
               // this node has not appeared and is not nan, write it!
               fprintf(stdout,"vt %12.7e %12.7e\n",curr_tri->texture[i]->uv.x,curr_tri->texture[i]->uv.y);
               // and set the index
               curr_tri->texture[i]->index = ++text_ct;
            }
         }
      }

      // and write the tri itself
      fprintf(stdout,"f");
      for (int i=0; i<3; i++) {
         fprintf(stdout," %d",curr_tri->node[i]->index);
         if (curr_tri->texture[i] != NULL) {
            fprintf(stdout,"/%d",curr_tri->texture[i]->index);
            if (curr_tri->norm[i] != NULL) fprintf(stdout,"/%d",curr_tri->norm[i]->index);
         } else {
            if (curr_tri->norm[i] != NULL) fprintf(stdout,"//%d",curr_tri->norm[i]->index);
         }
      }
      fprintf(stdout,"\n");

      tri_ct++;
      curr_tri = curr_tri->next_tri;
   }

   fprintf(stderr,"Wrote Wavefront ASCII information, %d triangles\n",tri_ct);
   return tri_ct;
}


/*
 * Write out a Wavefront-like .seg file
 *
 * compact notation is where each unique node is listed once, if this
 * is set to FALSE, preceding each triangle will be the three nodes' locations
 */
int write_seg(tri_pointer head, int argc, char **argv) {

   int tri_ct = 0;
   int node_ct = 0;
   tri_pointer curr_tri = head;

   fprintf(stderr,"Writing segments in Wavefront-like .seg format to stdout\n");
   fflush(stderr);
   fprintf(stdout,"# Segmented network written from rockconvert\n#");
   for (int i=0; i<argc; i++) {
      fprintf(stdout," %s",argv[i]);
   }
   fprintf(stdout,"\n");
   //fprintf(stdout,"d 3\n");
   //fprintf(stdout,"gr 0.01\n");

   // set all node indexes to -1
   curr_tri = head;
   while (curr_tri) {
      for (int i=0; i<3; i++) curr_tri->node[i]->index = -1;
      curr_tri = curr_tri->next_tri;
   }

   // march through all triangles, writing nodes and segments as they appear, 
   //    and setting indexes as they appear
   curr_tri = head;
   while (curr_tri) {

      // write any new nodes
      for (int i=0; i<3; i++) {
         if (curr_tri->node[i]->index == -1) {
            // this node has not appeared, write it!
            if (isnan(curr_tri->node[i]->loc.x) ||
                isnan(curr_tri->node[i]->loc.y) ||
                isnan(curr_tri->node[i]->loc.z)) {
              curr_tri->node[i]->loc.x = curr_tri->node[i%3]->loc.x;
              curr_tri->node[i]->loc.y = curr_tri->node[i%3]->loc.y;
              curr_tri->node[i]->loc.z = curr_tri->node[i%3]->loc.z;
            }
            fprintf(stdout,"v %12.7e %12.7e %12.7e\n",curr_tri->node[i]->loc.x,curr_tri->node[i]->loc.y,curr_tri->node[i]->loc.z);
            // and set the index
            curr_tri->node[i]->index = ++node_ct;
         }
      }

      // and write the three segments (these may not be unique!)
      for (int i=0; i<3; i++) {
         int ip1 = (i+1)%3;
         fprintf(stdout,"s %d %d\n", curr_tri->node[i]->index, curr_tri->node[ip1]->index);
      }

      tri_ct++;
      curr_tri = curr_tri->next_tri;
   }

   fprintf(stderr,"Wrote .seg ASCII information, %d segments\n",3*tri_ct);
   return tri_ct;
}


/*
 * Write out POV-readable syntax
 */
int write_pov (tri_pointer head, int keep_norms) {

   int i;
   int num = 0;
   tri_pointer curr_tri = head;

   fprintf(stderr,"Writing triangles to stdout");
   fflush(stderr);
   fprintf(stdout,"mesh {\n");

   /* run thru the triangle list and write them out */
   while (curr_tri) {

      if (keep_norms && curr_tri->norm[0] && curr_tri->norm[1] && curr_tri->norm[2]) {

         /* write the triangle with surface normals defined */
         fprintf(stdout,"smooth_triangle {\n");
         for (i=0; i<2; i++) {
            fprintf(stdout,"   <%g, %g, %g>",curr_tri->node[i]->loc.x,curr_tri->node[i]->loc.y,curr_tri->node[i]->loc.z);
            fprintf(stdout,", <%g, %g, %g>",curr_tri->norm[i]->norm.x,curr_tri->norm[i]->norm.y,curr_tri->norm[i]->norm.z);
            if (i!=2) fprintf(stdout,",");
            fprintf(stdout,"\n");
         }
         fprintf(stdout,"}\n");

      } else {

         /* write the triangle */
         fprintf(stdout,"triangle {");
         for (i=0; i<3; i++) {
            fprintf(stdout," <%g, %g, %g>",curr_tri->node[i]->loc.x,curr_tri->node[i]->loc.y,curr_tri->node[i]->loc.z);
         }
         fprintf(stdout,"}\n");

      }

      num++;
      if (num/DOTPER == (num+DPMO)/DOTPER)
         fprintf(stderr,".");
      curr_tri = curr_tri->next_tri;
   }
   fprintf(stderr,"\n");

   fprintf(stdout,"}\n");

   fprintf(stderr,"Wrote POV ASCII information, %d triangles\n",num);
   return num;
}


/*
 * Write out Radiance-readable syntax
 */
int write_rad (tri_pointer head, int keep_norms) {

   int i;
   int num = 0;
#ifdef ERODE
   double rad_start;
   node_ptr curr_node;
#endif
   tri_pointer curr_tri = head;

   fprintf(stderr,"Writing triangles to stdout");
   fflush(stderr);

   /* run thru the triangle list and write them out */
   while (curr_tri) {

      /* write the triangle */
      fprintf(stdout,"default polygon p%d\n0 0 9",num);
      for (i=0; i<3; i++) {
         fprintf(stdout," %g %g %g",curr_tri->node[i]->loc.x,curr_tri->node[i]->loc.y,curr_tri->node[i]->loc.z);
      }
      fprintf(stdout,"\n");

      num++;
      if (num/DOTPER == (num+DPMO)/DOTPER) fprintf(stderr,".");
      curr_tri = curr_tri->next_tri;
   }
   fprintf(stderr,"\n");

#ifdef ERODE
   /* if flow data is there, write it out, too */
   if (1==1) {
      i = 0;

      curr_node = node_head;
      while (curr_node) {
         // was >0.00001 and 0.02*sqrt
         if (curr_node->flow_rate > 0.0001) {
            rad_start = 0.01*sqrt(curr_node->flow_rate);
            if (curr_node->downstream) {
               /* rad_end = 0.02*sqrt(curr_node->downstream->flow_rate); */
               fprintf(stdout,"water sphere w_%d\n0 0 4 %lf %lf %lf %lf\n",i++,
                  curr_node->loc.x,curr_node->loc.y,curr_node->loc.z,rad_start);
               fprintf(stdout,"water cylinder w_%d\n0 0 7 %lf %lf %lf %lf %lf %lf %lf\n",i++,
                  curr_node->loc.x,curr_node->loc.y,curr_node->loc.z,
                  curr_node->downstream->loc.x,curr_node->downstream->loc.y,curr_node->downstream->loc.z,
                  rad_start);
            } else {
               fprintf(stdout,"water sphere w_%d\n0 0 4 %lf %lf %lf %lf\n",i++,
                  curr_node->loc.x,curr_node->loc.y,curr_node->loc.z,rad_start);
            }
         }
         curr_node = curr_node->next_node;
      }
   }
#endif

   fprintf(stderr,"Wrote Radiance ASCII information, %d triangles\n",num);
   return num;
}


/*
 * Write out Renderman-compliant syntax (RIB format for BMRT)
 */
int write_rib(tri_pointer head, int keep_norms) {

   int i;
   int num = 0;
   tri_pointer curr_tri = head;

   fprintf(stderr,"Writing triangles to stdout");
   fflush(stderr);

   /* run thru the triangle list and write them out */
   while (curr_tri) {

      /* write the triangle */
      fprintf(stdout,"Polygon \"P\" [");
      /* vertexes are ordered clockwise, not CCW like all other formats */
      for (i=2; i>-1; i--) {
         fprintf(stdout," %g %g %g",curr_tri->node[i]->loc.x,curr_tri->node[i]->loc.y,curr_tri->node[i]->loc.z);
      }

      if (keep_norms && curr_tri->norm[0] && curr_tri->norm[1] && curr_tri->norm[2]) {
         fprintf(stdout,"] \"N\" [");
         for (i=2; i>-1; i--) {
            fprintf(stdout," %g %g %g",curr_tri->norm[i]->norm.x,curr_tri->norm[i]->norm.y,curr_tri->norm[i]->norm.z);
         }
      }

      fprintf(stdout," ]\n");

      /* normals are stored like so: 
      Polygon "P" [ 0 1 0  0 8 0  4 4 0 ] "N" [ 1 0 0  1 0 0  0 1 0 ]
      */

      num++;
      if (num/DOTPER == (num+DPMO)/DOTPER) fprintf(stderr,".");
      curr_tri = curr_tri->next_tri;
   }
   fprintf(stderr,"\n");

   fprintf(stderr,"Wrote RIB ASCII information, %d triangles\n",num);
   return num;
}


/*
 * find_mesh_stats
 *
 * parse through a raw, tin, rad, or obj file and find the min and max
 * bounds, and the triangle and node count (if possible)
 */
int find_mesh_stats(char *infile, VEC *bmin, VEC *bmax, int *num_tris, int *num_nodes) {

   int i;
   int input_format = 0;	// integer flag for input file type
   char onechar,anotherchar;
   char twochar[2];
   char sbuf[128];
   char xs[20],ys[20],zs[20];
   char extension[4];		/* filename extension if infile */
   tri_pointer the_tri;
   FILE *ifp;

   (*num_nodes) = 0;
   (*num_tris) = 0;

   bmin->x = 9.9e+9;
   bmax->x = -9.9e+9;
   bmin->y = 9.9e+9;
   bmax->y = -9.9e+9;
   bmin->z = 9.9e+9;
   bmax->z = -9.9e+9;


   /* Determine the input file format from the .XXX extension, and read it */
   strncpy(extension,infile+strlen(infile)-3,4);
   if (strncmp(extension, "raw", 3) == 0)
      input_format = 1;
   else if (strncmp(extension, "tin", 1) == 0)
      input_format = 2;
   else if (strncmp(extension, "rad", 1) == 0)
      input_format = 3;
   else if (strncmp(extension, "obj", 1) == 0)
      input_format = 4;
   else {
      fprintf(stderr,"Input filename extension is assumed to be (%s)\n",extension);
      fprintf(stderr,"This is either not a supported input file format for this program,\n");
      fprintf(stderr,"   or you need to use a proper filename extension.\n");
      fprintf(stderr,"Supported input file formats are: .raw, .tin, .rad, and .obj\n");
      exit(0);
   }


   // Set up memory space for the working triangle, and the three nodes
   the_tri = alloc_new_tri();
   for (i=0; i<3; i++) the_tri->node[i] = (NODE *)malloc(sizeof(NODE));
   for (i=0; i<3; i++) the_tri->norm[i] = (NORM *)malloc(sizeof(NORM));

   // Open the file for reading
   ifp = fopen(infile,"r");
   if (ifp==NULL) {
      fprintf(stderr,"Could not open input file %s\n",infile);
      exit(0);
   }
   fprintf(stderr,"Opening file %s",infile);
   fflush(stderr);

   // for rad, tin, raw, use this method: ------------------------------------

   if (input_format < 4) {

   // as long as there are triangles available, operate
   while (get_tri(ifp,input_format,the_tri) == 1) {

      (*num_tris)++;

      // find the min/max bounds
      for (i=0; i<3; i++) {
         if (the_tri->node[i]->loc.x < bmin->x) bmin->x = the_tri->node[i]->loc.x;
         if (the_tri->node[i]->loc.x > bmax->x) bmax->x = the_tri->node[i]->loc.x;
         if (the_tri->node[i]->loc.y < bmin->y) bmin->y = the_tri->node[i]->loc.y;
         if (the_tri->node[i]->loc.y > bmax->y) bmax->y = the_tri->node[i]->loc.y;
         if (the_tri->node[i]->loc.z < bmin->z) bmin->z = the_tri->node[i]->loc.z;
         if (the_tri->node[i]->loc.z > bmax->z) bmax->z = the_tri->node[i]->loc.z;
      }

      // print a dot every DOTPER triangles
      if ((*num_tris)/DOTPER == ((*num_tris)+DPMO)/DOTPER) fprintf(stderr,".");
   }
   fprintf(stderr,"\n");

   *num_nodes = (*num_tris)*3;

   // for obj, use this method: ----------------------------------------------

   } else {

   VEC test;
   // read the input file and update statistics
   while (fread(&onechar,sizeof(char),1,ifp) == 1) {

      // only read nodes, we're putting them into bins on this pass
      if (onechar == 'v') {

         fread(&anotherchar,sizeof(char),1,ifp);
         // fprintf(stdout,"  (%c)\n",anotherchar); fflush(stdout);

         if (isspace(anotherchar)) {
            // read a vertex location
            fscanf(ifp,"%s %s %s",xs,ys,zs);
            test.x = atof(xs);
            test.y = atof(ys);
            test.z = atof(zs);
            if (test.x > bmax->x) bmax->x = test.x;
            if (test.y > bmax->y) bmax->y = test.y;
            if (test.z > bmax->z) bmax->z = test.z;
            if (test.x < bmin->x) bmin->x = test.x;
            if (test.y < bmin->y) bmin->y = test.y;
            if (test.z < bmin->z) bmin->z = test.z;
            (*num_nodes)++;
         } else {
            // if its not identifiable, skip it, do not scale, do not write
         }
      } else if (onechar == 'f') {
         // read a triangle line
         (*num_tris)++;
         // print a dot every DOTPER triangles
         if ((*num_tris)/DOTPER == ((*num_tris)+DPMO)/DOTPER) fprintf(stderr,".");
      }
      // if its not identifiable, skip it, do not scale, do not write
      fscanf(ifp,"%[^\n]",sbuf);      // read line beyond first char
      fscanf(ifp,"%[\n]",twochar);    // read newline
   }
   fprintf(stderr,"\n");

   }

   // why does it die on closing?
   //fclose(ifp);

   return(0);
}


/*
 * get_tri will read the appropriate format input file 
 * for the neccesary data on the next triangle
 */
int get_tri(FILE *input, int input_format, tri_pointer current_tri) {

   /* int have_normals = 0; */
   char first_token[32];
   char twochar[2];
   char sbuf[512];
   char d[18][32];
   char k[3][32];

   // free any memory associated with the normals
   for (int i=0; i<3; i++) {
      if (current_tri->norm[i] != NULL) free(current_tri->norm[i]);
   }

   /* read a line from the input file */
   while (fscanf(input,"%[^\n]",sbuf) != EOF) {

      if (input_format == 1) {
         /* Need to read from a .raw file, might have normal information */

         sscanf(sbuf,"%s",first_token);
         if (!(int)isdigit(first_token[0]) && first_token[0] != '+' && first_token[0] != '-') {
            /* read a comment line, or some other descriptor */
            fscanf(input,"%[\n]",twochar);	// read up to newline
            // fprintf(stdout,"%s\n",sbuf);	// write comment

         } else {
            /* read a triangle line */
            sscanf(sbuf,"%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s",d[0],d[1],d[2],d[3],d[4],d[5],d[6],d[7],d[8],d[9],d[10],d[11],d[12],d[13],d[14],d[15],d[16],d[17]);
            fscanf(input,"%[\n]",twochar);	/* read newline */

            /* assign values to triangle */
            current_tri->node[0]->loc.x = atof(d[0]);
            current_tri->node[0]->loc.y = atof(d[1]);
            current_tri->node[0]->loc.z = atof(d[2]);
            current_tri->node[1]->loc.x = atof(d[3]);
            current_tri->node[1]->loc.y = atof(d[4]);
            current_tri->node[1]->loc.z = atof(d[5]);
            current_tri->node[2]->loc.x = atof(d[6]);
            current_tri->node[2]->loc.y = atof(d[7]);
            current_tri->node[2]->loc.z = atof(d[8]);

            if ((int)isdigit(d[9][0]) || d[9][0] == '+' || d[9][0] == '-') {
              current_tri->norm[0] = (norm_ptr)malloc(sizeof(NORM));
              current_tri->norm[1] = (norm_ptr)malloc(sizeof(NORM));
              current_tri->norm[2] = (norm_ptr)malloc(sizeof(NORM));
              current_tri->norm[0]->norm.x = atof(d[9]);
              current_tri->norm[0]->norm.y = atof(d[10]);
              current_tri->norm[0]->norm.z = atof(d[11]);
              current_tri->norm[1]->norm.x = atof(d[12]);
              current_tri->norm[1]->norm.y = atof(d[13]);
              current_tri->norm[1]->norm.z = atof(d[14]);
              current_tri->norm[2]->norm.x = atof(d[15]);
              current_tri->norm[2]->norm.y = atof(d[16]);
              current_tri->norm[2]->norm.z = atof(d[17]);
              //fprintf(stderr,"  found norm %g %g\n",current_tri->norm[0].x,current_tri->norm[0].y);
            }
            fscanf(input,"%[\n]",twochar);	/* read newline */

            /* we have the triangle, now jump out of the "while" */
            return(1);
         }

      } else if (input_format == 2) {
         /* Need to read a .tin file, may have normal information */

         sscanf(sbuf,"%s",first_token);
         if (first_token[0] == 't') {
            /* read a triangle line */
            sscanf(sbuf,"%s %s %s %s %s %s %s %s %s %s",k[0],d[0],d[1],d[2],d[3],d[4],d[5],d[6],d[7],d[8]);
            fscanf(input,"%[\n]",twochar);	/* read newline */
            /* fprintf(stderr,"z-values are %s %s %s\n",d[2],d[5],d[8]); */

            /* assign values to triangle */
            current_tri->node[0]->loc.x = atof(d[0]);
            current_tri->node[0]->loc.y = atof(d[1]);
            current_tri->node[0]->loc.z = atof(d[2]);
            current_tri->node[1]->loc.x = atof(d[3]);
            current_tri->node[1]->loc.y = atof(d[4]);
            current_tri->node[1]->loc.z = atof(d[5]);
            current_tri->node[2]->loc.x = atof(d[6]);
            current_tri->node[2]->loc.y = atof(d[7]);
            current_tri->node[2]->loc.z = atof(d[8]);

            /* we have the triangle, now jump out of the "while" */
            return(1);

         } else if (first_token[0] == 'n') {
            /* read the list of surface normal vectors */
            sscanf(sbuf,"%s %s %s %s %s %s %s %s %s %s",k[0],d[0],d[1],d[2],d[3],d[4],d[5],d[6],d[7],d[8]);
            fscanf(input,"%[\n]",twochar);	/* read newline */
            /* fprintf(stderr,"z-values are %s %s %s\n",d[2],d[5],d[8]); */

            /* assign values to triangle */
            current_tri->norm[0] = (norm_ptr)malloc(sizeof(NORM));
            current_tri->norm[1] = (norm_ptr)malloc(sizeof(NORM));
            current_tri->norm[2] = (norm_ptr)malloc(sizeof(NORM));
            current_tri->norm[0]->norm.x = atof(d[0]);
            current_tri->norm[0]->norm.y = atof(d[1]);
            current_tri->norm[0]->norm.z = atof(d[2]);
            current_tri->norm[1]->norm.x = atof(d[3]);
            current_tri->norm[1]->norm.y = atof(d[4]);
            current_tri->norm[1]->norm.z = atof(d[5]);
            current_tri->norm[2]->norm.x = atof(d[6]);
            current_tri->norm[2]->norm.y = atof(d[7]);
            current_tri->norm[2]->norm.z = atof(d[8]);

         } else {
            // read a comment line, or some other descriptor
            fscanf(input,"%[\n]",twochar);	// read up to newline
            // fprintf(stdout,"%s\n",sbuf);	// write comment
         }

      } else if (input_format == 3) {
         /* Need to read a Radiance file, one written by rocktools, mind you.
          * default polygon p0
          * 0 0 9 0.602705 0.113831 0.380334 0.599856 0.118667 0.38346 0.596791 0.115323 0.378713
          */

         sscanf(sbuf,"%s",first_token);
         if (!isdigit(first_token[0]) && first_token[0] != '+' && first_token[0] != '-') {
            // read a comment line, or some other descriptor
            fscanf(input,"%[\n]",twochar);    // read up to newline
            // fprintf(stdout,"%s\n",sbuf);        // write comment

         } else {
            /* read a triangle line */
            sscanf(sbuf,"%s %s %s %s %s %s %s %s %s %s %s %s",k[0],k[1],k[2],d[0],d[1],d[2],d[3],d[4],d[5],d[6],d[7],d[8]);
            fscanf(input,"%[\n]",twochar);    /* read newline */

            /* assign values to triangle */
            current_tri->node[0]->loc.x = atof(d[0]);
            current_tri->node[0]->loc.y = atof(d[1]);
            current_tri->node[0]->loc.z = atof(d[2]);
            current_tri->node[1]->loc.x = atof(d[3]);
            current_tri->node[1]->loc.y = atof(d[4]);
            current_tri->node[1]->loc.z = atof(d[5]);
            current_tri->node[2]->loc.x = atof(d[6]);
            current_tri->node[2]->loc.y = atof(d[7]);
            current_tri->node[2]->loc.z = atof(d[8]);

            /* we have the triangle, now jump out of the "while" */
            return(1);
         }
      }
   }

   /* if we got this far, there weren't any more triangles to read */
   return(0);
}


/*
 * write_tri will write the current triangle to the appropriate
 * output file format specified.
 */
int write_tri(FILE *output, int output_format, tri_pointer current_tri) {

   if (output_format == 1) {

      /* write the triangle in .raw format */
      fprintf(output,"%g %g %g %g %g %g %g %g %g",
         current_tri->node[0]->loc.x,current_tri->node[0]->loc.y,current_tri->node[0]->loc.z,
         current_tri->node[1]->loc.x,current_tri->node[1]->loc.y,current_tri->node[1]->loc.z,
         current_tri->node[2]->loc.x,current_tri->node[2]->loc.y,current_tri->node[2]->loc.z);
      if (current_tri->norm[0] && current_tri->norm[1] && current_tri->norm[2]) {
         fprintf(output," %g %g %g %g %g %g %g %g %g\n",
           current_tri->norm[0]->norm.x,current_tri->norm[0]->norm.y,current_tri->norm[0]->norm.z,
           current_tri->norm[1]->norm.x,current_tri->norm[1]->norm.y,current_tri->norm[1]->norm.z,
           current_tri->norm[2]->norm.x,current_tri->norm[2]->norm.y,current_tri->norm[2]->norm.z);
      } else {
         fprintf(output,"\n");
      }

   } else if (output_format == 2) {

      /* write the triangle in .tin format */
      fprintf(output,"t %g %g %g %g %g %g %g %g %g\n",
         current_tri->node[0]->loc.x,current_tri->node[0]->loc.y,current_tri->node[0]->loc.z,
         current_tri->node[1]->loc.x,current_tri->node[1]->loc.y,current_tri->node[1]->loc.z,
         current_tri->node[2]->loc.x,current_tri->node[2]->loc.y,current_tri->node[2]->loc.z);
      if (current_tri->norm[0] && current_tri->norm[1] && current_tri->norm[2]) {
         fprintf(output,"n %g %g %g %g %g %g %g %g %g\n",
           current_tri->norm[0]->norm.x,current_tri->norm[0]->norm.y,current_tri->norm[0]->norm.z,
           current_tri->norm[1]->norm.x,current_tri->norm[1]->norm.y,current_tri->norm[1]->norm.z,
           current_tri->norm[2]->norm.x,current_tri->norm[2]->norm.y,current_tri->norm[2]->norm.z);
      }


   } else if (output_format == 3) {

      /* write the triangle in .rad format */
      fprintf(output,"default polygon p%d\n0 0 9 ",0);
      fprintf(output,"%g %g %g %g %g %g %g %g %g\n",
         current_tri->node[0]->loc.x,current_tri->node[0]->loc.y,current_tri->node[0]->loc.z,
         current_tri->node[1]->loc.x,current_tri->node[1]->loc.y,current_tri->node[1]->loc.z,
         current_tri->node[2]->loc.x,current_tri->node[2]->loc.y,current_tri->node[2]->loc.z);

   }

   return(1);
}

