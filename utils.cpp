#include "utils.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*
 * generate a polygon with num_verts vertices.
 * all the vertices are at (0, 0)
 */
Poly* generatePoly(int num_verts) {
  Poly *poly;

  if(num_verts<=0){
    return NULL;
  }
  poly = (Poly*)calloc(1, sizeof(Poly));
  poly->num_verts = num_verts;
  poly->verts = (Vertex*)calloc(num_verts, sizeof(Vertex));
  return poly;
}

/*
 * deallocate memory used by poly. must not use poly after that.
 * if delete_poly==1 then free memory of poly
 */
void deletePoly(Poly *poly, int delete_poly) {
  int i;
  if(poly == NULL){
    return;
  }
  if(poly->verts != NULL){
    free(poly->verts);
  }
  if(poly->triangles != NULL){
    for(i=0; i<poly->num_triangles; i++){
      deletePoly(poly->triangles+i, 0);
    }
    free(poly->triangles);
  }
  if(poly->normals != NULL){
    free(poly->normals);
  }
  if(poly->dists != NULL){
    free(poly->dists);
  }
  if(delete_poly){
    free(poly);
  }
}

/*
 * print details of a Poly object.
 * if str!=NULL, then write those details in str instead
 * (assuming str is big enough).
 * for debugging only.
 */
void printPoly(Poly *poly, char *str){
  int delete_str = 0;
  int i;
  if(str==NULL){
    str = (char*)malloc(1024);
    delete_str = 1;
  }
  if(poly==NULL){
    sprintf(str, "invalid poly (NULL)");
  }
  else{
    Vertex *v;
    char substr[256];
    sprintf(str, "%d vertices: ", poly->num_verts);
    for(v=poly->verts; v<poly->verts+poly->num_verts; v++){
      sprintf(substr, "(%+01.3f, %+01.3f)", v->x, v->y);
      strcat(str, substr);
      if(v<(poly->verts+(poly->num_verts-1))){
        strcat(str, ", ");
      }
    }
    sprintf(substr, "; area: %.3f", poly->area);
    strcat(str, substr);
    if(poly->normals!=NULL && poly->dists!=NULL){
      Vector *n;
      strcat(str, "; sides: ");
      for(n=poly->normals, i=0; i<poly->num_verts; i++, n++){
        sprintf(substr, "{(%+01.3f, %+01.3f), %+01.3f}", n->x, n->y, poly->dists[i]);
        strcat(str, substr);
        if(i<(poly->num_verts-1)){
          strcat(str, ", ");
        }
      }
    }
    if(poly->num_triangles>0){
      strcat(str, "; triangulation: ");
      for(i=0; i<poly->num_triangles; i++){
        printPoly(poly->triangles+i, substr);
        strcat(str, "{");
        strcat(str, substr);
        strcat(str, "}");
        if(i<(poly->num_triangles-1)){
          strcat(str, ", ");
        }
      }
    }
  }
  if(delete_str){
    strcat(str, "\n");
    printf(str);
    free(str);
  }
}

/*
 * reverse the order of vertices
 * e.g. if poly has 3 vertices v1, v2, and v3 stored in verts
 * then verts will contain v3, v2, v1
 */
void reverseVertices(Poly *poly){
  Vertex *v1, *v2;
  int i;
  float tmp;
  if(poly==NULL){
    return;
  }
  for(v1=poly->verts, v2=poly->verts+poly->num_verts-1, i=0; i<poly->num_verts/2; i++, v1++, v2--){
    tmp=v1->x;
    v1->x=v2->x;
    v2->x=tmp;
    tmp=v1->y;
    v1->y=v2->y;
    v2->y=tmp;
  }
}

/*
 * non convex polygon poly must have 4 vertices.
 * vertex index is inside triangle formed by the three other vertices.
 * cuts poly into two triangles and update poly->triangles.
 * (the two new polygons have "sorted" vertices)
 */
void triangulate13(Poly *poly, int index, int sign){
  Poly *p;
  int i[3], j;
  if(poly==NULL){
    return;
  }
  poly->num_triangles=2;
  poly->triangles = (Poly*)calloc(poly->num_triangles, sizeof(Poly));
  for(p=poly->triangles; p<poly->triangles+2; p++){
    p->num_verts=3;
    p->verts=(Vertex*)calloc(p->num_verts, sizeof(Vertex));
    i[0]=index;
    i[1]=(i[0]+2+sign*((int)(p-poly->triangles)))%4;
    i[2]=(i[0]+2+sign*((int)(p-poly->triangles)-1))%4;
    for(j=0; j<3; j++){
      p->verts[j].x=poly->verts[i[j]].x;
      p->verts[j].y=poly->verts[i[j]].y;
    }
  }
}

/*
 * non convex polygon poly must have 4 vertices.
 * the shape of the polygon looks like two triangles so that
 * vertex index and (index+1)%4 form the base of one triangle, the two other vertices form
 * the base of the second triangle, and the top of the two triangles is the intersection of
 * the lines made with (vertex index and index+3) and (vertex index+1 and index+2).
 */
void triangulate22(Poly *poly, int index){
  Poly *p;
  Point pts[4], pt;
  int i, ii[4];
  double det, coef1, coef2;
  if(poly==NULL){
    return;
  }
  ii[0]=index;
  ii[1]=(index+3)%4;
  ii[2]=(index+1)%4;
  ii[3]=(index+2)%4;
  for(i=0; i<4; i++){
    pts[i].x=poly->verts[ii[i]].x;
    pts[i].y=poly->verts[ii[i]].y;
  }
  det=(pts[0].x-pts[1].x)*(pts[3].y-pts[2].y)-(pts[3].x-pts[2].x)*(pts[0].y-pts[1].y);
  if(det==0) //check float comparisons here...
    return;
  coef1=pts[0].x*pts[1].y-pts[1].x*pts[0].y;
  coef2=pts[2].x*pts[3].y-pts[3].x*pts[2].y;
  pt.x=((pts[0].x-pts[1].x)*coef2+(pts[3].x-pts[2].x)*coef1)/det;
  pt.y=(coef1*(pts[3].y-pts[2].y)+coef2*(pts[0].y-pts[1].y))/det;
  poly->num_triangles=2;
  poly->triangles = (Poly*)calloc(poly->num_triangles, sizeof(Poly));
  for(p=poly->triangles; p<poly->triangles+2; p++){
    p->num_verts=3;
    p->verts=(Vertex*)calloc(p->num_verts, sizeof(Vertex));
    for(i=0; i<2; i++){
      p->verts[i].x=poly->verts[ii[2*i+(int)(p-poly->triangles)]].x;
      p->verts[i].y=poly->verts[ii[2*i+(int)(p-poly->triangles)]].y;
    }
    p->verts[2].x=pt.x;
    p->verts[2].y=pt.y;
  }
}

/*
 * the polygon poly must have 4 vertices.
 * check if poly is convex. if it's not then cut it into two triangles.
 * return 0 if poly is convex or could be cut into two triangles.
 * return -1 otherwise.
 */
int convexify4(Poly *poly){
  int signs[4], nb_signs[3];
  float coef;
  Vertex *v, *vm1, *vp1;
  int i;
  if(poly==NULL || poly->num_verts!=4)
    return -1;
  nb_signs[0]=nb_signs[1]=nb_signs[2]=0;
  // compute the sign of the sine of the 4 angles to find out the geometry of the polygon
  for(i=0, v=poly->verts; i<4; i++, v++){
    vm1=i==0?v+3:v-1;
    vp1=i==3?v-3:v+1;
    coef = (vm1->x-v->x)*(vp1->y-v->y)-(vp1->x-v->x)*(vm1->y-v->y);
    if(coef>0){
      signs[i]=1;
      nb_signs[1]++;
    }
    else if(coef<0){
      signs[i]=-1;
      nb_signs[2]++;
    }
    else{
      signs[i]=0;
      nb_signs[0]++;
    }
  }
  if(nb_signs[1]==4 || nb_signs[2]==4){
    if(nb_signs[2]==4){
      reverseVertices(poly);
    }
    return 0;
  }
  // break polygon into two triangles if possible
  if((nb_signs[1]==1 && nb_signs[2]==3) || (nb_signs[1]==3 && nb_signs[2]==1)){
    int sign=nb_signs[1]==1?1:-1;
    for(i=0; i<4; i++)
      if(signs[i]==sign)
        break;
    if(i==4) // this should not happen
      return -1;
    triangulate13(poly, i, signs[i]);
    return 0;
  }
  if(nb_signs[1]==2 && nb_signs[2]==2){
    for(i=0; i<4; i++)
      if(signs[i]==1 && signs[(i+1)%4]==1)
        break;
    if(i==4){ //don't think it's possible to get here... i.e. having signs like {-1, 1, -1, 1}
      return -1;
    }
    triangulate22(poly, i);
    return 0;
  }
  // nb_signs[0]>0. poly is either a triangle, a line, or a point.
  // TODO: take care of it here...
  return -1;
}

/*
 * compute the area of the polygon poly.
 * the vertices must be "sorted" clockwise.
 * if poly is not covex, then it's triangles must be set.
 * return the area of the poly (in addition of setting the area field in poly)
 */
float computeArea(Poly *poly){
  float area=0;
  Vertex *v1, *v2;
  int i;
  if(poly == NULL){
    return -1;
  }
  // check if triangulation set, in which case set the area of each triangle and get the sum.
  if(poly->num_triangles>0){
    for(area=0, i=0; i<poly->num_triangles; i++){
      area+=computeArea(poly->triangles+i);
    }
    poly->area=area;
    return area;
  }
  for(area=0, i=0; i<poly->num_verts; i++){
    v1=poly->verts+i;
    v2=poly->verts+((i+1)%poly->num_verts);
    area+=v2->x*v1->y-v1->x*v2->y;
  }
  area/=2;
  poly->area=area;
  return area;
}

/*
 * sets the normals and dists field of poly (defining the its sides).
 * if poly is not convex, then set the normals and dists of the triangles (must be set
 * beforehand)
 * also vertices must be sorted clockwise.
 */
void setSides(Poly *poly){
  Vertex *v1, *v2;
  float norm;
  int i;
  if(poly == NULL){
    return;
  }
  if(poly->num_triangles > 0){
    for(i=0; i<poly->num_triangles; i++){
      setSides(poly->triangles+i);
    }
    return;
  }
  poly->normals=(Vector*)calloc(poly->num_verts, sizeof(Vector));
  poly->dists=(float*)calloc(poly->num_verts, sizeof(float));
  for(i=0; i<poly->num_verts; i++){
    v1=poly->verts+i;
    v2=poly->verts+((i+1)%poly->num_verts);
    norm=sqrt((v2->x-v1->x)*(v2->x-v1->x)+(v2->y-v1->y)*(v2->y-v1->y));
    poly->normals[i].x=(v1->y-v2->y)/norm;
    poly->normals[i].y=(v2->x-v1->x)/norm;
    poly->dists[i]=poly->normals[i].x*v1->x+poly->normals[i].y*v1->y;
  }
}

/*
 * returns the area of the polygons poly1 and poly2.
 * the normals and dists fields of poly1 and poly2 must be set beforehand.
 */
float intersectionArea(Poly *poly1, Poly *poly2){
  Poly poly;
  Poly *p1, *p2;
  Poly *pp1, *pp2;
  Vector *normals, *n1, *n2, *n_end;
  Vertex *v;
  Point pt;
  float *dists, *angles;
  int num_p1, num_p2, num_dists, verts_allocated;
  float area, det, tmp;
  int i, j, k;
  if(poly1==NULL || poly2==NULL){
    return -1;
  }
  // set array of polynoms representing poly1
  if(poly1->num_triangles==0){
    p1=poly1;
    num_p1=1;
  }
  else{
    p1=poly1->triangles;
    num_p1=poly1->num_triangles;
  }
  // set array of polynoms representing poly2
  if(poly2->num_triangles==0){
    p2=poly2;
    num_p2=1;
  }
  else{
    p2=poly2->triangles;
    num_p2=poly2->num_triangles;
  }
  verts_allocated=4;
  poly.num_verts=0;
  poly.verts=(Vertex*)calloc(verts_allocated, sizeof(Vertex));
  poly.area=0;
  poly.num_triangles=0;
  poly.triangles=NULL;
  poly.dists=NULL;
  poly.normals=NULL;
  num_dists=0;
  normals=NULL;
  dists=NULL;
  angles=(float*)calloc(verts_allocated, sizeof(float));
  area=0;
  // check intersection of each pair of polygons
  for(pp1=p1; pp1<p1+num_p1; pp1++){
    for(pp2=p2; pp2<p2+num_p2; pp2++){
      // make set of normals and distances.
      if(pp1->num_verts+pp2->num_verts>num_dists){
        num_dists=pp1->num_verts+pp2->num_verts;
        normals=(Vector*)realloc(normals, num_dists*sizeof(Vector));
        dists=(float*)realloc(dists, num_dists*sizeof(float));
      }
      for(i=0; i<pp1->num_verts; i++){
        normals[i].x=pp1->normals[i].x;
        normals[i].y=pp1->normals[i].y;
        dists[i]=pp1->dists[i];
      }
      for(i=0, j=pp1->num_verts; i<pp2->num_verts; i++, j++){
        normals[j].x=pp2->normals[i].x;
        normals[j].y=pp2->normals[i].y;
        dists[j]=pp2->dists[i];
      }
      // for each pair of lines, compute intersection
      n_end=normals+pp1->num_verts+pp2->num_verts;
      poly.num_verts=0;
      for(n1=normals, i=0; n1<n_end-1; n1++, i++){
        for(n2=n1+1, j=i+1; n2<n_end; n2++, j++){
          det=n1->x*n2->y-n2->x*n1->y;
          if(fabs(det)<SMALL)
            continue;
          pt.x=(n2->y*dists[i]-n1->y*dists[j])/det;
          pt.y=(n1->x*dists[j]-n2->x*dists[i])/det;
          // check if intersection is inside pp1 and pp2
          for(k=0; k<pp1->num_verts; k++){
            if((pt.x*pp1->normals[k].x+pt.y*pp1->normals[k].y)>(pp1->dists[k]+10*SMALL))
              break;
          }
          if(k<pp1->num_verts)
            continue;
          for(k=0; k<pp2->num_verts; k++){
            if((pt.x*pp2->normals[k].x+pt.y*pp2->normals[k].y)>(pp2->dists[k]+10*SMALL))
              break;
          }
          if(k<pp2->num_verts)
            continue;
          // point pt is inside both pp1 and pp2: this point is a vertex.
          if(poly.num_verts>=verts_allocated){
            verts_allocated*=2;
            poly.verts=(Vertex*)realloc(poly.verts, verts_allocated*sizeof(Vertex));
            angles=(float*)realloc(angles, verts_allocated*sizeof(float));
          }
          poly.verts[poly.num_verts].x=pt.x;
          poly.verts[poly.num_verts].y=pt.y;
          poly.num_verts++;
        }
      }
      // order vertices clockwise (get point inside poly, compute bunch of angles and sort them)
      pt.x=pt.y=0;
      for(v=poly.verts; v<poly.verts+poly.num_verts; v++){
        pt.x+=v->x;
        pt.y+=v->y;
      }
      pt.x/=poly.num_verts;
      pt.y/=poly.num_verts;
      for(i=0, v=poly.verts; i<poly.num_verts; i++, v++){
        angles[i]=(float)acos((v->x-pt.x)/sqrt((v->x-pt.x)*(v->x-pt.x)+(v->y-pt.y)*(v->y-pt.y)));
        if((v->y-pt.y)<0)
          angles[i]=-angles[i];
      }
      for(i=0; i<poly.num_verts-1; i++){
        for(j=i+1; j<poly.num_verts; j++){
          if(angles[i]<angles[j]){
            tmp=angles[i];
            angles[i]=angles[j];
            angles[j]=tmp;
            tmp=poly.verts[i].x;
            poly.verts[i].x=poly.verts[j].x;
            poly.verts[j].x=tmp;
            tmp=poly.verts[i].y;
            poly.verts[i].y=poly.verts[j].y;
            poly.verts[j].y=tmp;
          }
        }
      }
      area+=computeArea(&poly);
    }
  }
  free(poly.verts);
  free(angles);
  free(normals);
  free(dists);
  return area;
}

/*
 * transformation of a ray of light through a lens located over the screen (position of the center
 * of the lens=(0, 0, heighth+radius)).
 * lens parameters (radius, height, indices, etc.) are stored in params.
 * ptin is the point the viewer is looking at. ptout is the point on the screen the viewer will
 * actually see due to the deflection of the ray through the lens.
 * return -1 if an error occured, 1 if the point is inside the circle (=lens projected on screen),
 * 0 if the point is outside (i.e. ptin=ptout)
 */
int lensTransform(void *params, Point *ptin, Point *ptout){
  float theta1, theta2, theta3, theta4;
  float l1, l2, l3, l4;
  float coef1, coef2, coef3;
  float s_radius; //radius of the sphere used to make the lens
  LensParams *lp;
  if(params==NULL || ptin==NULL || ptout==NULL){
    return -1;
  }
  lp=(LensParams*)params;
  l1=sqrt(ptin->x*ptin->x+ptin->y*ptin->y);
  if(l1>lp->radius)
    return 0;
  s_radius=(lp->radius*lp->radius+lp->thickness*lp->thickness)/(2*lp->thickness);
  switch(lp->lens_type){
    case ONE_CURVED_SIDE_TOP:
      theta1=(float)asin(l1/s_radius);
      theta2=(float)asin((l1/s_radius)*(lp->n1/lp->n2));
      theta3=(float)asin((lp->n2/lp->n1)*sin(theta2-theta1));
      l3=l1+lp->height*tan(theta3)+(lp->thickness-s_radius*(1-cos(theta1)))*tan(theta2-theta1);
      break;
    case ONE_CURVED_SIDE_BOTTOM:
      theta1=(float)asin(l1/s_radius);
      theta3=(float)asin((lp->n2/lp->n1)*l1/s_radius);
      l3=l1+(lp->height+s_radius*(1-cos(theta1)))*tan(theta1-theta3);
      break;
    case TWO_CURVED_SIDES:
      theta1=(float)asin(l1/s_radius);
      theta2=(float)asin((l1/s_radius)*(lp->n1/lp->n2))-theta1;
      l2=l1+(lp->thickness-s_radius*(1-cos(theta1)))*tan(theta2);
      coef1=l2*sin(theta2)+(s_radius-lp->thickness)*cos(theta2);
      coef2=lp->radius*lp->radius-l2*l2;
      coef3=sqrt(coef1*coef1+coef2)-coef1;
      l3=l2+sin(theta2)*coef3;
      theta3=(float)asin(l3/s_radius);
      theta4=(float)asin((sin(theta3-theta2)*(lp->n2/lp->n1)));
      l4=l3+(lp->height+lp->thickness-cos(theta2)*coef3)*tan(theta3-theta4);
      l3=l4;
      break;
    case SPHERE:
      theta1=(float)asin(l1/lp->radius);
      theta2=(float)asin((l1/lp->radius)*(lp->n1/lp->n2));
      l2=lp->radius*sin(2*theta2-theta1)+
          (lp->height+lp->radius*(1-cos(2*theta2-theta1)))*tan(2*(theta2-theta1));
      l3=l2;
      break;
    default:
      l3=l1;
      break;
  }
  if(fabs(l3)<10*SMALL){
    ptout->x=ptout->y=0;
  }
  else{
    ptout->x=ptin->x*(l3/l1);
    ptout->y=ptin->y*(l3/l1);
  }

  return 1;
}

/*
 *
 */
Filter* getFilter(int (*transform)(void*, Point*, Point*), void *params, int dim_x, int dim_y,
                  Point pt0){
  Filter *tg;
  Point *t_grid; //transformed grid
  Poly *cell_poly, *poly;
  Vertex *v;
  char *is_transformed;
  int dim;
  float min_x, max_x, min_y, max_y;
  float area;
  int i, j, k;
  int ii, jj;
  int index;
  int last, allocated;

  if(params==NULL)
    return NULL;
  tg=(Filter*)calloc(1, sizeof(Filter));
  tg->dim_x=dim_x;
  tg->dim_y=dim_y;
  tg->x0=pt0.x;
  tg->y0=pt0.y;
  tg->dx=tg->dy=1;
  dim=tg->dim_x*tg->dim_y;
  tg->num_pts=(int*)calloc(dim, sizeof(int));
  tg->pts=(Point**)calloc(dim, sizeof(Point*));
  tg->coefs=(float**)calloc(dim, sizeof(float*));

  // compute the transformed grid (used to get the polygons corresponding to the cells.
  t_grid=(Point*)calloc((tg->dim_x+1)*(tg->dim_y+1), sizeof(Point));
  is_transformed=(char*)calloc((tg->dim_x+1)*(tg->dim_y+1), sizeof(char));
  for(k=0, j=0; j<=tg->dim_y; j++){
    for(i=0; i<=tg->dim_x; i++, k++){
      t_grid[k].x=tg->x0+i*tg->dx-tg->dx/2;
      t_grid[k].y=tg->y0+j*tg->dy-tg->dy/2;
      is_transformed[k]=(char)transform(params, t_grid+k, t_grid+k);
      if(is_transformed[k]==-1)
        return NULL; //should never happen!
    }
  }

  //for each cell, get the corresponding polygon, transform it, get all the overlapping cells,
  //(number of overlaping cells=numPoints[k]) and compute percentage of area inside each of these
  // cells
  cell_poly=generatePoly(4);
  cell_poly->verts[0].x=0;
  cell_poly->verts[0].y=1;
  cell_poly->verts[1].x=1;
  cell_poly->verts[1].y=1;
  cell_poly->verts[2].x=1;
  cell_poly->verts[2].y=0;
  cell_poly->verts[3].x=0;
  cell_poly->verts[3].y=0;
  setSides(cell_poly); // verts won't be used from now on, just normals and dists. only update dists
  for(index=0, k=0, j=0; j<tg->dim_y; j++, index++){
    for(i=0; i<tg->dim_x; i++, k++, index++){
      if(!(is_transformed[index] &&
           is_transformed[index+1] &&
           is_transformed[index+tg->dim_x+1] &&
           is_transformed[index+tg->dim_x+2]))
        continue;
      if(!is_transformed[index] ||
         !is_transformed[index+1] ||
         !is_transformed[index+tg->dim_x+1] ||
         !is_transformed[index+tg->dim_x+2]){
        // cell contains border of transformation area.
        // TODO: handle this case
        continue;
      }
      poly=generatePoly(4);
      v=poly->verts;
      v[0].x=t_grid[index+tg->dim_x+1].x;
      v[0].y=t_grid[index+tg->dim_x+1].y;
      v[1].x=t_grid[index+tg->dim_x+2].x;
      v[1].y=t_grid[index+tg->dim_x+2].y;
      v[2].x=t_grid[index+1].x;
      v[2].y=t_grid[index+1].y;
      v[3].x=t_grid[index].x;
      v[3].y=t_grid[index].y;
      if(convexify4(poly)==-1){
        deletePoly(poly, 1);
        continue;
      }
      setSides(poly);
      computeArea(poly);
      min_x=max_x=v[0].x;
      min_y=max_y=v[0].y;
      for(ii=1; ii<4; ii++){
        if(v[ii].x<min_x)
          min_x=v[ii].x;
        else if(v[ii].x>max_x)
          max_x=v[ii].x;
        if(v[ii].y<min_y)
          min_y=v[ii].y;
        else if(v[ii].y>max_y)
          max_y=v[ii].y;
      }
      min_x=floorf(min_x);
      max_x=ceilf(max_x);
      min_y=floorf(min_y);
      max_y=ceilf(max_y);
      allocated=0;
      if((max_x-min_x) > 100 ||(max_y-min_y) > 100){
        tg->num_pts[k]=1;
        tg->coefs[k]=(float*)calloc(1, sizeof(float));
        tg->pts[k]=(Point*)calloc(1, sizeof(Point));
        tg->coefs[k][0]=1;
        tg->pts[k][0].x=floor((max_x-min_x)/2); //TODO: fix that!
        tg->pts[k][0].y=floor((max_y-min_y)/2);
        deletePoly(poly, 1);
        continue;
      }
      for(jj=(int)min_y; jj<(int)max_y; jj++){
        for(ii=(int)min_x; ii<(int)max_x; ii++){
          cell_poly->dists[0]=jj+1;
          cell_poly->dists[1]=ii+1;
          cell_poly->dists[2]=-jj;
          cell_poly->dists[3]=-ii;
          area=intersectionArea(poly, cell_poly);
          if(area<SMALL)
            continue;
          last=tg->num_pts[k]++;
          if(tg->num_pts[k]>allocated){
            allocated+=32;
            tg->coefs[k]=(float*)realloc(tg->coefs[k], allocated*sizeof(float));
            tg->pts[k]=(Point*)realloc(tg->pts[k], allocated*sizeof(Point));
          }
          tg->coefs[k][last]=area/poly->area;
          tg->pts[k][last].x=ii;
          tg->pts[k][last].y=jj;
        }
      }
      tg->coefs[k]=(float*)realloc(tg->coefs[k], tg->num_pts[k]*sizeof(float));
      tg->pts[k]=(Point*)realloc(tg->pts[k], tg->num_pts[k]*sizeof(Point));
      deletePoly(poly, 1);
    } //for(i)
  } //for(j)
  free(t_grid);
  free(is_transformed);
  return tg;
}

Filter* getSimpleFilter(int (*transform)(void*, Point*, Point*), void *params, int dim_x, int dim_y,
                        Point pt0){
  Filter *tg;
  Point *t_grid; //transformed grid
  char *is_transformed;
  int dim;
  int i, j, k;

  if(params==NULL)
    return NULL;
  tg=(Filter*)calloc(1, sizeof(Filter));
  tg->dim_x=dim_x;
  tg->dim_y=dim_y;
  tg->x0=pt0.x;
  tg->y0=pt0.y;
  tg->dx=tg->dy=1;
  dim=tg->dim_x*tg->dim_y;
  tg->num_pts=(int*)calloc(dim, sizeof(int));
  tg->pts=(Point**)calloc(dim, sizeof(Point*));
  tg->coefs=(float**)calloc(dim, sizeof(float*));

  // compute the transformed grid (used to get the polygons corresponding to the cells.
  t_grid=(Point*)calloc(tg->dim_x*tg->dim_y, sizeof(Point));
  is_transformed=(char*)calloc(tg->dim_x*tg->dim_y, sizeof(char));
  for(k=0, j=0; j<tg->dim_y; j++){
    for(i=0; i<tg->dim_x; i++, k++){
      t_grid[k].x=tg->x0+i*tg->dx;
      t_grid[k].y=tg->y0+j*tg->dy;
      is_transformed[k]=(char)transform(params, t_grid+k, t_grid+k);
      if(is_transformed[k]==-1)
        return NULL; //should never happen!
    }
  }

  for(k=0, j=0; j<tg->dim_y; j++){
    for(i=0; i<tg->dim_x; i++, k++){
      if(!is_transformed[k])
        continue;
      tg->num_pts[k]=1;
      tg->coefs[k]=(float*)calloc(1, sizeof(float));
      tg->pts[k]=(Point*)calloc(1, sizeof(Point));
      tg->coefs[k][0]=1;
      tg->pts[k][0].x=floor(t_grid[k].x);
      tg->pts[k][0].y=floor(t_grid[k].y);
    }
  }
  free(t_grid);
  free(is_transformed);
  return tg;
}

/*
 * release memory allocated for Filter.
 */
void freeFilter(Filter *filter){
  int k, dim;
  if(filter==NULL){
    return;
  }
  dim=filter->dim_x*filter->dim_y;
  for(k=0; k<dim; k++){
    free(filter->pts[k]);
    free(filter->coefs[k]);
  }
  free(filter->num_pts);
  free(filter->pts);
  free(filter->coefs);
  free(filter);
}

Filter* getLensFilter(LensParams *lp){
  int dim = 2*lp->radius+1;
  Point pt0;
  pt0.x=-lp->radius;
  pt0.y=-lp->radius;
  return getFilter(lensTransform, (void*)lp, dim, dim, pt0);
}

Filter* getSimpleLensFilter(LensParams *lp){
  int dim = 2*lp->radius+1;
  Point pt0;
  pt0.x=-lp->radius;
  pt0.y=-lp->radius;
  return getSimpleFilter(lensTransform, (void*)lp, dim, dim, pt0);
}
