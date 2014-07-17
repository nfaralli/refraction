#ifndef __UTILS_H__
#define __UTILS_H__

#define SMALL 1e-6

typedef struct Point_ {
  float x;
  float y;
} Point;
typedef Point Vertex;
typedef Point Vector;

/*
 * structure defining a polygon with num_verts vertices.
 * verts[i] is connected to verts[i+1] and verts[num_verts-1]
 * is connected to verts[0].
 */
typedef struct Poly_ {
  int num_verts;
  Vertex *verts;
  Vector *normals; // normals[i] and dists[i] define the side of the polygon between vertices i
  float *dists;    // and i+1. normals[i] is perpendicular to the side and is a unit vector.
  float area;
  struct Poly_ *triangles; //if polygon is not convex, contains triangulation of the polygon.
  int num_triangles;
} Poly;

typedef enum {
  ONE_CURVED_SIDE_TOP, // lens has one flat side, curved side located on top
  ONE_CURVED_SIDE_BOTTOM, // lens has one flat side, curved side located on bottom
  TWO_CURVED_SIDES, // both sides of lens are curved
  SPHERE, // lens is actually a sphere (<=> TWO_CURVED_SIDES where thickness=radius)
  NUM_LENS_TYPE
} LensType;

typedef struct LensParams_ {
  int radius; //lens radius in pixels
  float thickness; //lens thickness (half of thickness for TWO_CURVED_SIDES)
  float height; //height from the screen to the bottom of the sphere
  LensType lens_type;
  float n1; //index of air
  float n2; //index of glass
} LensParams;

typedef struct Filter_ {
  int dim_x, dim_y;
  float x0, y0; //coordinate of top left corner cell
  float dx, dy; //coordinate of cell (i, j) = (x0+i*dx, y0+j*dy), 0<=i<dim_x, 0<=j<dim-y
  int *num_pts; //array of dim_x*dim_y. numPoints[k]=0 <=> corresponding point isn't transformed
  Point **pts; //size of pts[k] = num_pts[k].
  float **coefs; //size of pts[k] = num_pts[k].
} Filter;

Filter* getLensFilter(LensParams *lp);
Filter* getSimpleLensFilter(LensParams *lp);
void freeFilter(Filter*);

#endif // __UTILS_H__
