/* C program that reads in a ppm P3 image file and extracts the  */
/* relevant data points from what it thinks might be a graph */
/* By Joshua Miller */
/* For use and reproduction by Andy Gillespie */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <getopt.h>
#define CUTOFF_LEVEL 150
typedef unsigned char byte;

/* PPM Image structure */
typedef struct in_memory_ppm {
  unsigned int w;
  unsigned int h;
  byte *reds;
  byte *greens;
  byte *blues;
} ppm;

typedef struct point_linked_list PLL;

struct point_linked_list{
  double x;
  double y;
  PLL * next;
};

/* Add plot point to list of found points */
PLL* pll_add(PLL * list, double x, double y){
  PLL *new = (PLL*) malloc(sizeof(PLL));
  new->x = x;
  new->y = y;
  new->next = list;
  return new;
}


void usage(){

  fprintf(stderr,"\nOGR is a software for extracting the data from data plots. ");
  fprintf(stderr,"Reads in an uncompressed PPM image and prints data points to stdout\n");
  fprintf(stderr,"\t\t\t----  By Joshua Miller ----\n");
  fprintf(stderr, "usage:\n");
  fprintf(stderr, "\togr [-x x][-X X][-y y][-Y Y][-r R][-o > imageOUT.ppm] < imageIN.ppm \n\n");
  fprintf(stderr, "\t\t <\t  imageIN.ppm (P3) is scanned plot to extract from\n\n");
  fprintf(stderr, "Optional arguments:\n");
  fprintf(stderr, "\t\t-r\tR is the radius of each data point in pixels [default:8]\n");
  fprintf(stderr, "\t\t-x\tx is the lowerbound x scale [default:0]\n");
  fprintf(stderr, "\t\t-X\tX is the upperbound x scale [default:1]\n");
  fprintf(stderr, "\t\t-y\ty is the lowerbound y scale [default:0]\n");
  fprintf(stderr, "\t\t-Y\tY is the upperbound y scale [default:1]\n");
  fprintf(stderr, "\t\t-o\t  Output post-processed image include [> imageOUT.ppm]\n");
  
  exit(2);
}

/* Reads ppm from stdin */
ppm* ppm_in(){

  /* Initialize ppm on heap */
  ppm *img = (ppm*) malloc(sizeof(ppm));

  /* Pass over ppm format */
  scanf("P3\n");

  /* Ignore comments */
  char c = getchar();
  while (c == '#') {
    while (getchar() != '\n') ;
    c = getchar();
  }
  ungetc(c, stdin);

  /* Read in dimensions */
  scanf("%u %u\n", &img->w, &img->h);
  if (img->w ==0 || img->h == 0){
    fprintf(stderr,"\nProblem reading in ppm image.\n");
    usage();
  }
  int s = img->w*img->h;


  /* Pass over color height */
  scanf("255\n");
  unsigned int i,r,g,b;

  /* Initialize RGB arrays */
  img->reds = (byte*) malloc(s*sizeof(byte));
  img->blues = (byte*) malloc(s*sizeof(byte));
  img->greens = (byte*) malloc(s*sizeof(byte));

  /* Read in RGB data */
  for (i=0; i < s; i++){
    scanf("%d %d %d\n", &r, &g, &b);
    img->reds[i] = (byte) r;
    img->greens[i] = (byte) g;
    img->blues[i] = (byte) b;
  }

  return img;

}

void ppm_out(ppm *img){

  /* Print formatting & size */
  printf("P3\n%d %d\n255\n", img->w,img->h);
  int s = img->w*img->h;
  int i = 0;
  
  /* Print RGB data */
  for (i=0; i < s; i++){
    printf("%d\t%d\t%d\n",img->reds[i], img->greens[i], img->blues[i]);
  }

}

ppm *ppm_clone(ppm *img){

  /* Allocate clone memory */
  ppm *clone = (ppm*) malloc(sizeof(ppm));
  int s = img->w*img->h;

  clone->w = img->w;
  clone->h = img->h;
  
  /* Initialize RGB arrays */
  clone->reds = (byte*) malloc(s*sizeof(byte));
  clone->blues = (byte*) malloc(s*sizeof(byte));
  clone->greens = (byte*) malloc(s*sizeof(byte));

  /* Copy RGB values */
  int i;
  for (i=0; i < s; i++){
    clone->reds[i] = img->reds[i];
    clone->greens[i] = img->greens[i];
    clone->blues[i] = img->blues[i];
  }
  return clone;
}

/* Frees PPM memory allocation */
void ppm_free(ppm *img){

  /* Free color arrays */
  free(img->reds);
  free(img->greens);
  free(img->blues);
 
  /* Free ppm pointer */
  free(img);
  return;
}

byte lum(byte r, byte g, byte b){
  return (byte) 0.2126*r + 0.7152*g + 0.0722*b;
}

/* Reduce noise for better edge detection */
void noise_reduce(ppm *img){
  int w = img->w;
  int h = img->h;
  int s = w*h;
  int i,x,y;
  int r = 3;
  long c = 0;

  /* Level the image with a cuttoff luminocity */
  for (i=r*w; i < s-r*w; i++){
    c = 0;
    int red = img->reds[i];
    if (red != 255){
      for (x = -r; x < r; x ++){
	for (y = -r; y < r; y ++){
	  c += img->reds[i+y*w + x];
	}
      }
    }

    img->reds[i] = img->greens[i] = img->blues[i] = (int) c/r/r;
  }
  return;
}

/* Imparts threshold contrast for edge detection */
void cutoff(ppm *img){

  byte r,g,b;
  int s = img->w*img->h;
  int i;

  /* Level the image with a cuttoff luminocity */
  for (i=0; i < s; i++){
    byte m = lum(img->reds[i], img->greens[i], img->blues[i]);
    if (m > CUTOFF_LEVEL){
      img->reds[i] = img->greens[i] = img->blues[i] = 0;
    } else { 
      img->reds[i] = img->greens[i] = img->blues[i] = 255;
    }

  }
  return;
}

/* Detects edges by gradient */
void edge_detect(ppm *img){
  cutoff(img);
  byte r,g,b;
  int w = img->w;
  int h = img->h;
  int s = w*h;
  int i;

  for (i=w; i < s-w; i++){
    int g1 = img->reds[i] - img->reds[i+1];
    int g2 = img->reds[i] - img->reds[i+w];
    int g = (abs(g1)>abs(g2)) ? g1 : g2;
    if (g>0){
      img->blues[i] = 255;
      img->reds[i] = img->greens[i] = 0;
    } else if (g< 0){
      img->greens[i] = 255;
      img->reds[i] = img->blues[i] = 0;
    } else {
      img->reds[i] = img->greens[i] = img->blues[i] = 0;
    }

  }

  return;
}

/* Offsets circle for hough transform */
void draw_circle(ppm *hough, int i, double R){
  double DEG_TO_RAD = 0.0174532925;
  int w = hough->w;
  int s = w*hough->h;
  int deg;
  for (deg = 0; deg < 360; deg ++){
    int x = R*cos(deg*DEG_TO_RAD);
    int y = R*sin(deg*DEG_TO_RAD);
    if (i + x >= 0 && i + x < s){
      if (i + y*w >= 0 && i + y*w  < s){
	if (hough->reds[i+x+y*w] < 255) hough->reds[i+x+y*w]++;
	if (hough->greens[i+x+y*w] < 255) hough->greens[i+x+y*w]++;
	if (hough->blues[i+x+y*w] < 255) hough->blues[i+x+y*w]++;
      }
    }
  }

  return;
}

/* Is the point a local max? */
int is_local_max(ppm *img, int i, int w){
  if (img->reds[i] < img->reds[i-1]) return 0;
  if (img->reds[i] < img->reds[i+1]) return 0;
  if (img->reds[i] < img->reds[i-w]) return 0;
  if (img->reds[i] < img->reds[i+w]) return 0;
  return 1;
}

/* Routine for automatic radius selection */
int count_circle(ppm* hough, int i, double R){
  int w = hough->w;
  double DEG_TO_RAD = 0.0174532925;
  int s = w*hough->h;
  int deg;
  int count = 0;
  for (deg = 0; deg < 360; deg ++){
    int x = R*cos(deg*DEG_TO_RAD);
    int y = R*sin(deg*DEG_TO_RAD);
    if (i + x >= 0 && i + x < s){
      if (i + y*w >= 0 && i + y*w  < s){
	count += hough->blues[i+x+y*w];
      }
    }
  }

  return count;
  
  
  
}

/* Gets data points from graph */
ppm *find_points(ppm *img, double r, PLL** centers){

  int w = img->w;
  int h = img->h;
  int s = w*h;
  int i;
  double maxR = 2;
  ppm *hough = ppm_clone(img);
  int max = 0;
  int count = 0;
  int maxAve = 0;

 /* TODO: Inclusion of automatic radius selection */

  /* { */
  /* /\* for (r = 2; r < 15; r += .5){ *\/ */
  /*   bzero(hough->reds, s*sizeof(byte)); */
  /*   bzero(hough->greens, s*sizeof(byte)); */
  /*   bzero(hough->blues, s*sizeof(byte)); */

  /*   ppm *h2 = ppm_clone(hough); */

  /*   for (i=w; i < s-w; i++){ */
  /*     int b = img->blues[i]; */
  /*     int g = img->greens[i]; */
  /*     if (b == 255 || g == 255) */
  /* 	draw_circle(hough, i, r); */
  /*   } */
    
  maxR = r;
  bzero(hough->reds, s*sizeof(byte));
  bzero(hough->greens, s*sizeof(byte));
  bzero(hough->blues, s*sizeof(byte));

  ppm *h2 = ppm_clone(hough);

  for (i=w; i < s-w; i++){
    int b = img->blues[i];
    int g = img->greens[i];
    if (b == 255 || g == 255)
      draw_circle(hough, i, maxR);
  }

  max = 0;
  for (i=w; i < s-w; i++){
    max = (hough->reds[i] > max) ? hough->reds[i] : max;
  }

  for (i=w; i < s-w; i++){
    if (hough->reds[i] > max*.8 && is_local_max(hough, i,w)){
      h2->reds[i] = h2->greens[i] = h2->blues[i]  =  255;
      *centers = pll_add(*centers, i % w, i/w);
    } else {
      h2->reds[i] = h2->greens[i] = h2->blues[i]  =  0;
    }
  }

  return hough;
}

/* main int,char *array -> int */
int main(int argc, char *argv[]){
  double r = 8;
  int opt;
  int out = 0;
  double lowerX = 0;
  double lowerY = 0;
  double upperX = 1;
  double upperY = 1;
  while ((opt = getopt (argc, argv, "x:X:y:Y:hor:")) != -1)
    switch (opt){
    case 'r':
      r = atof(optarg);
      fprintf(stderr, "Expected radius: %f\n",r);
      break;
    case 'o':
      out = 1;
      break;
    case 'h':
      usage();
      break;
    case 'x':
      lowerX = atof(optarg);
      fprintf(stderr, "Lowerbound x: %f\n", lowerX);
      break;
    case 'X':
      upperX = atof(optarg);
      fprintf(stderr, "Upperbound x: %f\n", upperX);
      break;
    case 'y':
      lowerY = atof(optarg);
      fprintf(stderr, "Lowerbound y: %f\n", lowerY);
      break;
    case 'Y':
      upperY = atof(optarg);
      fprintf(stderr, "Upperbound y: %f\n", upperY);
      break;
    default:
      fprintf(stderr, "Unknown command line arg. -h for help.");
      usage();
      exit(1);
    }

  ppm* img = ppm_in();
  edge_detect(img);

  PLL* centers = NULL;
  ppm *hough = find_points(img, r, &centers);
  int w = hough->w;
  int h = hough->h;
  while (centers){
    double x = centers->x/w*(upperX-lowerX)+lowerX;
    double y = (h -centers->y)/h*(upperY-lowerY)+lowerY;
    fprintf(stderr,"%f\t%f\n",x,y);
    centers = centers->next;
  }
  if (out) ppm_out(hough);

  ppm_free(img);
  
}

