/*ab model mer*/
#include<stdio.h>
#include<string.h>
#include <math.h>
#define MEMERR {fprintf(stderr,"\nNot sufficient memory\n");exit(1);}
#define MAXLENGTH 245 
#define ATOM MAXLENGTH 
#define NATOM 128
#define NTYPE 245
#define NRES NATOM 
#define MAXPAIR 275
#define MAXHBPAIR 81
#define urea 1
#define nourea 0
#define UREA nourea
#define RG 11.485
#define MG 3250.0

#define K 0.6/300


struct config{
double x;
double y;
double z;
};

struct res{
char type[4];
char type1[4];
char type2[4];
float eps;
char name[2][4];
int atomindex[2];
int typeindex[2];
int num;
struct config coor[2];
};


struct res native[MAXLENGTH];
/*struct config polymer[MAXLENGTH];*/
struct res maparray[210]; 
struct config polymer[MAXLENGTH];
struct config polymer1[MAXLENGTH];
