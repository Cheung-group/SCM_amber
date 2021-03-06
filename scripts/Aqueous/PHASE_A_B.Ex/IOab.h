/*ab model mer*/
#include<stdio.h>
#include<string.h>
#include <math.h>
#define MEMERR {fprintf(stderr,"\nNot sufficient memory\n");exit(1);}
#define MAXLENGTH  245 
#define ATOM MAXLENGTH 
#define NATOM  128
#define NTYPE  245
#define NRES NATOM 
#define MAXPAIR 275 /*132+hb4*/            /*total number of native contacts to be analyzed*/
/*#define MAXNNPAIR 2080*. /*132+hb4*/            /*total number of nn native contacts to be analyzed*/
#define MAXHBPAIR 81
#define STEP 5 
#define TERSTEP 5 
/*#define MAXBOX  1358*/
#define K 0.6/300

#define THRESLENGTH 1.94

struct config{
double x;
double y;
double z;
};

struct res{
char type[4];
char name[2][4];
int atomindex[2];
int typeindex[2];
int num;
int flag;
struct config coor[2];
double side;
};


struct res native[MAXLENGTH];
/*struct config polymer[MAXLENGTH];*/

struct config polymer[MAXLENGTH];
struct config polymer1[MAXLENGTH];
struct config npolymer[MAXLENGTH];
struct config centercoor[1];


FILE *openfile(char *filename, char *mode);
int nextline(FILE *fp);
int get_Ntot(FILE *fp);



int get_Ntot(FILE *fp)
{

        /* Get total number of atoms. */
        int Ntot=0;
        int Q;

        do{
        fscanf(fp,"%d",&Q);
        Ntot++;
        }while( (nextline(fp)) != EOF );

        rewind(fp);

        printf("\nTotal number of Q read: %d.\n",Ntot-1);

        return Ntot-1;
}



int nextline(FILE *fp)
{
        int ch;
        while ((ch=fgetc(fp)) != '\n'){ if (ch == EOF) return EOF;}
        return 1;
}

FILE *openfile(char *filename, char *mode)
{
        FILE *fp;

        if((fp=fopen(filename, mode))==NULL){
        fprintf(stderr,"\nCan't open file %s.\n", filename);
        exit(1);}
        return fp;
}

