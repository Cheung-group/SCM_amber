/*****************************************
frcfield_ab.c

This is the program for building the frcfield for explicit beta model

the command line: a.out input.crd typeindex pairindex


9.17.01

This code can be used to generate force field for (Bulk, Sph and DB crowders, and urea)

****************************************/

#include "IOablinux.h"
#include "ctype.h"
#define PI 180
#define LJ 6 
#define YES 1
#define NO 0 
#define VdW10_12 NO
#define HB10_12 NO 
#define HETEROGO NO 
#define CONVERT 180/3.1415927

#define SIX(a) a*a*a*a*a*a
#define FIV(a) a*a*a*a*a
#define TEN(a) FIV(a)*FIV(a)
#define TWV(a) SIX(a)*SIX(a)
#define HB 0.6 
#define HHB 0.6 
#define GoEn 0.6 
#define HGoEn 0.6 

struct list{
  int index;
  int a1;
  int a2;
  int t1;
  int t2;
  int res1;
  int res2;
  double dist;
  char type[4];
  int minicore;
 };

struct list GoPair[MAXPAIR];
int nextline(FILE *fp); 
int get_Ntot(FILE *fp);
FILE *openfile(char *filename, char *mode);

double di,modulo;
int di1;
char *inputfile,*indexfile,*pairfile,*mapfile,*disulfile, corefile[50];
FILE *fp1,*fp2,*fp3,*fp4,*fp5,*fp6,*fp7;
int Ntot, Nds;
int dsA[NATOM], dsB[NATOM];

struct config vec,vec1, vec2;

/*even number; cak odd number, cb*/

void read_crd(FILE *fp);
void read_mapfile(void);
void read_disulfile(void);
void read_atomindex(void);
void read_pairindex(FILE *fp,struct list Pair[MAXPAIR]);
double distance(struct config poly1, struct config poly2);
double calan(struct config A1, struct config A2, struct config A3);
void caldi(double *di,struct config A1,struct config A2,struct config A3,struct config A4);
void checkatomtype(char c[4],int *t);

char xx[2];
char xy[2];

main(int argc, char *argv[])
{

int ci1, ci2, ci3, ci4, ia, ja;
int i,j,k,N,m,n,l;
int i1,i2,i3,i4;
int a;
int aa, bb;
char c1[4],c2[4],type[4];
float epsilon1,epsilon2,epsilon3,epsilon4,epsilon5,epsilon6,av_epsilon;
double zero;
double A,B,dummy;
if(argc<6){printf("Usage: a.out input.crd typeindex(INDEX_ATOM) pairindex(INDEX_allres.pair) mapfile disulfile\n"); exit(1);}
/*typeindex(INDEX_atom) is the type of residue.
if it is "a" for gly, "b" for non-gly*/
/*pairindex (INDEX_allres.pair)*/

        inputfile=argv[1];
        indexfile=argv[2];
        pairfile=argv[3];
        mapfile=argv[4];
        disulfile=argv[5];

        fp1=openfile(inputfile,"rt");
        fp2=openfile("frc.go","wt");
        fp3=openfile(indexfile,"rt");
        fp4=openfile(pairfile,"rt");
        fp6=openfile(mapfile,"rt");
        fp7=openfile(disulfile,"rt");

	
        /*outfile file*/

 for (i=0;i<NATOM;i++) {

        for(k=0;k<2;k++)
        {
        if(k==0)
                {
                ia=(int)(i/36);
                ja=(i%36);
                ci1=65+ia;
                if(ja<10) ci2=48+ja;
                else      ci2=65+ja-10;
	
                sprintf(native[i].name[0],"%c%c",ci1,ci2); /*ca*/ 
	 if (strcmp(native[i].name[0],"DU")==0)	
                {sprintf(native[i].name[0],"Z9"); /*ca*/ }
	 if (strcmp(native[i].name[0],"EP")==0)	
                {sprintf(native[i].name[0],"ZP"); /*ca*/ }
		
                }
        if(k==1)
                {
		ia=(int)(i/36);
                ja=(i%36);
                ci3=78+ia;
                if(ja<10) ci4=48+ja;
                else      ci4=65+ja-10;
                sprintf(native[i].name[1],"%c%c",ci3,ci4);
                }

	}
  }

	/*crowding agent is BY*/

                sprintf(xx,"ZY");
                sprintf(xy,"ZZ");


       /*read index crd*/
        read_atomindex();
	
	/* read BTmap file   */
	read_mapfile();

	/* read disulfide bond list file   */
	read_disulfile();

# if LJ == 12
        read_pairindex(fp4,GoPair);
# endif

        read_crd(fp1);

#if HETEROGO == YES
	printf("Input heterogeneous Go pair list>>");
	scanf("%s",corefile);
	//printf("%s\n",corefile);
	fp5=openfile(corefile,"rt");
	N=get_Ntot(fp5);
	for(i=0;i<N;i++)
	{
	fscanf(fp5,"%d%s%d%d",&j,type,&m,&n);
	j=j-1; /*become C index*/
	if(GoPair[j].res1!=m || GoPair[j].res2!=n){printf(" error matching minicore pairs, %d %d\n", m,n);exit(1);}
	GoPair[j].minicore=1;
	}	

# endif

/*atom mass*/
fprintf(fp2,"\n");

for(i=0;i<NATOM;i++){
        fprintf(fp2,"%s 50.0\n",native[i].name[0]);
        fprintf(fp2,"%s 50.0\n",native[i].name[1]);
}


        fprintf(fp2,"%s %.1f\n",xx,MG);
        fprintf(fp2,"%s %.1f\n",xy,MG);

fprintf(fp2,"\n\n");

        /*bond length*/

        /*between calphas*/     
        for(i=0;i<NATOM-1;i++)
        {
        fprintf(fp2,"%s %s    120.00      %lf\n",native[i].name[0],native[i+1].name[0],distance(native[i].coor[0],native[i+1].coor[0]));
        }

        /*between calpha-cbeta only for non-GLY residues*/        
        for(i=0;i<NATOM;i++)
        {       
        /*if(i!=4 && i!=14 && i!= 33)    */
	if (strcmp(native[i].type,"GLY")!=0)
        {
        fprintf(fp2,"%s %s    120.00      %lf\n",native[i].name[0],native[i].name[1],distance(native[i].coor[0],native[i].coor[1]));
        }
        }

	fprintf(fp2,"%s %s    5000.0      %.5lf\n",xx,xy,2.0*RG); //crowders
//disulfide bond length
	for(i=0;i<Nds;i++){
		aa = dsA[i]; bb = dsB[i];
		fprintf(fp2,"%s %s    120.00      %lf\n",native[aa].name[1],native[bb].name[1],distance(native[aa].coor[1],native[bb].coor[1]));
	}
        fprintf(fp2,"\n");

/*atom angle*/
        /*angle on calpha chain, ca-ca-ca*/       

for(i=0;i<NRES-2;i++)
        {
        fprintf(fp2,"%s %s %s     24.00    %.3lf\n",native[i].name[0],native[i+1].name[0],native[i+2].name[0],calan(native[i].coor[0],native[i+1].coor[0],native[i+2].coor[0]));
        }

        /*angle on calpha cbeta, cb-ca-ca for non-gly */
        /*first one*/   

if(strcmp(native[0].type,"GLY")!=0) {
        fprintf(fp2,"%s %s %s     24.00    %.3lf\n",native[0].name[1],native[0].name[0],native[1].name[0],calan(native[0].coor[1],native[0].coor[0],native[1].coor[0]));
	}

        /*middle one*/
for(i=1;i<NATOM-1;i++){
if(strcmp(native[i].type,"GLY")!=0) {
        fprintf(fp2,"%s %s %s     24.00    %.3lf\n",native[i].name[1],native[i].name[0],native[i-1].name[0],calan(native[i].coor[1],native[i].coor[0],native[i-1].coor[0]));
        
        fprintf(fp2,"%s %s %s     24.00    %.3lf\n",native[i].name[1],native[i].name[0],native[i+1].name[0],calan(native[i].coor[1],native[i].coor[0],native[i+1].coor[0]));
        } }

        /*last one*/

	i=NATOM-1;
     if(strcmp(native[i].type,"GLY")!=0) {
        fprintf(fp2,"%s %s %s     24.00    %.3lf\n",native[i].name[1],native[i].name[0],native[i-1].name[0],calan(native[i].coor[1],native[i].coor[0],native[i-1].coor[0]));
    }
//disulfide bond angle potential
for(i=0;i<Nds;i++){
	aa = dsA[i]; bb = dsB[i];
fprintf(fp2,"%s %s %s     24.00    %.3lf\n",native[aa].name[0],native[aa].name[1],native[bb].name[1],calan(native[aa].coor[0],native[aa].coor[1],native[bb].coor[1]));
fprintf(fp2,"%s %s %s     24.00    %.3lf\n",native[aa].name[1],native[bb].name[1],native[bb].name[0],calan(native[aa].coor[1],native[bb].coor[1],native[bb].coor[0]));
}
        fprintf(fp2,"\n");

/*dihedral*/
 /*dihedral only on calpha*/

//printf("%s %s\n",native[1].type,maparray[100].type1);

/*if(strcmp(native[1].type,maparray[100].type1)==0)
{ printf("Hello\n");}*/



for(i=0;i<NATOM-3;i=i+1)
{

 	caldi(&di,native[i].coor[0],native[i+1].coor[0],native[i+2].coor[0],native[i+3].coor[0]);	
	/*** di*3 is important*/


 /* Antonios changed here, sequence dependence in the dihedral's epilon */



	
 for (j=0;j<210;j++)   //big for loop
{
	if (((strcmp(native[i].type,maparray[j].type1)==0) && (strcmp(native[i+1].type,maparray[j].type2)==0)) || ((strcmp(native[i].type,maparray[j].type2)==0) && (strcmp(native[i+1].type,maparray[j].type1)==0))) {
//	printf("%s\t%s\t%f\n",native[i].type,native[i+1].type,maparray[j].eps);
	epsilon1 = maparray[j].eps;
	}

if (((strcmp(native[i].type,maparray[j].type1)==0) && (strcmp(native[i+2].type,maparray[j].type2)==0)) || ((strcmp(native[i].type,maparray[j].type2)==0) && (strcmp(native[i+2].type,maparray[j].type1)==0))) {
//        printf("%s%s\n",native[i].type,native[i+2].type);
	epsilon2 = maparray[j].eps;
	}


if (((strcmp(native[i].type,maparray[j].type1)==0) && (strcmp(native[i+3].type,maparray[j].type2)==0)) || ((strcmp(native[i].type,maparray[j].type2)==0) && (strcmp(native[i+3].type,maparray[j].type1)==0))) {
//        printf("%s%s\n",native[i].type,native[i+3].type);
	epsilon3 = maparray[j].eps;
        }

if (((strcmp(native[i+1].type,maparray[j].type1)==0) && (strcmp(native[i+2].type,maparray[j].type2)==0)) || ((strcmp(native[i+1].type,maparray[j].type2)==0) && (strcmp(native[i+2].type,maparray[j].type1)==0))) {
//        printf("%s%s\n",native[i+1].type,native[i+2].type);
	epsilon4 = maparray[j].eps;
        }

if (((strcmp(native[i+1].type,maparray[j].type1)==0) && (strcmp(native[i+3].type,maparray[j].type2)==0)) || ((strcmp(native[i+1].type,maparray[j].type2)==0) && (strcmp(native[i+3].type,maparray[j].type1)==0))) {
  //      printf("%s%s\n",native[i+1].type,native[i+3].type);
	epsilon5 = maparray[j].eps;
        }


if (((strcmp(native[i+2].type,maparray[j].type1)==0) && (strcmp(native[i+3].type,maparray[j].type2)==0)) || ((strcmp(native[i+2].type,maparray[j].type2)==0) && (strcmp(native[i+3].type,maparray[j].type1)==0))) {
//        printf("%s%s\n",native[i+2].type,native[i+3].type);
	epsilon6 = maparray[j].eps;
        }


av_epsilon = (epsilon1 + epsilon2 + epsilon3 + epsilon4 + epsilon5 + epsilon6)/6.0;

	} // end of for loop

// End of Antonios

//printf("%f\n",av_epsilon);
	 
	/*thirumalai dihedral*/
/*	di1= (int)di % 360;
	di1+= (di>=0 ? 0 : 359);
	modulo = (di1< 180 ? 90 : -90);
	printf("%d di %lf, di1 %d modulo %lf\n",i, di,di1, modulo);
	modulo = modulo+di;
       	fprintf(fp2,"%2s-%2s-%2s-%2s%4d%15.2lf%15.2lf%15.2lf\n", native[i].name[0],native[i+1].name[0],native[i+2].name[0],native[i+3].name[0],1,0.5, modulo,-1.0);
       	fprintf(fp2,"%2s-%2s-%2s-%2s%4d%15.2lf%15.2lf%15.2lf\n", native[i].name[0],native[i+1].name[0],native[i+2].name[0],native[i+3].name[0],1,1.0, di,1.0); 
*/
	 
if(UREA){
       	fprintf(fp2,"%2s-%2s-%2s-%2s%4d%15.2lf%15.2lf%15.2lf\n", native[i].name[0],native[i+1].name[0],native[i+2].name[0],native[i+3].name[0],1,1*av_epsilon, di*3,-3.0);
       	fprintf(fp2,"%2s-%2s-%2s-%2s%4d%15.2lf%15.2lf%15.2lf\n", native[i].name[0],native[i+1].name[0],native[i+2].name[0],native[i+3].name[0],1,2*av_epsilon, di,1.0); 
       	} 
else{
	fprintf(fp2,"%2s-%2s-%2s-%2s%4d%15.2lf%15.2lf%15.2lf\n", native[i].name[0],native[i+1].name[0],native[i+2].name[0],native[i+3].name[0],1,1*GoEn, di*3,-3.0);
       	fprintf(fp2,"%2s-%2s-%2s-%2s%4d%15.2lf%15.2lf%15.2lf\n", native[i].name[0],native[i+1].name[0],native[i+2].name[0],native[i+3].name[0],1,2*GoEn, di,1.0); 
	}
}

/*cb ca ca cb, assigned but the value is 0 */
for(i=0;i<NATOM-1;i=i+1)
{
if(strcmp(native[i].type,"GLY")!=0 && strcmp(native[i+1].type,"GLY")!=0) {

        caldi(&di,native[i].coor[1],native[i].coor[0],native[i+1].coor[0],native[i+1].coor[1]);
        fprintf(fp2,"%2s-%2s-%2s-%2s%4d%15.2lf%15.2lf%15.2lf\n", native[i].name[1],native[i].name[0],native[i+1].name[0],native[i+1].name[1],1,0.0, di,1.0);

}}

/*cb ca ca ca, assigned but the value is 0*/
for(i=0;i<NATOM-2;i=i+1)
{
if(strcmp(native[i].type,"GLY")!=0) {
	
        caldi(&di,native[i].coor[1],native[i].coor[0],native[i+1].coor[0],native[i+2].coor[0]);
        fprintf(fp2,"%2s-%2s-%2s-%2s%4d%15.2lf%15.2lf%15.2lf\n", native[i].name[1],native[i].name[0],native[i+1].name[0],native[i+2].name[0],1,0.0, di,1.0);
} }

/*ca ca ca cb, assigned but the value is 0*/

for(i=0;i<NATOM-2;i=i+1)
{
if(strcmp(native[i+2].type,"GLY")!=0) {
        caldi(&di,native[i].coor[0],native[i+1].coor[0],native[i+2].coor[0],native[i+2].coor[1]);
        fprintf(fp2,"%2s-%2s-%2s-%2s%4d%15.2lf%15.2lf%15.2lf\n", native[i].name[0],native[i+1].name[0],native[i+2].name[0],native[i+2].name[1],1,0.0, di,1.0);
}}
// Disulfide bond dihedral
for(i=0;i<Nds;i++){
	aa = dsA[i]; bb = dsB[i];
caldi(&di,native[aa-1].coor[0],native[aa].coor[0],native[aa].coor[1],native[bb].coor[1]);
fprintf(fp2,"%2s-%2s-%2s-%2s%4d%15.2lf%15.2lf%15.2lf\n", native[aa-1].name[0],native[aa].name[0],native[aa].name[1],native[bb].name[1],1,0.0, di,1.0);
caldi(&di,native[aa+1].coor[0],native[aa].coor[0],native[aa].coor[1],native[bb].coor[1]);
fprintf(fp2,"%2s-%2s-%2s-%2s%4d%15.2lf%15.2lf%15.2lf\n", native[aa+1].name[0],native[aa].name[0],native[aa].name[1],native[bb].name[1],1,0.0, di,1.0);
caldi(&di,native[aa].coor[0],native[aa].coor[1],native[bb].coor[1],native[bb].coor[0]);
fprintf(fp2,"%2s-%2s-%2s-%2s%4d%15.2lf%15.2lf%15.2lf\n", native[aa].name[0],native[aa].name[1],native[bb].name[1],native[bb].name[0],1,0.0, di,1.0);
caldi(&di,native[aa].coor[1],native[bb].coor[1],native[bb].coor[0],native[bb-1].coor[0]);
fprintf(fp2,"%2s-%2s-%2s-%2s%4d%15.2lf%15.2lf%15.2lf\n", native[aa].name[1],native[bb].name[1],native[bb].name[0],native[bb-1].name[0],1,0.0, di,1.0);
caldi(&di,native[aa].coor[1],native[bb].coor[1],native[bb].coor[0],native[bb+1].coor[0]);
fprintf(fp2,"%2s-%2s-%2s-%2s%4d%15.2lf%15.2lf%15.2lf\n", native[aa].name[1],native[bb].name[1],native[bb].name[0],native[bb+1].name[0],1,0.0, di,1.0);
}
zero=0;

fprintf(fp2,"\nX -X -C -O          10.5         180.          2.\n\n");

#if LJ==6
fprintf(fp2,"  HW  OW%15.2lf%15.2lf\n\n",zero,zero);
#elif LJ==10

      for(i=0;i<MAXPAIR;i++)
      {

	#if VdW10_12 == YES
        if(strcmp(GoPair[i].type,"b")==0 )
	         {
		
            B=6*GoEn*(TEN(GoPair[i].dist));
            A=5*GoEn*TWV(GoPair[i].dist);
	    dummy=GoEn;
		#if HETEROGO == YES
		  if(GoPair[i].minicore==1){
		    B=6*HGoEn*(TEN(GoPair[i].dist));
		    A=5*HGoEn*TWV(GoPair[i].dist);
		    dummy=HGoEn;
		  }
		#endif

            checkatomtype(c1,&GoPair[i].t1);
            checkatomtype(c2,&GoPair[i].t2);


		#if HETEROGO == YES
         // printf("%d  %s  %s  %lf minicore=%d\n",i+1,c1,c2,dummy,GoPair[i].minicore);
		#endif
          fprintf(fp2,"  %s  %s%20.2lf%20.2lf\n",c1,c2,A,B);
        }
	#endif	
 
	#if HB10_12 == YES 

        if(strcmp(GoPair[i].type,"h")==0 )
         {

            B=6*HB*(TEN(GoPair[i].dist));
            A=5*HB*TWV(GoPair[i].dist);
		dummy=HB;
		#if HETEROGO == YES
		  if(GoPair[i].minicore==1){
		    B=6*HHB*(TEN(GoPair[i].dist));
		    A=5*HHB*TWV(GoPair[i].dist);
			dummy=HHB;
		  }
		#endif

	    checkatomtype(c1,&GoPair[i].t1);
	    checkatomtype(c2,&GoPair[i].t2);
	   	   

		#if HETEROGO == YES
          //printf("%d  %s  %s %lf minicore=%d\n",i+1,c1,c2,dummy,GoPair[i].minicore);;
		#endif
          fprintf(fp2,"  %s  %s%20.2lf%20.2lf\n",c1,c2,A,B);
        }
	#endif
   }


#endif

fprintf(fp2,"N   NA\nC   C*\n\n");
fprintf(fp2,"MOD4      RE\n");

l=1;
for(i=0;i<NATOM;i++)
{
fprintf(fp2,"  %s          %d.000  0.6000             OPLS\n",native[i].name[0],l);
l++;
}

for(i=0;i<NATOM;i++)
{
fprintf(fp2,"  %s          %d.000  0.6000             OPLS\n",native[i].name[1],l);
l++;
}

/*crowding*/
fprintf(fp2,"  %s          0.0000  0.0000             OPLS\n",xx);
fprintf(fp2,"  %s          1.0000  0.0000             OPLS\n",xy);

fprintf(fp2,"\nEND\nEND");

fclose(fp2);

}

void checkatomtype(char c[4],int *t)
{
int i,j;

    for(i=0;i<NRES;i++)
	for(j=0;j<native[i].num;j++)	
	{
	if(native[i].typeindex[j]==(*t+1))
	 	{
		strcpy(c,native[i].name[j]);		
		}
	}
}

void read_atomindex(void)	
{
int i,j;
int N;
int count=0;
int typecount=0;
float xx;
/*accumalate atomic index*/
        N=get_Ntot(fp3);
        if(N!=NATOM){printf("wrong Natom %d\n",N);}
        for(i=0;i<NATOM;i++)
        {
                fscanf(fp3,"%d%s%f",&j,native[i].type,&xx);
                if(strcmp(native[i].type,"GLY")!=0) 
                {native[i].num=2;}
                if(strcmp(native[i].type,"GLY")==0) 
                {native[i].num=1;}

                 for(j=0;j<native[i].num;j++)
                 {
                 count++; /*fortran index*/
                 native[i].atomindex[j]=count;
                /*printf("%d %d %s %d\n",i,j,native[i].type,native[i].atomindex[j]);*/
                 }

                if(strcmp(native[i].type,"GLY")==0) 
                {
                 native[i].atomindex[1]=native[i].atomindex[0];
                } 
        }       
                if(count!=MAXLENGTH){printf("wrong with maxlength %d\n",count);exit(1);} 
                /*printf("%d %s\n",i,native[i].type);*/


        for(i=0;i<NATOM;i++)
        {
                         for(j=0;j<native[i].num;j++)
                        {
                         {
                         typecount++; /*f index*/
                        native[i].typeindex[j]=typecount;
                       // printf("hello %d %d %s %d %d\n",i,j,native[i].name[j],native[i].typeindex[j],native[i].atomindex[j]);
                         
                         }
                	}

                if(strcmp(native[i].type,"GLY")==0) 
                {
                 native[i].typeindex[1]=native[i].typeindex[0];
                } 


        }       
                if(count!=MAXLENGTH){printf("wrong with maxlength %d\n",count);exit(1);} 
                /*printf("%d %s\n",i,native[i].type);*/



}

void read_crd(FILE *fp)	
{

        int i,j,num,count;
	double x,y,z;
	char word[10];
        /****Read in native.crd****/

	fscanf(fp,"%s",word);
	fscanf(fp,"%d",&Ntot);
       
	count = 0;
        for(j=0;j<NATOM;j++)
	{
	  num=native[j].num;
	  for(i=0;i<num;i++)	
          {
		  fscanf(fp,"%lf", &native[j].coor[i].x);
		  fscanf(fp,"%lf", &native[j].coor[i].y);
		  fscanf(fp,"%lf", &native[j].coor[i].z);
		  count++;
	  }
	        if(strcmp(native[j].type,"GLY")==0) /*if gly for simplicity*/ 
		{
		native[j].coor[num].x = native[j].coor[num-1].x ;
                native[j].coor[num].y = native[j].coor[num-1].y ;
                native[j].coor[num].z = native[j].coor[num-1].z ;
			}
		
	}

	if(count!=Ntot){printf("check read_crd Ntot not match count=%d\n",count);exit(1);}

 rewind(fp);
       fscanf(fp,"%s",word);
       fscanf(fp,"%d",&num);

        for(j=0;j<MAXLENGTH;j++)
        {
                  fscanf(fp,"%lf", &polymer1[j].x);
                  fscanf(fp,"%lf", &polymer1[j].y);
                  fscanf(fp,"%lf", &polymer1[j].z);
        }

#if LJ == 12
          for(j=0;j<MAXPAIR;j++)
          {
         GoPair[j].dist=distance(polymer1[GoPair[j].a1], polymer1[GoPair[j].a2]);
	//printf("dist %d %lf \n",j+1,GoPair[j].dist);
          }
#endif
}



double distance(struct config poly1, struct config poly2)
  {
double dist;

dist = ( (poly1.x-poly2.x)*(poly1.x-poly2.x)+(poly1.y-poly2.y)*(poly1.y-poly2.y) + (poly1.z-poly2.z)*(poly1.z-poly2.z) );


return sqrt(dist);
  }


double calan(struct config A1, struct config A2, struct config A3)
{

#define vector(a,b,c) ( (c).x=(b).x-(a).x, (c).y =(b).y-(a).y, (c).z=(b).z-(a).z )
#define inner_product(a,b,c) (c=(a).x*(b).x+(a).y*(b).y+(a).z*(b).z)

int i,j,k;
static int count=0;
double prod=0;
double param=0;
double angle=0;

/* j,a1; j+1,a2; j+2,a3;*/

 prod=0;
  
    vector(A2,A1,vec1), vector(A2,A3,vec2);
   
     /*printf("%lf, %lf, %lf\t %lf, %lf, %lf\t",vec1.x,vec1.y,vec1.z,vec2.x,vec2.y,vec2.z);*/

    inner_product(vec1, vec2, prod);  

    prod=prod/(distance(A1,A2)*distance(A2,A3));


    /*printf("%lf\n",acos(prod)*CONVERT);*/

                                                 
    if (prod <=1 && prod >= -1){
      angle=acos(prod)*CONVERT;
	return angle;
    }
    else{printf("acos is out of range. %d th prod is %lf\n",count,(prod));exit(1);return 0;}

}



void read_pairindex(FILE *fp,struct list Pair[MAXPAIR])
{
int i,j,m,n;
int N;
int num;
char type[4];
int count;
count=0;
        N=get_Ntot(fp);
        printf("read_pairindex %d\n",N);
        for(i=0;i<N;i++)
        {
#if HETEROGO == YES
		Pair[i].minicore=0;
#endif
                fscanf(fp,"%d%s%d%d",&j,type,&m,&n);
                strcpy(Pair[i].type,type);
                /*printf("%s\n",Pair[i].type);*/
		Pair[i].res1=m;
		Pair[i].res2=n;
	
                if(strcmp(type,"b")==0)
                {
                 /*cb-cb*/
                 num=1;

                 /*this is fortran atomic index*/
                 Pair[i].index=j;
                 /*fortran atomic index to C atomic index*/
                 Pair[i].a1=native[m].atomindex[num]-1;
                 Pair[i].a2=native[n].atomindex[num]-1;

     
//printf("%d %s %d %d %d %d \n",i+1,type,Pair[i].a1,Pair[i].a2,m,n); 


                 /*fortran type index to C type index*/
                 Pair[i].t1=native[m].typeindex[num]-1;
                 Pair[i].t2=native[n].typeindex[num]-1;
/* printf("%d %s %d %d \n",i+1,type,GoPair[i].t1,GoPair[i].t2); */
                        }
                 else if(strcmp(type,"h")==0)
                {
                /*H-bond*/
                 Pair[i].index=j;
                 /*fortran atomic index to C atomic index*/
                 Pair[i].a1=native[m].atomindex[0]-1;
                 Pair[i].a2=native[n].atomindex[0]-1;

//printf("%d %s %d %d %d %d\n",i+1,type,Pair[i].a1,Pair[i].a2,m,n); 

                 /*fortran type index to C type index*/
                 Pair[i].t1=native[m].typeindex[0]-1;
                 Pair[i].t2=native[n].typeindex[0]-1;
                 /*printf("%d %s %d %d \n",i+1,type,Pair[i].t1,Pair[i].t2); */
                        
                } 
                else{printf("non type matched\n");exit(1);}
        
        }

}
 
FILE *openfile(char *filename, char *mode)
{
        FILE *fp;
 
        if((fp=fopen(filename, mode))==NULL){
        fprintf(stderr,"\nCan't open file %s.\n", filename);
        exit(1);}
        return fp;
}


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


void caldi(double *di,struct config A1,struct config A2,struct config A3,struct config A4)
{
#define THRESANGLE 360
#define PI 180
#define END -1
#define AMBERIN 1
	 
#define unitvec(a,b)    ((a).x/=b,(a).y/=b,(a).z/=b )
#define VECTOR(a,b,c) ( (c).x=(a).x-(b).x, (c).y =(a).y-(b).y, (c).z=(a).z-(b).z )
#define INNER_product(a,b) ((a).x*(b).x+(a).y*(b).y+(a).z*(b).z)
#define cross_product(a,b,c) (c.x=(a).y*(b).z-(b).y*(a).z, c.y=-(a).x*(b).z+(b).x*(a).z, (c).z=(a).x*(b).y-(b).x*(a).y )
  typedef struct {
  double x;
  double y;
  double z;
  } con;

	con vec, vec1, vec2, nvec,nvec2;
	double theta1, theta2, theta12, diangle;
	double Dist[4]={0};
	int N;

static int count =  0;
double d;
double dist1, dist2;
int j,k;
double prod=0;
double prod1=0;
double param=0;
double angle=0;

	/*calc_dist(poly, Dist);*/

	Dist[1]=distance(A1,A2);
	Dist[2]=distance(A2,A3);
	Dist[3]=distance(A3,A4);
/*
printf("%lf %lf %lf\n",Dist[1],Dist[2],Dist[3]);
printf("%lf %lf %lf\n",A1.x,A1.y,A1.z);*/

	    VECTOR(A1,A2,vec1); VECTOR(A4,A3,vec2); VECTOR(A2,A3,vec);

    		 /*printf("%lf, %lf, %lf\t %lf, %lf, %lf\t,%lf, %lf, %lf\n ",
			vec1.x,vec1.y,vec1.z,vec2.x,vec2.y,vec2.z,vec.x,vec.y,vec.z);*/
	    unitvec(vec1,Dist[1]);
	    unitvec(vec,Dist[2]);
	    unitvec(vec2,Dist[3]);

	    theta1=acos(INNER_product(vec,vec1));
	    theta2=acos(INNER_product(vec,vec2));

/*	    	printf("theta1 =%lf, theta2=%lf\n", theta1, theta2);*/

	    diangle=(INNER_product(vec1,vec2)-cos(theta1)*cos(theta2) )/(sin(theta1)*sin(theta2)); /*absolute value*/

	    /*printf("%lf \n",diangle);*/
	    cross_product(vec1,vec2,nvec);
	    prod1=INNER_product(nvec,vec); /*sign*/


/*for amber parameter file format*/

		if (diangle<=1 && diangle>= -1){
		      diangle=acos(diangle)*CONVERT;
				if(prod1<0) 
					{
					//printf("%.2lf\n",diangle-PI);
					/*DI[j]=diangle-PI;*/
					*di=diangle-PI;
					}
				else
					{
					//printf("%.2lf\n",-diangle-PI);
					/*DI[j]=-diangle-PI;*/
					*di=-diangle-PI;
					}
			

		/*	*di=diangle;*/

		    }

		else{printf("acos is out of range. %d th prod is %lf\n",count,diangle);exit(1);}


}

void read_mapfile(void)

{

int i,j;

for (i=0;i<210;i++)
	{
	fscanf(fp6,"%s %s %f",maparray[i].type1,maparray[i].type2,&maparray[i].eps);

/* Check of BTmap is read correctly OK!
	printf("%s   %s    %f   \n",maparray[i].type1,maparray[i].type2,maparray[i].eps); */

	}


}

//reading the list of disulfide bonds
void read_disulfile(void)
{
	int i, A, B;

	fscanf(fp7,"%d",&Nds);
	printf("%d\n",Nds);
	for (i=0;i<Nds;i++)
	{
		fscanf(fp7,"%d %d",&A,&B);
		dsA[i]=A;
		dsB[i]=B;
	}

}
