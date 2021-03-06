/*********************************************
*countHB.c
*This program is to evaluate the angles from the *output coordinates from amber.
* The command line : a.out mdcrd Ndat input.crd
*7.14.98 by Margaret Cheung
**********************************************/

#include "IOab.1.h"
#define YES 1
#define NO 0
#define PSEUDO YES 
#define NONNATIVE NO
#define THRESANGLE 360
#define CONVERT 180/3.1415927
#define PI 180 
#define PIRAD 3.1415927
#define pt000 0.0001
#define END -1
#define vector(a,b,c) ( (c).x=(b).x-(a).x, (c).y =(b).y-(a).y, (c).z=(b).z-(a).z )
#define inner_product(a,b,c) (c=(a).x*(b).x+(a).y*(b).y+(a).z*(b).z)
#define distance1(a,b,c) (c= ((a).x-(b).x)* ((a).x-(b).x) +   ((a).y-(b).y)* ((a).y-(b).y) +  ((a).z-(b).z)* ((a).z-(b).z) )
#define GoEn 0.6
#define SIX(a) a*a*a*a*a*a
#define TWV(a) SIX(a)*SIX(a)
#define COSPHIA 0.893373
#define Bfactor 1
FILE *openfile(char *filename, char *mode);
int nextline(FILE *fp);
int get_Ntot(FILE *fp);
                                                                                                                                                             
struct bond{
double HB;
int hi;
int hi1;
int hi2;
int hi3;
double GoDist;
};

struct list{
  int index;
  int a1;
  int a2;
  int t1;
  int t2;
  int alpha1;
  int alpha2;
  double dist;
  char type[4];
 };

struct list GoPair[MAXPAIR];
struct config vec1, vec2;
struct bond pseudoHB[MAXHBPAIR+1];

void read_crd(FILE *fp);
void read_atomindex(void);
void read_pairindex(void);
double calan(struct config A1, struct config A2, struct config A3);
void caldi(double *di,struct config A1,struct config A2,struct config A3,struct config A4);
double distance(struct config poly1, struct config poly2);
void countHB(void);


double angle[MAXLENGTH];
int map[MAXLENGTH][MAXLENGTH];
char inputfile[20], nativefile[20],atomfile[20],pairfile[20];

int HBNtot;
int Ntot;
int Ndat;
FILE *fp1, *fp2, *fp3, *fp4, *fp5,*fp6,*fp7;

main()
{

int i,j,k;
double prod=0;
double param=0;
double angle=0;
double amp;
double di=0.0;
char word[20];
double phi,theta;
int xi,xi1,xi2,xi3; 
int count=0;
/*if(argc<4){printf("Usage: a.out mdcrd Ndat input.crd\n"); exit(1);}*/

printf("a.out: native.crd atomfile pairfile Ndat md.crd>>\n");
scanf("%s%s%s%d%s",nativefile,atomfile,pairfile,&Ndat,inputfile);
printf("native.crd: %s, atomfile: %s, pairfile %s\n",nativefile,atomfile,pairfile); 

fp1=openfile(nativefile,"rt");
fp2=openfile(atomfile,"rt");
fp3=openfile(pairfile,"rt");
fp4=openfile(inputfile,"rt");
fp5=openfile("countHB","wt");

/*
#if PSEUDO == YES 
fp6=openfile("pseudo.inp","wt");
printf("input amp=(1)\n");
scanf("%lf",&amp);
	amp=Bfactor;
	fprintf(fp6,"%d %lf\n",MAXHBPAIR,amp);
#else 
fp6=openfile("hbond.inp","wt");
	fprintf(fp6,"%d\n",MAXHBPAIR);
#endif
*/

fp6=openfile("pseudo.inp","wt");
amp=Bfactor;

        read_atomindex();
        read_pairindex();
        read_crd(fp1);


#if NONNATIVE == NO

fprintf(fp6,"%d %lf\n",MAXHBPAIR,amp);
for(j=0;j<MAXPAIR;j++){

/*xi2<xi1*/

if( strcmp(GoPair[j].type,"h")==0) {
	count++;
	/*printf("%d %d\n",j,count);*/
	xi2= GoPair[j].alpha1;
	xi1= GoPair[j].alpha2;
/*if no gly*/
	/*if(strcmp(native[GoPair[j].alpha1].type,"a")==0 || strcmp(native[GoPair[j].alpha2].type,"a")==0 ) {printf("contact with gly %d\n",j+1);continue;}*/
	if(xi2>=xi1){printf("xi2 has to be smaller than xi1:pair %d\n",GoPair[j].index); exit(1); }


	/*vector is referenced by C-terminus ,pointed to N*/
	/*r1=xi1-xi,r2=xi2-xi1,r3=xi3-xi2*/
	/*costheta=-r2*r3/|r2||r3|*/
	/*cosphi=-r1*r2/|r1||r2| */
	/* 0<theta,phi<PI*/

	/* we first dont consider the terminal HB, will check here later*/ 
	xi3=xi2-1;
	xi=xi1+1;

printf("before %d %d %d %d\n",xi,xi1,xi2,xi3);	
	if(xi1< (NATOM-1) && xi2>0 )
	{
	xi3=xi2-1;
        xi=xi1+1;

        }
        else if(xi1== (NATOM-1) && xi2>0 )
         {
         
        xi3=xi2-1;
        xi=xi1-1;

         }
        else if (xi2==0 && xi1< (NATOM-1))
        {

        xi3=xi2+1;
        xi=xi1+1;

        }       

	else if (xi2==0 && xi1== (NATOM-1))
	{
	xi3=xi2+1;         xi=xi1-1;
	}
  
printf("after %d %d %d %d\n",xi,xi1,xi2,xi3);	
	/*need F index*/
	fprintf(fp6,"%d %d %d %d \n",native[xi].atomindex[0],native[xi1].atomindex[0],native[xi2].atomindex[0],native[xi3].atomindex[0]);
//	theta=calan(native[xi1].coor[0],native[xi2].coor[0],native[xi3].coor[0]);	
//	phi=calan(native[xi].coor[0],native[xi1].coor[0],native[xi2].coor[0]);

/*printf("phi %lf theta %lf \n",phi,theta);*/
	if(phi<0)
	{phi=-phi;}
	else if(phi>PIRAD)
	{phi=2*PIRAD-phi;}

	if(theta<0)
	{theta=-theta;}
	else if(theta>PIRAD)
	{theta=2*PIRAD-theta;}

	/*fprintf(fp6,"%lf %lf\n",theta,phi);	*/

      caldi(&di,native[xi3].coor[0],native[xi2].coor[0],native[xi1].coor[0],native[xi].coor[0]);	
	pseudoHB[count-1].HB=di;
	pseudoHB[count-1].hi=native[xi].atomindex[0];
	pseudoHB[count-1].hi1=native[xi1].atomindex[0];
	pseudoHB[count-1].hi2=native[xi2].atomindex[0];
	pseudoHB[count-1].hi3=native[xi3].atomindex[0];
	pseudoHB[count-1].GoDist=GoPair[j].dist;
	printf("pseudoHB %d %d %d %d\n",pseudoHB[count-1].hi,pseudoHB[count-1].hi1,pseudoHB[count-1].hi2,pseudoHB[count-1].hi3);
	}


	} /*end of j*/

/*good above*/

	countHB();
#endif

#if NONNATIVE == YES 

	count=0;
	for(i=0;i<NATOM;i++){
	for(j=i+4;j<NATOM;j++){
	count++;	
	}}
	fprintf(fp6,"%d %lf\n",count,amp);


	for(i=0;i<NATOM;i++){
	for(j=i+4;j<NATOM;j++){

        xi2= i;
        xi1= j;
        if(xi2>=xi1){printf("xi2 has to be smaller than xi1:pair %d\n",GoPair[j].index); exit(1); }


        /*vector is referenced by C-terminus ,pointed to N*/

        /* we first dont consider the terminal HB, will check here later*/
        xi3=xi2-1;
        xi=xi1+1;

        if(xi1< (NATOM-1) && xi2>0 )
        {
        xi3=xi2-1;
        xi=xi1+1;

        }
        else if(xi1== (NATOM-1) && xi2>0 )
         {

        xi3=xi2-1;
        xi=xi1-1;

         }
        else if (xi2==0 && xi1< (NATOM-1))
        {

        xi3=xi2+1;
        xi=xi1+1;

        }

        else if (xi2==0 && xi1== (NATOM-1))
        {
        xi3=xi2+1;         xi=xi1-1;
        }
 
        /*need F index*/
        fprintf(fp6,"%d %d %d %d \n",native[xi].atomindex[0],native[xi1].atomindex[0],native[xi2].atomindex[0],native[xi3].atomindex[0]);
	}
	}
#endif



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
      /*angle=acos(prod)*CONVERT;*/
      angle=acos(prod);
        return angle;
    }
    else{printf("acos is out of range. %d th prod is %lf\n",count,(prod));exit(1);return 0;}

}


void read_atomindex(void)       
{
int i,j;
int N;
int count=0;
int typecount=0;

/*accumalate atomic index*/
        N=get_Ntot(fp2);
        if(N!=NATOM){printf("wrong Natom %d\n",N);}
        for(i=0;i<NATOM;i++)
        {
                fscanf(fp2,"%d%s%lf",&j,native[i].type,&native[i].side);
                if(strcmp(native[i].type,"GLY")!=0) 
                {native[i].num=2;}
                if(strcmp(native[i].type,"GLY")==0) 
                {native[i].num=1;}

                 for(j=0;j<native[i].num;j++)
                 {
                 count++; /*fortran index*/
                 native[i].atomindex[j]=count;
                printf("%d %d %s %d\n",i,j,native[i].type,native[i].atomindex[j]);
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
                        printf("%d %d %s %d %d\n",i,j,native[i].name[j],native[i].typeindex[j],native[i].atomindex[j]);
                         
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

void read_pairindex(void)
{
int i,j,m,n;
int N;
int num;
char type[4];
int count;
count=0;
HBNtot=0;
        N=get_Ntot(fp3);
        if(N!=MAXPAIR){printf("wrong maxpair %d",N);exit(1);}
        for(i=0;i<N;i++)
        {
                fscanf(fp3,"%d%s%d%d",&j,type,&m,&n);
//		printf("%d %d\n", m, n);
                strcpy(GoPair[i].type,type);
	   	GoPair[i].alpha1=m;
		GoPair[i].alpha2=n;
                /*printf("%s\n",GoPair[i].type);*/
//printf("%d %d\n", GoPair[i].alpha1, GoPair[i].alpha2);
                if(strcmp(type,"b")==0)
                {
                 /*cb-cb*/
                 num=1;

                 /*this is fortran atomic index*/
                 GoPair[i].index=j;
                 /*fortran atomic index to C atomic index*/
                 GoPair[i].a1=native[m].atomindex[num]-1;
                 GoPair[i].a2=native[n].atomindex[num]-1;

     
printf("%d %s %d %d \n",i+1,type,GoPair[i].a1,GoPair[i].a2); 


                 /*fortran type index to C type index*/
                 GoPair[i].t1=native[m].typeindex[num]-1;
                 GoPair[i].t2=native[n].typeindex[num]-1;
/* printf("%d %s %d %d \n",i+1,type,GoPair[i].t1,GoPair[i].t2); */
                        }
                 else if(strcmp(type,"h")==0)
                {
                /*H-bond*/

	 	 HBNtot++;

                 GoPair[i].index=j;
                 /*fortran atomic index to C atomic index*/
                 GoPair[i].a1=native[m].atomindex[0]-1;
                 GoPair[i].a2=native[n].atomindex[0]-1;
printf("%d %s %d %d \n",i+1,type,native[m].atomindex[0],native[n].atomindex[0]);
/*printf("%d %s %d %d \n",i+1,type,GoPair[i].a1,GoPair[i].a2); */

                 /*fortran type index to C type index*/
                 GoPair[i].t1=native[m].typeindex[0]-1;
                 GoPair[i].t2=native[n].typeindex[0]-1;
                 printf("%d %s %d %d \n",i+1,type,GoPair[i].t1,GoPair[i].t2); 
                        
                } 
                else{printf("non type matched\n");exit(1);}
        }

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
                if(strcmp(native[j].type,"a")==0) /*if gly for simplicity*/ 
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


          for(j=0;j<MAXPAIR;j++)
          {
         GoPair[j].dist=distance(polymer1[GoPair[j].a1], polymer1[GoPair[j].a2]);
          }

}





double distance(struct config poly1, struct config poly2)
  {
double dist;
dist = 0;
dist = ( (poly1.x-poly2.x)*(poly1.x-poly2.x)+(poly1.y-poly2.y)*(poly1.y-poly2.y) + (poly1.z-poly2.z)*(poly1.z-poly2.z) );

return sqrt(dist);
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

FILE *openfile(char *filename, char *mode)
{
        FILE *fp;

        if((fp=fopen(filename, mode))==NULL){
        fprintf(stderr,"\nCan't open file %s.\n", filename);
        exit(1);}
        return fp;
}



void caldi(double *di,struct config A1,struct config A2,struct config A3,struct config A4)
{
#define THRESANGLE 360
	/*#define CONVERT 1*/
	/*#define PI 3.1415926536*/
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
		/*      diangle=acos(diangle)*CONVERT;
				if(prod1<0) 
					{
					**printf("%.2lf\n",diangle-PI);**
					**DI[j]=diangle-PI;**
					*di=diangle-PI;
					}
				else
					{
					**printf("%.2lf\n",-diangle-PI);**
					**DI[j]=-diangle-PI;**
					*di=-diangle-PI;
					}
		*/	

			*di=diangle;

		    }

		else{printf("acos is out of range. %d th prod is %lf\n",count,diangle);exit(1);}


}

void countHB(void)
{

int i,j;
int xi,xi1,xi2,xi3;
int countline;
int HB;
char word[20];
double di,dist,ener,dist6,dist12,angleHB,dummy;
fscanf(fp4,"%s",word);

countline = 0;


while(countline<Ndat){
  /*read in*/
	HB=0;
	   countline++;
	for(j=0;j<Ntot;j++){
	fscanf(fp4,"%lf", &polymer[j].x);
	fscanf(fp4,"%lf", &polymer[j].y);
	fscanf(fp4,"%lf", &polymer[j].z);

	/*printf("%lf %lf %lf\n",polymer[j].x,polymer[j].y,polymer[j].z);*/

	}
                /*strip out crd of crowding agent and x y z of the box size*/
/*	for(j=Ntot;j<MAXBOX+1;j++){
	fscanf(fp2,"%lf", &dummy);
	fscanf(fp2,"%lf", &dummy);
	fscanf(fp2,"%lf", &dummy);
	}
*/


	for(i=0;i<7;i++)
	{


	xi=pseudoHB[i].hi;
	xi1=pseudoHB[i].hi1;
	xi2=pseudoHB[i].hi2;
	xi3=pseudoHB[i].hi3;

      /*printf("new HB %d %d %d %d\n",xi,xi1,xi2,xi3);*/
      /*printf("new HB %lf %lf %lf\n",polymer[xi3].x,polymer[xi3].y,polymer[xi3].z);*/

      caldi(&di,polymer[xi3],polymer[xi2],polymer[xi1],polymer[xi]);
      dist=distance(polymer[xi2], polymer[xi1]);
      dist=pseudoHB[i].GoDist/dist;
      dist6=SIX(dist);
      dist12=dist6*dist6;
      ener = -2*GoEn*dist6+GoEn*dist12; 

/*      angleHB= (1+cos(di)) * (1-cos(di)) * (1-(cos(di)/COSPHIA));
      angleHB= angleHB*angleHB*Bfactor + 1;

      
      if(angleHB>0){
      angleHB=1/angleHB;}
      else{angleHB=0;}

      ener=ener*angleHB;
      ener=ener/GoEn;
      ener=fabs(ener);
      if(ener>0.15){HB++;}
*/

      /*printf("i %d,ener %lf, angle %lf\n",i,ener,angleHB);*/

	}/*end of i*/

	fprintf(fp5,"%d\n",HB);

	}/*end of while*/

}
