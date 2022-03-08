/*********************************************
*triple.c
*This program is to evaluate the triple scaler product from the 
*output coordinates from amber.

(AxB)*C=(BxC)*A=(AxC)*B=det|A B C|

* The command line : a.out mdcrd Ndat input.crd *7.5.02 by Margaret Cheung
**********************************************/

#include "IOab.h"
#define YES 1
#define NO 0
#define vector(a,b,c) ( (c).x=(b).x-(a).x, (c).y =(b).y-(a).y, (c).z=(b).z-(a).z )
#define inner_product(a,b,c) (c=(a).x*(b).x+(a).y*(b).y+(a).z*(b).z)
#define distance1(a,b,c) (c= ((a).x-(b).x)* ((a).x-(b).x) +   ((a).y-(b).y)* ((a).y-(b).y) +  ((a).z-(b).z)* ((a).z-(b).z) )


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

struct config vec1, vec2, vec3;

void read_crd(FILE *fp);
void read_atomindex(void);
void read_pairindex(void);
double calan(struct config A1, struct config A2, struct config A3);
double distance(struct config poly1, struct config poly2);

char inputfile[20], nativefile[20],atomfile[20],pairfile[20];

int Ntot;
double triple;
FILE *fp1, *fp2, *fp3, *fp4, *fp5,*fp6;
int HBNtot=0;
main()
{
double kchi;
int i,j,k;
int Ndat;
char word[20];
int ca,cb,c,n; 
int count=0;
printf("a.out: native.crd atomfile pairfile kchi>>\n");
scanf("%s%s%s%lf",nativefile,atomfile,pairfile,&kchi);
printf("native.crd: %s, atomfile: %s, pairfile %s, kchi %lf\n",nativefile,atomfile,pairfile,kchi); 



fp1=openfile(nativefile,"rt");
fp2=openfile(atomfile,"rt");
fp3=openfile(pairfile,"rt");
fp6=openfile("triple.inp","wt");
       /*read index crd*/
        read_atomindex();

        /*read_pairindex();*/
        read_crd(fp1);

	count=0;
	/*first and last dont have chiral problem*/
	for(i=1;i<NATOM-1;i++)
	{
		if(strcmp(native[i].type,"GLY")!=0 ) {
		count++;	
		}
	}

	fprintf(fp6,"%d %lf\n",count,kchi);

	for(i=1;i<NATOM-1;i++)
	{
	/*ca, cb, n, c*/	
	if(strcmp(native[i].type,"GLY")!=0  ){

	vector(native[i].coor[0],native[i].coor[1],vec1);
	vector(native[i].coor[0],native[i-1].coor[0],vec2);
	vector(native[i].coor[0],native[i+1].coor[0],vec3);
	triple= (vec1.x)*(vec2.y)*(vec3.z)
	       +(vec1.z)*(vec2.x)*(vec3.y)
	       +(vec1.y)*(vec2.z)*(vec3.x)
	       -(vec1.z)*(vec2.y)*(vec3.x)
	       -(vec1.y)*(vec2.x)*(vec3.z)
	       -(vec1.x)*(vec2.z)*(vec3.y);

	fprintf(fp6,"%d %d %d %d %lf\n",native[i].typeindex[0],native[i].typeindex[1],native[i-1].typeindex[0],native[i+1].typeindex[0], triple);
	
	}}
}

double calan(struct config A1, struct config A2, struct config A3)
{

#define vector(a,b,c) ( (c).x=(b).x-(a).x, (c).y =(b).y-(a).y, (c).z=(b).z-(a).z )
#define inner_product(a,b,c) (c=(a).x*(b).x+(a).y*(b).y+(a).z*(b).z)
#define CONVERT 180/3.1415927

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
double dummy;

/*accumalate atomic index*/
        N=get_Ntot(fp2);
        if(N!=NATOM){printf("wrong Natom %d\n",N);}
        for(i=0;i<NATOM;i++)
        {
                fscanf(fp2,"%d%s%lf",&j,native[i].type,&dummy);
                native[i].num=2;

                if(strcmp(native[i].type,"GLY")==0 ) 
                {native[i].num=1;}

                 for(j=0;j<native[i].num;j++)
                 {
                 count++; /*fortran index*/
                 native[i].atomindex[j]=count;
                printf("%d %d %s %d\n",i,j,native[i].type,native[i].atomindex[j]);
                 }

                if(strcmp(native[i].type,"GLY")==0  ) 
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

                if(strcmp(native[i].type,"GLY")==0 ) 
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
        N=get_Ntot(fp3);
        if(N!=MAXPAIR){printf("wrong maxpair %d",N);exit(1);}
        for(i=0;i<N;i++)
        {
                fscanf(fp3,"%d%s%d%d",&j,type,&m,&n);
                strcpy(GoPair[i].type,type);
	   	GoPair[i].alpha1=m;
		GoPair[i].alpha2=n;
                /*printf("%s\n",GoPair[i].type);*/

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

printf("%d %s %d %d \n",i+1,type,GoPair[i].a1,GoPair[i].a2); 

                 /*fortran type index to C type index*/
                 GoPair[i].t1=native[m].typeindex[0]-1;
                 GoPair[i].t2=native[n].typeindex[0]-1;
                 /*printf("%d %s %d %d \n",i+1,type,GoPair[i].t1,GoPair[i].t2); */
                        
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
                if(strcmp(native[j].type,"GLY")==0 ) /*if gly for simplicity*/ 
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


     /*     for(j=0;j<MAXPAIR;j++)
          {
         GoPair[j].dist=distance(polymer1[GoPair[j].a1], polymer1[GoPair[j].a2]);
          }*/

}





double distance(struct config poly1, struct config poly2)
  {
double dist;
dist = 0;
dist = ( (poly1.x-poly2.x)*(poly1.x-poly2.x)+(poly1.y-poly2.y)*(poly1.y-poly2.y) + (poly1.z-poly2.z)*(poly1.z-poly2.z) );

return sqrt(dist);
  }











