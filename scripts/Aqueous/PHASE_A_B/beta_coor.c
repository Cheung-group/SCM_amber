/**************************************************************************
beta_coor.c
	This program is to evaluate the xyz coordinate of beta atom in an 
	alpha-beta model 
	the command line is : order mdcrd Ndat native.crd ratio

5.11.01 by Margaret Cheung

make sure the C, N terminal atoms were gone.

****************************************************************************/

#include "IOab.h" 
#define INIT 0
#define ATOMINDEX 0 /*the first number of for the atom index in pdb*/
struct beta{
double x[3];
double mass;
char type[4];
int number;
};

struct RES{
char amino[4];
int num;
struct beta sidechain[15];
struct beta betamodel[4]; /*N:0,Ca:1,Cb:2,C:3*/
double dist;
};

struct RES CI2[NATOM];

void read_pdb(FILE *fp);
int amino_type(char AMIN[4]);
double findmass(char type[4]);
void amino_data(char TYPE[4],char AMIN[4],int *NUM,double *X,double *Y, double *Z,int *RES,int *ATO, int *atomin);
void COM(FILE *fp);

double distance(double x1[3], double x2[3]);

char inputfile[30], outfile[20], atomfile[20];
int Ntot;

FILE *fp1, *fp2, *fp3, *fp4,*fp5;

main(void)
{

int i,j,k;

/*****Read in from the command line and openfiles*****/

	printf("Usage: a.out file.pdb\n"); 

	scanf("%s",inputfile);
	fp1= fopen(inputfile,"r");
	sprintf(outfile,"%s.beta",inputfile);
	fp2= fopen(outfile,"wt");

	/****Read in pdb file****/

	Ntot=get_Ntot(fp1);

	read_pdb(fp1);
	printf("done read pdb\n");
        fclose(fp1);

	COM(fp2);
	fclose(fp2);

	sprintf(atomfile,"%s.atom",inputfile);
	fp5= fopen(atomfile,"wt");
	for(i=0;i<NATOM;i++)
	{
	printf ("%d %s %f\n",i,CI2[i].amino,CI2[i].dist);
	fprintf (fp5,"%d %s %f\n",i,CI2[i].amino,CI2[i].dist);
		}
}


void read_pdb(FILE *fp)
{
        int i,j,line,subline,init;
	
	char atom[5];
	char type[4];
	char amin[4];	
	int num,atomindex;
	double x,y,z;
	line=0;
	j=0;	
	init=0;

	/****Read in pdb****/
	while(line<Ntot)	
	{

	  	
	  fscanf(fp,"%s%d%s%s%*3c%d%lf%lf%lf%*26c\n",atom,&atomindex,type,amin,&num,&x,&y,&z);

	  amino_data(type,amin,&num,&x,&y,&z,&j,&init,&atomindex);
	  printf("%s %s %s %d %lf %lf %lf %d\n\n",atom,type,amin,num,x,y,z,atomindex);

         /*match the amino type to determine read how many lines*/
 
	subline=amino_type(amin); 
	printf("%s %d\n",amin,subline);
	
		for(i=1;i<subline;i++)	
		  {	
		fscanf(fp,"%s%d%s%s%*3c%d%lf%lf%lf%*26c\n",atom,&atomindex,type,amin,&num,&x,&y,&z);
	  printf("%s %s %s %d %lf %lf %lf\n\n",atom,type,amin,num,x,y,z);
		amino_data(type,amin,&subline,&x,&y,&z,&j,&i,&atomindex);
	         		
		}
	
	line=line+subline;
	j++;
	}
}

double findmass(char type[4])
{

	
	if( strcmp(type,"N")==0 || strcmp(type,"O")==0 || strcmp(type,"CA")==0 ||  strcmp(type,"C")==0){return 0;}
      else if( strcmp(type,"ND2")==0 || strcmp(type,"ND1")==0  || strcmp(type,"NZ")==0 || strcmp(type,"NE1")==0 ||  strcmp(type,"NE2")==0){return 14;}
      else if( strcmp(type,"NH1")==0 || strcmp(type,"NH2")==0 || strcmp(type,"NE")==0 ||  strcmp(type,"ND")==0){return 14;}
      else if( strcmp(type,"OG")==0 || strcmp(type,"OE1")==0 || strcmp(type,"OE2")==0 ||  strcmp(type,"OG1")==0){return 16;}
      else if( strcmp(type,"OD1")==0 || strcmp(type,"OD2")==0 || strcmp(type,"OH")==0){return 16;}
      else if( strcmp(type,"CB")==0 || strcmp(type,"CG")==0 || strcmp(type,"CE")==0 ||  strcmp(type,"CG2")==0){return 12;}
      else if( strcmp(type,"CD")==0 || strcmp(type,"CD1")==0 || strcmp(type,"CD2")==0 ||  strcmp(type,"CZ2")==0){return 12;}
      else if( strcmp(type,"CZ3")==0 || strcmp(type,"CH2")==0 || strcmp(type,"CG1")==0 ||  strcmp(type,"CZ")==0){return 12;}
      else if( strcmp(type,"CE3")==0 || strcmp(type,"CE2")==0 || strcmp(type,"CE1")==0 ){return 12;}
      else if( strcmp(type,"SD")==0 || strcmp(type,"SG")==0 ){return 32;}
      else if( strcmp(type,"OXT")==0 ){return 0;}
	else { printf("can not find mass of %s\n",type); exit(1);return 0;}	

}

void amino_data(char TYPE[4],char AMIN[4],int *NUM,double *X,double *Y, double *Z,int *RES,int *ATO,int *atomin)
{

      	strcpy(CI2[*RES].amino,AMIN);
      	CI2[*RES].num=*NUM;	
      	CI2[*RES].sidechain[*ATO].x[0]=*X;	
      	CI2[*RES].sidechain[*ATO].x[1]=*Y;	
      	CI2[*RES].sidechain[*ATO].x[2]=*Z;	
	CI2[*RES].sidechain[*ATO].number=*atomin-ATOMINDEX;
      	strcpy(CI2[*RES].sidechain[*ATO].type,TYPE);	
	
      	CI2[*RES].sidechain[*ATO].mass=findmass(TYPE);	

	printf("%s %d %lf %lf %lf %s %lf\n",CI2[*RES].amino,CI2[*RES].num,CI2[*RES].sidechain[*ATO].x[0],CI2[*RES].sidechain[*ATO].x[1],CI2[*RES].sidechain[*ATO].x[2],CI2[*RES].sidechain[*ATO].type,CI2[*RES].sidechain[*ATO].mass);
}

void COM(FILE *fp)
{

int i,j,k;
double X[3];
double M;
int count;
double dist;
count=0;
	/*print n-ca-cb-c*/
	for(i=0;i<NATOM;i++)
	{

		for(k=0;k<(CI2[i].num);k++)
		{
	
	        if( strcmp(CI2[i].sidechain[k].type,"CA")==0)	
			{
			for(j=0;j<3;j++) {
			fprintf(fp,"%lf ",CI2[i].sidechain[k].x[j]);
		 	CI2[i].betamodel[1].x[j]=CI2[i].sidechain[k].x[j];	
			printf("ca %lf ",CI2[i].betamodel[1].x[j]);
				}
		        fprintf(fp,"\n");
		        printf("\n");
			}
		}
	  	   
	/*print cb not for gly*/
	 if(strcmp(CI2[i].amino,"GLY")!=0){
	   for(j=0;j<3;j++)
		{	
		X[j]=0;
		M=0;
		/*printf("\n%d\t", CI2[i].num);*/
			for(k=0;k<(CI2[i].num);k++)
			{
			X[j]+=(CI2[i].sidechain[k].x[j])*(CI2[i].sidechain[k].mass);
			M+=CI2[i].sidechain[k].mass;
			/*printf("%lf %lf\t",CI2[i].sidechain[k].x[j],CI2[i].sidechain[k].mass);*/
			}
		X[j]=X[j]/M;
 		fprintf(fp,"%lf ",X[j]);
		CI2[i].betamodel[2].x[j]=X[j];	
 		printf("cb %lf ",CI2[i].betamodel[2].x[j]);
		}
		fprintf(fp,"\n");
		printf("\n");
	  
		} /*end of nogly*/

	 if(strcmp(CI2[i].amino,"GLY")==0)
	 { 
		 
		 for(j=0;j<3;j++){
		 CI2[i].betamodel[2].x[j]=CI2[i].betamodel[1].x[j];	
		 printf("%d %lf ",i,CI2[i].betamodel[2].x[j]);}
				}
		
		fprintf(fp,"\n");
		printf("\n");

	    }	
printf("\n\n");
	for(i=0;i<NATOM;i++){
	dist=distance(CI2[i].betamodel[2].x,CI2[i].betamodel[1].x);
	printf("dist %d %lf %lf %lf\n",i,CI2[i].betamodel[2].x[0],CI2[i].betamodel[1].x[0],dist);	
	CI2[i].dist=dist;
	}
}

int amino_type(char AMIN[4])
{

	if( strcmp(AMIN,"ASN")==0 ){return 8;}
	else if( strcmp(AMIN,"ASP")==0 ){return 8;}
	else if( strcmp(AMIN,"LEU")==0 ){return 8;}
	else if( strcmp(AMIN,"LYS")==0 ){return 9;}
	else if( strcmp(AMIN,"THR")==0 ){return 7;}
	else if( strcmp(AMIN,"GLU")==0 ){return 9;}
	else if( strcmp(AMIN,"TRP")==0 ){return 14;}
	else if( strcmp(AMIN,"PRO")==0 ){return 7;}
	else if( strcmp(AMIN,"VAL")==0 ){return 7;}
	else if( strcmp(AMIN,"GLY")==0 ){return 4;}
	else if( strcmp(AMIN,"SER")==0 ){return 6;}
	else if( strcmp(AMIN,"ALA")==0 ){return 5;}
	else if( strcmp(AMIN,"ILE")==0 ){return 8;}
	else if( strcmp(AMIN,"MET")==0 ){return 8;}
	else if( strcmp(AMIN,"TYR")==0 ){return 12;}
	else if( strcmp(AMIN,"ARG")==0 ){return 11;}
	else if( strcmp(AMIN,"PHE")==0 ){return 11;}
	else if( strcmp(AMIN,"GLN")==0 ){return 9;}
	else if( strcmp(AMIN,"HIS")==0 ){return 10;}
	else if( strcmp(AMIN,"CYS")==0 ){return 6;}
	else {	printf("No such a type: %s\n",AMIN);exit(1);return 0;}


}


double distance(double x1[3], double x2[3])
  {
double dist;


dist = ( (x1[0]-x2[0])*(x1[0]-x2[0])+(x1[1]-x2[1])*(x1[1]-x2[1]) + (x1[2]-x2[2])*(x1[2]-x2[2]) );


return sqrt(dist);
  }

