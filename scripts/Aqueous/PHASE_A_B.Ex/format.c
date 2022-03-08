/* format.c
	
	This program converts the xyz file to amber inputcrd format.

	a.out input.crd
*/
	
#include <stdio.h>


char *inputfile;

main(int argc, char *argv[])
{

FILE *fp1, *fp2;
float d;
int i,j,k;
int Ntot;
int count;
char word[1];

 if(argc<2){printf("Usage: a.out input.crd \n"); exit(1);}

	inputfile=argv[1];

fp1= fopen(inputfile,"r");

	fscanf(fp1,"%s",word);
	printf("%s\n",word);
	fscanf(fp1,"%d",&Ntot);
	printf("   %d\n",Ntot);

for(i=1;i<(Ntot*3+1);i++)
{
                fscanf(fp1,"%f", &d);
                printf("%12.7f",d);
	if(i>0 && i%6==0){ printf("\n");}
}

fclose(fp1);

}

