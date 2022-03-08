#include<stdio.h>
#include<string.h>
#include<dirent.h>
int readstring(char *);
main (int argc, char **argv){
	int i,sec_arg,k1,k2,nalig,kpr,kk,nn,ys,xc,yc,ycn;
	DIR *d;
	FILE *f,*f1,*f2;
	char fname[50],fnpdb[80],fnpdb2[80],fnpdb1[15],com[512],nlig[4],chnm[2],
	fnpr[80],fhtml[80],fhtml1[16],ftxt[80],ftxt1[15],flgfpdb[11],
	fnn[256],ss[256],ftt[5],namedir[3],str[8],abc[80],sst[35],
	s19[20],*progpath,*libpath,*pwd,fnpdb3[15],fnpdb4[15],*http,
	tz[150];


/*	if(argc!=4){
		printf("The program requires three arguments!!!\n");
		exit(1);
	}

	if((libpath=(char *)getenv("COFLIBPATH"))==NULL)
        {
		printf("You have to define COFLIBPATH.\n");
		exit(1);
        }
*/
	sec_arg=readstring(argv[2]);
	chnm[0]=' ';
	chnm[1]='\0';
	if (argc == 4) {
	sprintf(chnm,"%s",argv[3]);}

	rescprep_(argv[1],&sec_arg,&chnm,&k1,&k2);
	resccal_(&sec_arg,&chnm,&k1,&k2);
}int readstring(char *S){
	int i,l=strlen(S);
	int sum=0;
	for(i=0;i<l;i++)
		if(S[i]>='0'&&S[i]<='9')

			sum=sum*10+S[i]-'0';
		else{
			printf("The first argument should be numerical\n");
			exit(1);
		}	

	return(sum);
}
