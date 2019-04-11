#include <stdio.h>
#include "dsgd.h"

#define ROW 13
#define COL 13

void itoa(int num, char* str)
{
        int i=0,j;
        int t;
        char str1[10];
        while(1)
        {
                if(num==0) break;
                t=num-num/10*10;
                num=(num-t)/10;
                str1[i]=(char)(t+48);
                i++;
        }
        for(j=0;j<i;j++) str[j]=str1[i-j-1];
        str[j]='\0';
}

int Input_File(char* str, int *u_entries, int *v_entries)
{
	FILE *fp;
	int user,item;
	int rating;
	long temp;
	int num=0;
	if((fp=fopen(str,"r"))==NULL)
	{
		printf("cannot open this file!\n");
		exit(0);
	}
    while(fscanf(fp,"%d %d %d %ld",&user,&item,&rating,&temp)!=EOF)
	{
		u_entries[user-1]++;
		v_entries[item-1]++;
		num++;
	}
	fclose(fp);
	return num;
}

void Make_Block(int *u_entries, int *v_entries, int *u_block, int *v_block, int user, int item, int aver)
{
	int i,j;
	int t0,t1,temp;
	u_block[0]=0;
	i=1;
	j=0;
	temp=0;
	t0=aver;
	while(j<user)
	{
		temp+=u_entries[j];
		t1=(temp>aver)?(temp-aver):(aver-temp);
		if(t1<=t0)
		{
			t0=t1;
			j++;
		}
		else
		{
			u_block[i]=j;
			i++;
			temp=0;
			t0=aver;
		}
		if(i==ROW)
		{
			u_block[i]=user;
			break;
		}
	}
	v_block[0]=0;
	i=1;
	j=0;
	temp=0;
	t0=aver;
	while(j<item)
	{
		temp+=v_entries[j];
		t1=(temp>aver)?(temp-aver):(aver-temp);
		if(t1<=t0)
		{
			t0=t1;
			j++;
		}
		else
		{
			v_block[i]=j;
			i++;
			temp=0;
			t0=aver;
		}
		if(i==COL)
		{
			v_block[i]=item;
			break;
		}
	}
}

void File_Block_Get(Block *rl_tag, int rl_size, RatingMatrix *rm_buffer) //Get the corresponding blocks in files for each thread
{
	FILE *fp;
	int i;
	int sz=0;
	char str[20],seq[3];
	for(i=0;i<rl_size;i++)
	{
//			int d=0;//
        strcpy(str,"./data_files/r");
/*		seq[0]=(char)(rl_tag[i].row+48);
		seq[1]=(char)(rl_tag[i].col+48);
		seq[2]='\0';*/
		itoa(rl_tag[i].row*COL+rl_tag[i].col,seq);
		strcat(str,seq);
		if((fp=fopen(str,"r"))==NULL)
		{
			printf("cannot open this file!\n");
			exit(0);
		}
		while(fscanf(fp,"%d %d %lf",&rm_buffer[sz].user,&rm_buffer[sz].item,&rm_buffer[sz].rating)!=EOF) sz++;//,d++;//
//		printf("%s:%d ",seq,d);//
		fclose(fp);
	}
}

void File_Separate(char *str, int *u_block, int *v_block, Block *r_tag) //Separate the Rating Matrix into blocks and store them in files
{
	FILE *fp,*fq[400];
	int i,j;
	int r,c,p;
	int user,item,rating;
	long temp;
	int size;
	char s[400][20];
	char seq[10];
	if((fp=fopen(str,"r"))==NULL)
	{
		printf("cannot open this file!\n");
		exit(0);
	}
	size=ROW*COL;
	for(i=0;i<size;i++)
	{
		strcpy(s[i],"./data_files/r");
/*		seq[0]=(char)(i/COL+48);
		seq[1]=(char)(i-(i/COL)*COL+48);
		seq[2]='\0';*/
		itoa(i,seq);
		strcat(s[i],seq);
		if((fq[i]=fopen(s[i],"w"))==NULL)
		{
			printf("cannot write to this file!\n");
			exit(0);
		}
	}
    while(fscanf(fp,"%d %d %d %ld",&user,&item,&rating,&temp)!=EOF)
	{
		r=user-1;
		c=item-1;
		for(i=0;i<=ROW;i++) if(r<u_block[i]) break;
		for(j=0;j<=COL;j++) if(c<v_block[j]) break;
		p=(i-1)*COL+j-1;
		fprintf(fq[p],"%d    %d    %d\n",r,c,rating);
		r_tag[p].entries++;
	}
	Store_Rtag(r_tag);
	fclose(fp);
	for(i=0;i<size;i++) fclose(fq[i]);
}

void Store_Rtag(Block *r_tag)
{
	FILE *fp=NULL;
	int size=ROW*COL;
	int i;
	if((fp=fopen("rtag.txt","w"))==NULL)
	{
		printf("cannot write to this file!\n");
		exit(0);
	}
	for(i=0;i<size;i++) fprintf(fp,"%d %d %d\n",r_tag[i].row,r_tag[i].col,r_tag[i].entries);
	fclose(fp);
}

void Get_Rtag(Block *r_tag)
{
	FILE *fp=NULL;
	int size=ROW*COL;
	int i;
	if((fp=fopen("rtag.txt","r"))!=NULL)
	{
		for(i=0;i<size;i++) fscanf(fp,"%d %d %d",&r_tag[i].row,&r_tag[i].col,&r_tag[i].entries);
	}
	fclose(fp);
}
