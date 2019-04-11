#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "ccd.h"

#define ROW 12
#define COL 12

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
/*	int i;
    int row_number,col_number;
    row_number=user/ROW;
    if(user>row_number*ROW) row_number++;
    col_number=item/COL;
    if(item>col_number*COL) col_number++;
    u_block[0]=0;
    v_block[0]=0;
    u_block[ROW]=user;
    v_block[COL]=item;
    for(i=1;i<ROW;i++) u_block[i]=u_block[i-1]+row_number;
    for(i=1;i<COL;i++) v_block[i]=v_block[i-1]+col_number;*/
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

void File_Block_Get(Block rl_tag, char *type, RatingMatrix *rm_buffer) //Get the corresponding blocks in files for each thread
{
	FILE *fp;
	int sz=0;
	char str[20],seq[10];
	if(strcmp(type,"row")==0)
	{
		strcpy(str,"./data_files/r");
		itoa(rl_tag.row,seq);
		strcat(str,seq);
		if((fp=fopen(str,"r"))==NULL)
                {
                        printf("cannot open this file!\n");
                        exit(0);
                }
                while(fscanf(fp,"%d %d %lf",&rm_buffer[sz].user,&rm_buffer[sz].item,&rm_buffer[sz].rating)!=EOF) sz++;
                fclose(fp);
	}
	else if(strcmp(type,"col")==0)
	{
	        strcpy(str,"./data_files/c");
                itoa(rl_tag.col,seq);
		strcat(str,seq);
                if((fp=fopen(str,"r"))==NULL)
                {
                        printf("cannot open this file!\n");
                        exit(0);
                }
                while(fscanf(fp,"%d %d %lf",&rm_buffer[sz].user,&rm_buffer[sz].item,&rm_buffer[sz].rating)!=EOF) sz++;
                fclose(fp);	
	}
}

void File_Separate(char *str1, char *str2, int *u_block, int *v_block, Block *r_rtag, Block *r_ctag) //Separate the Rating Matrix into blocks and store them in files
{
	FILE *fp1,*fp2,*fq[50],*fr[50];
	int i,j;
	int r,c,p,q;
	int user,item,rating;
	long temp;
	int size;
	char s[50][30];
	char seq[10];
	if((fp1=fopen(str1,"r"))==NULL)
	{
		printf("cannot open this file!\n");
		exit(0);
	}
	if((fp2=fopen(str2,"r"))==NULL)
	{
		printf("cannot open this file!\n");
		exit(0);
	}
	size=ROW;
	for(i=0;i<size;i++)
	{
		strcpy(s[i],"./data_files/r");
		itoa(i,seq);
		strcat(s[i],seq);
		if((fq[i]=fopen(s[i],"w"))==NULL)
		{
			printf("cannot write to this file!\n");
			exit(0);
		}
	}
	for(i=0;i<size;i++)
	{
		strcpy(s[i],"./data_files/c");
		itoa(i,seq);
		strcat(s[i],seq);
		if((fr[i]=fopen(s[i],"w"))==NULL)
		{
			printf("cannot write to this file!\n");
			exit(0);
		}
	}
	while(fscanf(fp1,"%d %d %d %ld",&user,&item,&rating,&temp)!=EOF)
	{
		r=user-1;
		c=item-1;
		for(i=0;i<=ROW;i++) if(r<u_block[i]) break;
		p=i-1;
		fprintf(fq[p],"%d    %d    %d\n",r,c,rating);
		r_rtag[p].entries++;
	}
	while(fscanf(fp2,"%d %d %d",&user,&item,&rating)!=EOF)
	{
		r=user-1;
		c=item-1;
		for(j=0;j<=COL;j++) if(r<v_block[j]) break;
		q=j-1;
		fprintf(fr[q],"%d    %d    %d\n",r,c,rating);
		r_ctag[q].entries++;
	}
	Store_Rtag("r_rtag.txt",r_rtag);
	Store_Rtag("r_ctag.txt",r_ctag);
	fclose(fp1);
	fclose(fp2);
	for(i=0;i<size;i++)
	{
		fclose(fq[i]);
		fclose(fr[i]);
	}
}

void Store_Rtag(char *str, Block *r_tag)
{
	FILE *fp=NULL;
	int size=ROW;
	int i;
	if((fp=fopen(str,"w"))==NULL)
	{
		printf("cannot write to this file!\n");
		exit(0);
	}
	for(i=0;i<size;i++) fprintf(fp,"%d %d %d\n",r_tag[i].row,r_tag[i].col,r_tag[i].entries);
	fclose(fp);
}

void Get_Rtag(char *str, Block *r_tag)
{
	FILE *fp=NULL;
	int size=ROW;
	int i;
	if((fp=fopen(str,"r"))!=NULL)
	{
		for(i=0;i<size;i++) fscanf(fp,"%d %d %d",&r_tag[i].row,&r_tag[i].col,&r_tag[i].entries);
	}
	fclose(fp);
}

void File_Transpose(RatingMatrixTemp *p, RatingMatrixTemp *pt, int size)
{
	int i,j;
	int psize;
	int begin,end;
	RatingMatrixTemp temp;
	for(i=0;i<size;i++)
	{
		pt[i].user=p[i].item;
		pt[i].item=p[i].user;
		pt[i].rating=p[i].rating;
	}
	for(i=0;i<size;i++)
		for(j=size-1;j>i;j--)
		{
			if(pt[j].user < pt[j-1].user)
			{
				temp=pt[j];
				pt[j]=pt[j-1];
				pt[j-1]=temp;
			}
		}
	end=begin=0;
	while(1)
	{
		while(pt[end].item==pt[begin].item) end++;
		for(i=begin;i<end;i++)
			for(j=end-1;j>i;j--)
			{
				if(pt[j].item < pt[j-1].item)
				{
					temp=pt[j];
					pt[j]=pt[j-1];
					pt[j-1]=temp;
				}
			}
		begin=end;
		if(end>=size) break;
	}
}

void Index_Match(int *u_entries, int *v_entries, RatingMatrix *rm_rbuffer, RatingMatrix *rm_cbuffer, Block rl_rtag, Block rl_ctag, int *u_block, int *v_block, int *r_to_c, int *c_to_r)
{
	int i,j;
	int c0,c1,c2;
	c0=rl_ctag.col;
	c1=c0+1;
	for(i=0;i<rl_rtag.entries;i++)
	{
		if(rm_rbuffer[i].item>=v_block[c0] && rm_rbuffer[i].item<v_block[c1])
		{
			c2=rm_rbuffer[i].item;
			j=0;
			while(rm_cbuffer[j].user!=c2) j+=v_entries[rm_cbuffer[j].user];
			while(1)
			{
				if(rm_cbuffer[j].item==rm_rbuffer[i].user) break;
				else j++;
			}
//			if(j==rl_ctag.entries) printf("fail!\n");
			r_to_c[i]=j;
			c_to_r[j]=i;
/*			printf("%d %d\n",rm_rbuffer[i].user,rm_rbuffer[i].item);
			printf("%d %d\n",rm_cbuffer[j].user,rm_cbuffer[j].item);
			printf("\n");*/
		}
	}
}

void File_Sort(int id, Block rl_ctag) //execute in several processes
{
	FILE *fq,*fr;
	int j;
	char s[30];
	char seq[10];
	RatingMatrixTemp *rmt1,*rmt2;
        strcpy(s,"./data_files/cc");
        itoa(id,seq);
        strcat(s,seq);
        if((fq=fopen(s,"r"))==NULL)
        {
                printf("cannot write to this file!\n");
                exit(0);
        }

        strcpy(s,"./data_files/c");
//	itoa(id,seq);
        strcat(s,seq);
        if((fr=fopen(s,"w"))==NULL)
        {
                printf("cannot write to this file!\n");
                exit(0);
        }

        rmt1=(RatingMatrixTemp *)calloc(rl_ctag.entries, sizeof(RatingMatrixTemp));
        rmt2=(RatingMatrixTemp *)calloc(rl_ctag.entries, sizeof(RatingMatrixTemp));
        for(j=0;j<rl_ctag.entries;j++) fscanf(fq, "%d %d %d", &rmt1[j].user, &rmt1[j].item, &rmt1[j].rating);
        File_Transpose(rmt1, rmt2, rl_ctag.entries);
        for(j=0;j<rl_ctag.entries;j++) fprintf(fr, "%d             %d              %d\n", rmt2[j].user, rmt2[j].item, rmt2[j].rating);
        free(rmt1);
        free(rmt2);

        fclose(fq);
        fclose(fr);
        strcpy(s,"./data_files/cc");
//	itoa(id,seq);
        strcat(s,seq);
        remove(s);
}
