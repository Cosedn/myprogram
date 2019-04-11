#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include "ccd.h"
#include "mpi.h"

#define USERS 6040 //row
#define ITEMS 3952 //col
#define ROW 12
#define COL 12 //ROW==COL
#define ATTR 5 //attributes

int main(int argc, char *argv[])
{
	FILE *record;
	int i,j,k;
	int t,tg;
	int schedule;
	RatingMatrix *rm_rbuffer, *rm_cbuffer;
	double *u_buffer,*v_buffer,*r_rbuffer,*r_cbuffer;
	double u_temp[USERS]={0},v_temp[ITEMS]={0};
//	int *r_to_c,*c_to_r;
	Block r_rtag[50],r_ctag[50],rl_rtag,rl_ctag;
	int r_list[200];
	int attr=ATTR;
	int rm_size=0,r_size=ROW;
	int u_entries[USERS]={0},v_entries[ITEMS]={0};
	int u_block[ROW+1]={0},v_block[COL+1]={0};
	int u_bsize[ROW]={0},v_bsize[COL]={0};
	int aver;
	
	int loop_continue,loop_counts;
	double sub_rmse,curr_rmse,rmse=100000000.0;
	double rmse_array[ROW]={0};
	int num;
	int begin,end,row0,col0,u0,v0;
	double sum1,sum2;
	struct timeval tv1,tv2;

	double e;
	double lmd1=0.05;
	double lmd2=0.05;
	double temp;
	double t0;
	int r0,c0,r1,c1;
	
	MPI_Status status[400];
	MPI_Status status1[ROW];
	MPI_Request request[400];
	MPI_Request request1[ROW];
	int myid,numprocs;
	
	for(i=0;i<r_size;i++) 
	{
		r_rtag[i].row=i;
		r_rtag[i].col=0;
		r_rtag[i].entries=0;
		r_ctag[i].row=0;
		r_ctag[i].col=i;
		r_ctag[i].entries=0;
	}
	
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	if(myid==0)
	{
		num=Input_File("../../datasets/original/ratings1m.dat", u_entries, v_entries);
		aver=num/ROW;
		Make_Block(u_entries, v_entries, u_block, v_block , USERS, ITEMS, aver);
/*		printf("%d\n",aver);
		for(i=0;i<=ROW;i++) printf("%d ",u_block[i]);
		printf("\n");
		for(i=0;i<=COL;i++) printf("%d ",v_block[i]);
		printf("\n");*/
/*		for(i=0;i<ROW;i++)
		{
			t=0;
			for(j=u_block[i];j<u_block[i+1];j++) t+=u_entries[j];
			printf("%d ",t);
		}
		printf("\n");
		for(i=0;i<COL;i++)
		{
			t=0;
			for(j=v_block[i];j<v_block[i+1];j++) t+=v_entries[j];
			printf("%d ",t);
		}
		printf("\n");*/
		File_Separate("../../datasets/original/ratings1m.dat","../../datasets/transpose/ratings1mc.dat", u_block, v_block, r_rtag, r_ctag); //store r_tag into a file meanwhile
		printf("File Separate Finish!\n");
		
		for(i=0;i<ROW;i++) u_bsize[i]=u_block[i+1]-u_block[i];
		for(i=0;i<COL;i++) v_bsize[i]=v_block[i+1]-v_block[i];
		Get_Rtag("r_rtag.txt",r_rtag);
		Get_Rtag("r_ctag.txt",r_ctag);
/*		for(i=0;i<r_size;i++) printf("%d %d %d\n",r_rtag[i].row,r_rtag[i].col,r_rtag[i].entries);
                for(i=0;i<r_size;i++) printf("%d %d %d\n",r_ctag[i].row,r_ctag[i].col,r_ctag[i].entries);*/
/*		for(i=0;i<ROW;i++) printf("%d ",u_bsize[i]);
		printf("\n");
		for(i=0;i<COL;i++) printf("%d ",v_bsize[i]);
		printf("\n");*/
		
		for(i=0;i<r_size;i++)
		{
			MPI_Send(&r_rtag[i], sizeof(Block), MPI_BYTE, i+1, i+1, MPI_COMM_WORLD);
			MPI_Send(&r_ctag[i], sizeof(Block), MPI_BYTE, i+1, i+1+numprocs, MPI_COMM_WORLD);
		}
		MPI_Bcast(u_bsize, ROW, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(v_bsize, COL, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(u_block, ROW+1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(v_block, COL+1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(u_entries, USERS, MPI_INT, 0, MPI_COMM_WORLD);
                MPI_Bcast(v_entries, ITEMS, MPI_INT, 0, MPI_COMM_WORLD);
		if((record=fopen("record.txt","w"))==NULL)
		{
			printf("cannot write to this file!\n");
			exit(0);
		}
	}
	else
	{
		MPI_Recv(&rl_rtag, sizeof(Block), MPI_BYTE, 0, myid, MPI_COMM_WORLD, &status[myid]);
		MPI_Recv(&rl_ctag, sizeof(Block), MPI_BYTE, 0, myid+numprocs, MPI_COMM_WORLD, &status[2*myid]);
		MPI_Bcast(u_bsize, ROW, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(v_bsize, COL, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(u_block, ROW+1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(v_block, COL+1, MPI_INT, 0, MPI_COMM_WORLD);
                MPI_Bcast(u_entries, USERS, MPI_INT, 0, MPI_COMM_WORLD);
                MPI_Bcast(v_entries, ITEMS, MPI_INT, 0, MPI_COMM_WORLD);
		rm_rbuffer=(RatingMatrix *)calloc(rl_rtag.entries, sizeof(RatingMatrix));
		rm_cbuffer=(RatingMatrix *)calloc(rl_ctag.entries, sizeof(RatingMatrix));
		r_rbuffer=(double *)calloc(rl_rtag.entries, sizeof(double));
		r_cbuffer=(double *)calloc(rl_ctag.entries, sizeof(double));
//		r_to_c=(int *)calloc(rl_rtag.entries, sizeof(int));
//		c_to_r=(int *)calloc(rl_ctag.entries, sizeof(int));

//		printf("%d:%d %d\n",myid,rl_rtag.row,rl_ctag.col);	
		File_Block_Get(rl_rtag, "row", rm_rbuffer);
		File_Block_Get(rl_ctag, "col", rm_cbuffer);

//		Index_Match(u_entries, v_entries, rm_rbuffer, rm_cbuffer, rl_rtag, rl_ctag, u_block, v_block, r_to_c, c_to_r);

		u_buffer=(double *)calloc(USERS*attr,sizeof(double));
		v_buffer=(double *)calloc(ITEMS*attr,sizeof(double));
		
		j=ITEMS*attr;
		for(i=0;i<j;i++) v_buffer[i]=1;
		for(i=0;i<rl_rtag.entries;i++) r_rbuffer[i]=rm_rbuffer[i].rating;
		for(i=0;i<rl_ctag.entries;i++) r_cbuffer[i]=rm_cbuffer[i].rating;
/*		//initialize r_begin_point and r_end_point
		r_begin_point=0;
		for(i=0;i<u_block[myid-1];i++) r_begin_point+=u_entries[i];
		r_end_point=r_begin_point;
		for(i=u_block[myid-1];i<u_block[myid];i++) r_end_point+=u_entries[i];
		//initialize c_begin_point and c_end_point
		c_begin_point=0;
                for(i=0;i<v_block[myid-1];i++) c_begin_point+=v_entries[i];
                c_end_point=c_begin_point;
                for(i=v_block[myid-1];i<v_block[myid];i++) r_end_point+=v_entries[i];*/ //delete
	}

	MPI_Barrier(MPI_COMM_WORLD);
	
	if(myid==0)
	{
		loop_counts=0;
		gettimeofday(&tv1,NULL);
		while(1)
		{
			loop_continue=0;
			sub_rmse=0;
			curr_rmse=0;
			schedule=0;
			MPI_Bcast(&loop_continue, 1, MPI_INT, 0, MPI_COMM_WORLD);
			for(i=0;i<ROW;i++) MPI_Irecv(&rmse_array[i], 1, MPI_DOUBLE, i+1, i+201, MPI_COMM_WORLD, &request1[i]);
			for(i=0;i<ROW;i++) MPI_Wait(&request1[i], &status1[i]);
			for(i=0;i<ROW;i++)
			{
				//printf("%d %.4f\n",i,rmse_array[i]);
				curr_rmse+=rmse_array[i]/num*r_rtag[i].entries;
			}
			curr_rmse=sqrt(curr_rmse);
			fprintf(record,"%.4f\n",curr_rmse);
			if(curr_rmse<rmse && loop_counts<10000 ) rmse=curr_rmse;//
			else
			{
				loop_continue=1;
				MPI_Bcast(&loop_continue, 1, MPI_INT, 0, MPI_COMM_WORLD);
				break;
			}
			loop_counts++;
		}
		gettimeofday(&tv2,NULL);
		fprintf(record, "loop counts:%d\n", loop_counts);
		fprintf(record, "average time per loop:%.4fms\n",((double)(tv2.tv_sec-tv1.tv_sec)*1000+(double)(tv2.tv_usec-tv1.tv_usec)/1000)/loop_counts);
		fclose(record);
	}
	else
	{
		while(1)
		{
			MPI_Bcast(&loop_continue, 1, MPI_INT, 0, MPI_COMM_WORLD);
			if(loop_continue==1) break;
			sub_rmse=0;

			for(i=0;i<attr;i++)
			{
				row0=i*USERS;
				col0=i*ITEMS;

				//printf("%d %d\n",u_block[myid-1],rm_rbuffer[0].user);

				for(j=0;j<USERS;j++) u_temp[j]=u_buffer[row0+j];
				for(j=0;j<ITEMS;j++) v_temp[j]=v_buffer[col0+j];
//				for(j=0;j<USERS;j++) u_temp[j]=0;
//                              for(j=0;j<ITEMS;j++) v_temp[j]=0;

				end=begin=0;
				while(1)
				{
					while(rm_rbuffer[end].user==rm_rbuffer[begin].user) end++;
					u0=row0+rm_rbuffer[begin].user;
					for(j=begin;j<end;j++) 
					{
						v0=col0+rm_rbuffer[j].item;
						r_rbuffer[j]+=u_buffer[u0]*v_buffer[v0];
					}
					begin=end;
					if(end>=rl_rtag.entries) break;
				}
				end=begin=0;
                                while(1)
                                {
                                        while(rm_cbuffer[end].user==rm_cbuffer[begin].user) end++;
                                        v0=col0+rm_cbuffer[begin].user;
                                        for(j=begin;j<end;j++)
                                        {
                                                u0=row0+rm_cbuffer[j].item;
                                                r_cbuffer[j]+=u_buffer[u0]*v_buffer[v0];
                                        }
                                        begin=end;
                                        if(end>=rl_ctag.entries) break;
                                }
				/*for(j=0;j<20;j++) printf("%.4f ",r_rbuffer[j]);
				printf("\n");
				for(j=0;j<20;j++) printf("%.4f ",r_cbuffer[j]);
                                printf("\n");*/

/*				for(j=0;j<rl_rtag.entries;j++) r_rbuffer[j]+=u_buffer[row0+rm_rbuffer[j].user]*v_buffer[col0+rm_rbuffer[j].item];
				for(j=0;j<rl_ctag.entries;j++) r_cbuffer[j]+=u_buffer[row0+rm_cbuffer[j].item]*v_buffer[col0+rm_cbuffer[j].user];*/
				for(k=0;k<1;k++)
				{
					end=begin=0;
					while(1)
					{
						while(rm_rbuffer[end].user==rm_rbuffer[begin].user) end++;
						sum1=sum2=0;
						for(j=begin;j<end;j++)
						{
							sum1+=r_rbuffer[j]*v_temp[rm_rbuffer[j].item];
							sum2+=v_temp[rm_rbuffer[j].item]*v_temp[rm_rbuffer[j].item];
//							sum1+=r_rbuffer[j]*v_buffer[col0+rm_rbuffer[j].item];
//							sum2+=v_buffer[col0+rm_rbuffer[j].item]*v_buffer[col0+rm_rbuffer[j].item];
						}
						u_temp[rm_rbuffer[begin].user]=sum1/(lmd1+sum2);
                                		begin=end;
                                		if(end>=rl_rtag.entries) break;
					}
					//printf("%d %.4f\n",myid,u_temp[u_block[0]]);
					for(j=1;j<numprocs;j++) if(j!=myid) MPI_Irecv(&u_temp[u_block[j-1]], u_bsize[j-1], MPI_DOUBLE, j, myid+j*numprocs, MPI_COMM_WORLD, &request[myid+j*numprocs]);
					for(j=1;j<numprocs;j++) if(j!=myid) MPI_Isend(&u_temp[u_block[myid-1]], u_bsize[myid-1], MPI_DOUBLE, j, j+myid*numprocs, MPI_COMM_WORLD, &request[j+myid*numprocs]);
					for(j=1;j<numprocs;j++) if(j!=myid) MPI_Wait(&request[myid+j*numprocs], &status[myid+j*numprocs]);
					//printf("%d %.4f\n",myid,u_temp[u_block[0]]);

			        	end=begin=0;
                        		while(1)
                        		{
                                		while(rm_cbuffer[end].user==rm_cbuffer[begin].user) end++;
                                        	sum1=sum2=0;
                                        	for(j=begin;j<end;j++)
						{
                                                	sum1+=r_cbuffer[j]*u_temp[rm_cbuffer[j].item];
                                                	sum2+=u_temp[rm_cbuffer[j].item]*u_temp[rm_cbuffer[j].item];
//							sum1+=r_cbuffer[j]*u_buffer[row0+rm_cbuffer[j].item];
//							sum2+=u_buffer[row0+rm_cbuffer[j].item]*u_buffer[row0+rm_cbuffer[j].item];
                                        	}
                                        	v_temp[rm_cbuffer[begin].user]=sum1/(lmd2+sum2);
                                		begin=end;
                                		if(end>=rl_ctag.entries) break;
                                	}
					for(j=1;j<numprocs;j++) if(j!=myid) MPI_Irecv(&v_temp[v_block[j-1]], v_bsize[j-1], MPI_DOUBLE, j, myid+j*numprocs, MPI_COMM_WORLD, &request[myid+j*numprocs]);
					for(j=1;j<numprocs;j++) if(j!=myid) MPI_Isend(&v_temp[v_block[myid-1]], v_bsize[myid-1], MPI_DOUBLE, j, j+myid*numprocs, MPI_COMM_WORLD, &request[j+myid*numprocs]);
                                        for(j=1;j<numprocs;j++) if(j!=myid) MPI_Wait(&request[myid+j*numprocs], &status[myid+j*numprocs]);
				}

/*				for(j=0;j<20;i++) printf("%.4f ",u_temp[j]);
				printf("\n");
				for(j=0;j<20;i++) printf("%.4f ",v_temp[j]);
				printf("\n");*/

				end=begin=0;
                                while(1)
                                {
                                        while(rm_rbuffer[end].user==rm_rbuffer[begin].user) end++;
                                        for(j=begin;j<end;j++)
					{
						r_rbuffer[j]-=u_temp[rm_rbuffer[j].user]*v_temp[rm_rbuffer[j].item];
//						sub_rmse+=r_rbuffer[j]*r_rbuffer[j]/rl_rtag.entries;
					}
                                        begin=end;
                                        if(end>=rl_rtag.entries) break;
                                }
                                end=begin=0;
                                while(1)
                                {
                                        while(rm_cbuffer[end].user==rm_cbuffer[begin].user) end++;
                                        for(j=begin;j<end;j++)
					{
						r_cbuffer[j]-=u_temp[rm_cbuffer[j].item]*v_temp[rm_cbuffer[j].user];
		//				sub_rmse+=r_cbuffer[j]*r_cbuffer[j]/rl_ctag.entries;
					}
                                        begin=end;
                                        if(end>=rl_ctag.entries) break;
                                }
/*				for(j=0;j<rl_rtag.entries;j++)
				{
					r_rbuffer[j]-=u_temp[rm_rbuffer[j].user]*v_temp[rm_rbuffer[j].item];
					sub_rmse+=r_rbuffer[j]*r_rbuffer[j]/rl_rtag.entries;
				}
                                for(j=0;j<rl_ctag.entries;j++) r_cbuffer[j]-=u_temp[rm_cbuffer[j].item]*v_temp[rm_cbuffer[j].user];*/
				for(j=0;j<USERS;j++) u_buffer[row0+j]=u_temp[j];
				for(j=0;j<ITEMS;j++) v_buffer[col0+j]=v_temp[j];
			}

			//printf("%d %.4f\n",myid,sub_rmse);//
			for(i=0;i<rl_rtag.entries;i++) sub_rmse+=r_rbuffer[i]*r_rbuffer[i]/rl_rtag.entries;
			MPI_Isend(&sub_rmse, 1, MPI_DOUBLE, 0, myid+200, MPI_COMM_WORLD, &request1[myid-1]);
		}
	}

	MPI_Finalize();
	
	if(myid!=0)
	{
		free(u_buffer);
		free(v_buffer);
		free(rm_rbuffer);
		free(rm_cbuffer);
		free(r_rbuffer);
		free(r_cbuffer);
	}
	
	return 0;
}
