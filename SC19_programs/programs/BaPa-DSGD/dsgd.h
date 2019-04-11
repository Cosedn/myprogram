#ifndef DSGD_H
#define DSGD_H

typedef struct
{
	int user;
	int item;
	double rating;
}RatingMatrix;

typedef struct
{
	int row;
	int col;
	int entries;
}Block;

int Input_File(char *, int *, int *);
void Make_Block(int *, int *, int *, int *, int, int, int);
void File_Block_Get(Block *, int, RatingMatrix *);
void File_Separate(char *, int *, int *, Block *);
void Store_Rtag(Block *);
void Get_Rtag(Block *);

#endif
