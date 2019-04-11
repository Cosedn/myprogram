#ifndef CCD_H
#define CCD_H

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

typedef struct
{
	int user;
	int item;
	int rating;
}RatingMatrixTemp;

int Input_File(char *, int *, int *);
void Make_Block(int *, int *, int *, int *, int, int, int);
void File_Block_Get(Block, char *, RatingMatrix *);
void File_Separate(char *, char *, int *, int *, Block *, Block *);
void Store_Rtag(char *, Block *);
void Get_Rtag(char *, Block *);
void File_Transpose(RatingMatrixTemp *, RatingMatrixTemp *, int);
void Index_Match(int *, int *, RatingMatrix *, RatingMatrix *, Block, Block, int *, int *, int *, int *);
void File_Sort(int, Block);

#endif
