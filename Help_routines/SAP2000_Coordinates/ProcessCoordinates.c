#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

typedef struct Substructure {
     int JointID;
     float X, Y, Z;
} Substructure_t;

typedef struct JointInfo {
     int JointID;
     int x, y, z, rotx, roty, rotz;
} JointInfo_t;

typedef struct SubEntry {
     int EqNum, PrevEqNum, Column, Drum;
} SubEntry_t;

int find_value( const float a, const float b, Substructure_t *const Sub, const int length, const int startIndex );
void Read_TXE_File (const char* FileName, const int NumEntries, JointInfo_t *const Info);
void Generate_SubEntries ( const int NumEntries, const JointInfo_t const* Info, const int NumJointSub, const int NumDrums, const Substructure_t *const Sub, SubEntry_t *const SubEntries );
void Reorder_SubEntries ( const int NumJointSub, SubEntry_t *const SubEntries );
void Print_CoupleNodes ( const int NumJointSub, const int NumDrums, const SubEntry_t *const SubEntries );

int main (int argc, char **argv)
{
     FILE *TheFile;
     char *FileName;
     
     char *line = NULL;
     char dummy[40];
     size_t len = 0;
     ssize_t read;

     int JointID, i, j, NumJointSub, pos, NumEntries = 10;
     float XCoord, YCoord, ZCoord;
     float XValues[2] = {0.0, 446.0}, ZValues[4] = {459.0, 591.0, 723.0, 855.0};
     Substructure_t *Sub, temp;
     JointInfo_t *JointEq;
     SubEntry_t *SubEntries;
     

     TheFile = fopen("Test.s2k", "r");
     if (TheFile == NULL){
	  exit(EXIT_FAILURE);
     }

     /* Get the number of substructures */
     NumJointSub = 0;
     while ((read = getline(&line, &len, TheFile)) != -1){
	  sscanf( line, "%d %s %s %s %s %s %s %f %f %f", &JointID, dummy, dummy, dummy, dummy, dummy, dummy, &XCoord, &YCoord, &ZCoord );
	  if( ZCoord > 327.0 ){
	       /* This is a substructure */
	       NumJointSub = NumJointSub + 1;
	       printf( "Joint %d is in Coord (%f,%f,%f)\n", JointID, XCoord, YCoord, ZCoord );
	  }
     }

     fclose( TheFile );

     if (NumJointSub != 2*4){
	  printf("Error. Number of substructures is %d but it should be %d.\n", NumJointSub, 4*2);
	  exit( EXIT_FAILURE);
     };

     Sub = (Substructure_t *) calloc( (size_t) NumJointSub, sizeof(Substructure_t) );
     i = 0;
     /* Save all substructures of the file */

     TheFile = fopen("Test.s2k", "r");
     if (TheFile == NULL){
	  exit(EXIT_FAILURE);
     }
     while ((read = getline(&line, &len, TheFile)) != -1){
	  sscanf( line, "%d %s %s %s %s %s %s %f %f %f", &JointID, dummy, dummy, dummy, dummy, dummy, dummy, &XCoord, &YCoord, &ZCoord );
	  if( ZCoord > 327.0 ){
	       /* This is a substructure */
	       Sub[i].JointID = JointID;
	       Sub[i].X = XCoord;
	       Sub[i].Y = YCoord;
	       Sub[i].Z = ZCoord;
	       i = i + 1;
	  }
	      
     }

     /* Reorder the substructures */
     for( i = 0; i < 2; i++ ){
	  for(j = 0; j < 4; j++){
	       pos = find_value(XValues[i], ZValues[j], Sub, NumJointSub, i*4 + j);
	       temp = Sub[i*4 + j];
	       Sub[i*4 + j] = Sub[pos];
	       Sub[pos] = temp;
	  }
     }

     for ( i = 0; i < NumJointSub; i++){
	  printf("Joint number: %d\n", Sub[i].JointID);
     }
     
     /* Read the TXE file with Joint/equation information */
     JointEq = (JointInfo_t *) calloc( (size_t) NumEntries, sizeof(JointInfo_t) );
     Read_TXE_File ("Joints info.TXE", NumEntries, JointEq);

     /* Generate substructure entries */
     SubEntries = (SubEntry_t *) calloc( (size_t) NumJointSub, sizeof(SubEntry_t) );
     Generate_SubEntries ( NumEntries, JointEq, NumJointSub, 4, Sub, SubEntries );

     /* Reorder entries */
     Reorder_SubEntries ( NumJointSub, SubEntries );

     /* Print Couple_Nodes.txt */
     Print_CoupleNodes ( NumJointSub, 4, SubEntries );

     free( JointEq );
     free( Sub );
     free( SubEntries );
     free( line );

     return 0;
}

int find_value( const float a, const float b, Substructure_t *const Sub, const int length, const int startIndex )
{
     int pos;
     bool found = false;

     pos = startIndex;
     while( pos < length && !found ){
	  if( Sub[pos].X == a && Sub[pos].Z == b ){
	       found = true;
	  } else {
	       pos = pos + 1;
	  }
     };
     
     return pos;
}

void Read_TXE_File (const char* FileName, const int NumEntries, JointInfo_t *const Info)
{
     FILE *InFile;
     char Header[100];
     int i;
     
     InFile = fopen( FileName, "r" );

     if (InFile == NULL){
	  printf( "Could not open %s\n", FileName );
	  exit(EXIT_FAILURE);
     }
     
     fgets( Header, (size_t) 100, InFile );
     while( fscanf( InFile, "%i %i %i %i %i %i %i", &Info[i].JointID, &Info[i].x, &Info[i].y, &Info[i].z, &Info[i].rotx, &Info[i].roty, &Info[i].rotz ) != EOF ){
	  printf("Read Joint %i %i %i %i %i %i %i\n", Info[i].JointID, Info[i].x, Info[i].y, Info[i].z, Info[i].rotx, Info[i].roty, Info[i].rotz );
	  i = i + 1;
     }

     fclose( InFile );
}

void Generate_SubEntries ( const int NumEntries, const JointInfo_t const* Info, const int NumJointSub, const int NumDrums, const Substructure_t *const Sub, SubEntry_t *const SubEntries )
{
     int i, Column, Drum;
     Column = 0; Drum = 0;
     
     for (i = 0; i < NumJointSub; i++){
	  SubEntries[i].EqNum = Info[Sub[i].JointID - 1].x;
	  if ( i % NumDrums == 0 ){
	       SubEntries[i].PrevEqNum = -1;
	       Column = Column + 1;
	       Drum = 1;
	  } else {
	       SubEntries[i].PrevEqNum = Info[Sub[i-1].JointID - 1].x;
	       Drum = Drum + 1;
	  }
	  SubEntries[i].Column = Column;
	  SubEntries[i].Drum = Drum;
	  printf("Eq number %i has as previous %i.\n", SubEntries[i].EqNum, SubEntries[i].PrevEqNum );
     }	  
}

void Reorder_SubEntries ( const int NumJointSub, SubEntry_t *const SubEntries )
{
     int j,bound = NumJointSub - 1;
     SubEntry_t temp;
     int i;

     for( i = 0 ; i < NumJointSub - 2 ; i++ ){
	 for ( j = 0 ; j < bound ; j++ ) {
	      if ( SubEntries[j].EqNum > SubEntries[j + 1].EqNum ) {
		   temp = SubEntries[j];
		   SubEntries[j] = SubEntries[j + 1];
		   SubEntries[j + 1] = temp;
	      }            
	 }
	 bound = j-1;
     }	  
}

void Print_CoupleNodes ( const int NumJointSub, const int NumDrums, const SubEntry_t *const SubEntries )
{
     FILE *TheFile;
     int i;

     TheFile = fopen("Couple Nodes.txt", "w");

     fprintf( TheFile, "%d\n", NumJointSub );
     
     for (i = 0; i < NumJointSub; i++ ){	  
	  fprintf(TheFile, "Sim_StoneDrums, 1 %d, 1 %d, Column: %d, Drum: %d;\n", SubEntries[i].EqNum, SubEntries[i].PrevEqNum, SubEntries[i].Column, SubEntries[i].Drum);
     }

     fclose( TheFile );
};