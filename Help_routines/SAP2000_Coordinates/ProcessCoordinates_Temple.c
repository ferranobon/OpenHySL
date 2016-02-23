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

int find_value( const float a, const float b, const float c, Substructure_t *const Sub, const int length, const int startIndex );
void Read_TXE_File (const char* FileName, const int NumEntries, JointInfo_t *const Info);
void Generate_SubEntries ( const int NumEntries, const JointInfo_t const* Info, const int NumJointSub, const int NumDrums, const Substructure_t *const Sub, SubEntry_t *const SubEntries );
void Reorder_SubEntries ( const int NumJointSub, SubEntry_t *const SubEntries );
void Print_CoupleNodes ( const int NumJointSub, const int NumDrums, const SubEntry_t *const SubEntries );
void Reorder_JointEntries ( const int NumEntries, JointInfo_t *const JointEntries );

int main (int argc, char **argv)
{
     FILE *TheFile;
     char *FileName;
     
     char *line = NULL;
     char dummy[40];
     size_t len = 0;
     ssize_t read;

     int JointID, i, j, k, NumJointSub, pos, NumEntries = 396;
     float XCoord, YCoord, ZCoord;
     float XValues[6] = {0.0, 4.3, 8.76, 13.22, 17.68, 21.98};
     float YValues[14] = {0.0, 4.3, 8.76, 13.22, 17.68, 22.14, 26.6, 31.06, 35.52, 39.98, 44.44, 48.9, 53.36, 57.66};
     float ZValues[9] = {4.59, 5.91, 7.23, 8.55, 9.87, 11.21, 12.12, 12.56, 12.785};

     Substructure_t *Sub, temp;
     JointInfo_t *JointEq;
     SubEntry_t *SubEntries;
     

     TheFile = fopen("Coordinates.txt", "r");
     if (TheFile == NULL){
	  exit(EXIT_FAILURE);
     }

     /* Get the number of substructures */
     NumJointSub = 0;
     while ((read = getline(&line, &len, TheFile)) != -1){
	  sscanf( line, "%d %s %s %s %s %s %s %f %f %f", &JointID, dummy, dummy, dummy, dummy, dummy, dummy, &XCoord, &YCoord, &ZCoord );
	  if( (ZCoord > 3.931) && (ZCoord < 14.28) ){
	       /* This is a substructure */
	       NumJointSub = NumJointSub + 1;
	       printf( "Joint %d is in Coord (%f,%f,%f)\n", JointID, XCoord, YCoord, ZCoord );
	  }
     }

     fclose( TheFile );

     if (NumJointSub != (2*14*9 + 2*4*9)){
	  printf("Error. Number of substructures is %d but it should be %d.\n", NumJointSub, 2*14*9 + 2*4*9);
	  exit( EXIT_FAILURE);
     };

     Sub = (Substructure_t *) calloc( (size_t) NumJointSub, sizeof(Substructure_t) );
     i = 0;
     /* Save all substructures of the file */

     TheFile = fopen("Coordinates.txt", "r");
     if (TheFile == NULL){
	  exit(EXIT_FAILURE);
     }
     while ((read = getline(&line, &len, TheFile)) != -1){
	  sscanf( line, "%d %s %s %s %s %s %s %f %f %f", &JointID, dummy, dummy, dummy, dummy, dummy, dummy, &XCoord, &YCoord, &ZCoord );
	  if( (ZCoord > 3.931) && (ZCoord < 14.28) ){
	       /* This is a substructure */
	       Sub[i].JointID = JointID;
	       Sub[i].X = XCoord;
	       Sub[i].Y = YCoord;
	       Sub[i].Z = ZCoord;
	       i = i + 1;
	  }
	      
     }

     /* Reorder the substructures */
     int Column = 0;
     for( i = 0; i < 6; i++ ){
	  for(j = 0; j < 14; j++){
	       if( (i != 0 || i != 5) && (j != 0 || j != 13)){
		    
	       } else {
		    for (k = 0; k < 9; k++){
			 pos = find_value(XValues[i], YValues[j], ZValues[k], Sub, NumJointSub, Column*9 + k);
			 temp = Sub[Column*9 + k];
			 Sub[Column*9 + j] = Sub[pos];
			 Sub[pos] = temp;
		    }
	       }
	       Column = Column + 1;
	  }
     }

     for ( i = 0; i < NumJointSub; i++){
	  printf("Joint number: %d\n", Sub[i].JointID);
     }
    
     /* Read the TXE file with Joint/equation information */
     JointEq = (JointInfo_t *) calloc( (size_t) NumEntries, sizeof(JointInfo_t) );
     Read_TXE_File ("Joints info.TXE", NumEntries, JointEq);
     Reorder_JointEntries( NumEntries, JointEq);

     /* Generate substructure entries */
     SubEntries = (SubEntry_t *) calloc( (size_t) NumJointSub, sizeof(SubEntry_t) );
     Generate_SubEntries ( NumEntries, JointEq, NumJointSub, 9, Sub, SubEntries );

     /* Reorder entries */
     Reorder_SubEntries ( NumJointSub, SubEntries );

     /* Print Couple_Nodes.txt */
     Print_CoupleNodes ( NumJointSub, 9, SubEntries );

     free( JointEq );
     free( Sub );
     free( SubEntries );
     free( line );

     return 0;
}

int find_value( const float a, const float b, const float c, Substructure_t *const Sub, const int length, const int startIndex )
{
     int pos;
     bool found = false;

     pos = startIndex;
     while( pos < length && !found ){
	  if( (Sub[pos].X == a) && (Sub[pos].Y == b) && (Sub[pos].Z == c) ){
	       found = true;
	  } else {
	       pos = pos + 1;
	  }
     };

     if(!found){
	  printf("Value %f %f %f not found.\n", a, b, c);
	  exit(EXIT_FAILURE);
     }
     
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
	  //printf("Read Joint %i %i %i %i %i %i %i\n", Info[i].JointID, Info[i].x, Info[i].y, Info[i].z, Info[i].rotx, Info[i].roty, Info[i].rotz );
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
	  if (SubEntries[i].EqNum == 0 || SubEntries[i].EqNum == -1){
	       printf("Error assigning Eq. num to joint entry %d since %d\n", Sub[i].JointID, Info[Sub[i].JointID - 1].JointID);
	       exit(EXIT_FAILURE);
	  }
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

void Reorder_JointEntries ( const int NumEntries, JointInfo_t *const JointEntries )
{
     int j,bound = NumEntries - 1;
     JointInfo_t temp;
     int i;

     for( i = 0 ; i < NumEntries - 2 ; i++ ){
	 for ( j = 0 ; j < bound ; j++ ) {
	      if ( JointEntries[j].JointID > JointEntries[j + 1].JointID ) {
		   temp = JointEntries[j];
		   JointEntries[j] = JointEntries[j + 1];
		   JointEntries[j + 1] = temp;
	      }            
	 }
	 bound = j-1;
     }

     for (i = 0; i < NumEntries; i++){
	  printf("Joint label %d in position %d\n", JointEntries[i].JointID, i);
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
