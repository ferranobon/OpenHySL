#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include <getopt.h>


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
void Read_TXE_File ( FILE *TheFile, const int NumEntries, JointInfo_t *const Info);
void Generate_SubEntries ( const int NumEntries, const JointInfo_t const* Info, const int NumJointSub, const int NumDrums, const Substructure_t *const Sub, SubEntry_t *const SubEntries );
void Reorder_SubEntries ( const int NumJointSub, SubEntry_t *const SubEntries );
void Print_CoupleNodes ( const int NumJointSub, const int NumDrums, const SubEntry_t *const SubEntries, FILE *TheFile );
void Reorder_JointEntries ( const int NumEntries, JointInfo_t *const JointEntries );
int Find_EqPos ( const int JointID, const int NumEntries, const JointInfo_t const* Info );
void print_help( char **argv, const struct option *long_opts );

int main (int argc, char **argv)
{
     FILE *InFile, *TXEFile;
     FILE *OutFile;
     char *InFileName;
     
     char *line = NULL;
     char dummy[40];
     size_t len = 0;
     ssize_t read;

     int JointID, i, j, k, NumJointSub, pos, numEntries;// = 396;
     float XCoord, YCoord, ZCoord;
     float XValues[6] = {0.0, 4.3, 8.76, 13.22, 17.68, 21.98};
     float YValues[14] = {0.0, 4.3, 8.76, 13.22, 17.68, 22.14, 26.6, 31.06, 35.52, 39.98, 44.44, 48.9, 53.36, 57.66};
     float ZValues[9] = {4.59, 5.91, 7.23, 8.55, 9.87, 11.21, 12.12, 12.56, 12.785};

     Substructure_t *Sub, temp;
     JointInfo_t *JointEq;
     SubEntry_t *SubEntries;

     int rc, idx;

     struct option long_options[] = 
     {
	  {"input-file", required_argument, 0, 'i'},
	  {"output-file", required_argument, 0, 'o'},
	  {"num-entries", required_argument, 0, 'n'},
	  {"txe-file", required_argument, 0, 't'},
	  {"help", no_argument, 0, 'h'},
	  {0, 0, 0, 0}
     };

     if( argc <= 1 && argc != 4 ){
	  fprintf( stderr, "Invalid number of arguments.\n" );
	  print_help( argv, long_options );
	  return( EXIT_FAILURE );
     }	  

     while( (rc = getopt_long( argc, argv, "i:o:n:t:h", long_options, &idx )) != -1 ){

	  switch( rc ){
	  case 'i':
	       InFile = fopen( optarg, "r" );
	       InFileName = optarg;
	       if( InFile == NULL ){
		    fprintf( stderr, "Could not open file %s.\n", optarg );
		    fprintf( stderr, "Exiting.\n" );
		    exit( EXIT_FAILURE );
	       }
	       break;
	  case 'o':
	       OutFile = fopen( optarg, "w" );
	       if( OutFile == NULL ){
		    fprintf( stderr, "Could not open file %s.\n", optarg );
		    fprintf( stderr, "Exiting.\n" );
		    exit( EXIT_FAILURE );
	       }
	       break;
	  case 'n':
	       numEntries = atoi( optarg );
	       if( numEntries <= 0 ){
		    fprintf( stderr, "The number of entries is equal or less than 0" );
		    fprintf( stderr, "Exiting.\n" );
		    exit( EXIT_FAILURE );
	       }
	       break;
	  case 't':
	       TXEFile = fopen( optarg, "r" );
	       if( TXEFile == NULL ){
		    fprintf( stderr, "Could not open file %s.\n", optarg );
		    fprintf( stderr, "Exiting.\n" );
		    exit( EXIT_FAILURE );
	       }
	       break;
	  case 'h':
	       print_help( argv, long_options );
	       exit( EXIT_FAILURE );
	       break;
	  case '?':
	       /* Long options already prints an error message telling that there is an unrecognised option */
	       break;
	  case ':':
	       /* Long options already prints an error message telling that the option requieres an argument */
	       break;
	  }
     }

     /* Get the number of joints with substructures */
     NumJointSub = 0;
     while ((read = getline(&line, &len, InFile)) != -1){
	  sscanf( line, "%d %s %s %s %s %s %s %f %f %f", &JointID, dummy, dummy, dummy, dummy, dummy, dummy, &XCoord, &YCoord, &ZCoord );
	  if( (ZCoord > 3.931) && (ZCoord < 14.28) ){
	       /* This is a substructure */
	       NumJointSub = NumJointSub + 1;
	       printf( "Joint %d is in Coord (%f,%f,%f)\n", JointID, XCoord, YCoord, ZCoord );
	  }
     }

     fclose( InFile );

     /* Sanity check */
     if (NumJointSub != (2*14*9 + 2*4*9)){
	  printf("Error. Number of substructures is %d but it should be %d.\n", NumJointSub, 2*14*9 + 2*4*9);
	  exit( EXIT_FAILURE);
     };

     Sub = (Substructure_t *) calloc( (size_t) NumJointSub, sizeof(Substructure_t) );
     i = 0;
     /* Save all substructures of the file */

     InFile = fopen( InFileName, "r");
     if( InFile == NULL ){
	  fprintf( stderr, "Could not open file %s.\n", InFileName );
	  fprintf( stderr, "Exiting.\n" );
	  exit( EXIT_FAILURE );
     }
     while ((read = getline(&line, &len, InFile)) != -1){
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

     fclose(InFile);
     
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
    
     /* Read the TXE file with Joint/equation information */
     JointEq = (JointInfo_t *) calloc( (size_t) numEntries, sizeof(JointInfo_t) );   /* Memory allocation */
     Read_TXE_File ( TXEFile, numEntries, JointEq);                                  /* Read entries */
     fclose ( TXEFile );                                                             /* Close the TXE file */
     Reorder_JointEntries( numEntries, JointEq);                                     /* Reorder entries */
     
     
     
     /* Generate substructure entries */
     SubEntries = (SubEntry_t *) calloc( (size_t) NumJointSub*3, sizeof(SubEntry_t) );  /* Memory allocation */
     Generate_SubEntries ( numEntries, JointEq, NumJointSub*3, 9, Sub, SubEntries );    /* Substructure entries generation */
     Reorder_SubEntries ( NumJointSub*3, SubEntries );                                  /* Reorder entries */

     /* Sanity check. No double entries */
     for (i = 0; i < NumJointSub*3 -1; i++ ){                    
	  if( SubEntries[i].EqNum == SubEntries[i+1].EqNum ){
	       printf("Critical Error: Entry %d has the same equation assigned (Eq %d) as entry %d\n", i, SubEntries[i].EqNum, i+1 );
	       exit( EXIT_FAILURE );
	  }
     }

     /* Print Couple_Nodes.txt */
     Print_CoupleNodes ( NumJointSub*3, 9, SubEntries, OutFile );    
     fclose( OutFile );

     /* Memory deallocation */
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

void Read_TXE_File (FILE *TheFile, const int numEntries, JointInfo_t *const Info)
{
     char Header[100];
     int i;
     
     fgets( Header, (size_t) 100, TheFile );
     while( fscanf( TheFile, "%i %i %i %i %i %i %i", &Info[i].JointID, &Info[i].x, &Info[i].y, &Info[i].z, &Info[i].rotx, &Info[i].roty, &Info[i].rotz ) != EOF ){
	  i = i + 1;
     }
     
}

void Generate_SubEntries ( const int numEntries, const JointInfo_t const* Info, const int NumJointSub, const int NumDrums, const Substructure_t *const Sub, SubEntry_t *const SubEntries )
{
     int i, pos, posprev, Column, Drum;
     Column = 0; Drum = 0;
     
     for (i = 0; i < NumJointSub/3; i++){
	  pos = Find_EqPos( Sub[i].JointID, numEntries, Info);
	  SubEntries[i].EqNum = Info[pos].x;
	  SubEntries[NumJointSub/3 + i].EqNum = Info[pos].y;
	  SubEntries[NumJointSub*2/3 + i].EqNum = Info[pos].z;
	  if (SubEntries[i].EqNum == 0 || SubEntries[i].EqNum == -1){
	       printf("Error assigning Eq. num to joint entry %d since %d\n", Sub[i].JointID, Info[pos].JointID);
	       exit(EXIT_FAILURE);
	  }
	  if (SubEntries[NumJointSub/3 + i].EqNum == 0 || SubEntries[NumJointSub/3 + i].EqNum == -1){
	       printf("Error assigning Eq. num to joint entry %d since %d\n", Sub[i].JointID, Info[pos].JointID);
	       exit(EXIT_FAILURE);
	  }
	  if (SubEntries[NumJointSub*2/3 + i].EqNum == 0 || SubEntries[NumJointSub*2/3 + i].EqNum == -1){
	       printf("Error assigning Eq. num to joint entry %d since %d\n", Sub[i].JointID, Info[pos].JointID);
	       exit(EXIT_FAILURE);
	  }	    
	  if ( i % NumDrums == 0 ){
	       SubEntries[i].PrevEqNum = -1;
	       SubEntries[NumJointSub/3 + i].PrevEqNum = -1;
	       SubEntries[NumJointSub*2/3 + i].PrevEqNum = -1;
	       Column = Column + 1;
	       Drum = 1;
	  } else {
	       posprev = Find_EqPos( Sub[i-1].JointID, numEntries, Info);
	       SubEntries[i].PrevEqNum = Info[posprev].x;
	       SubEntries[NumJointSub/3 + i].PrevEqNum = Info[posprev].y;
	       SubEntries[NumJointSub*2/3 + i].PrevEqNum = Info[posprev].z;
	       Drum = Drum + 1;
	  }
	  SubEntries[i].Column = Column;
	  SubEntries[i].Drum = Drum;
	  SubEntries[NumJointSub/3 + i].Column = Column;
	  SubEntries[NumJointSub/3 + i].Drum = Drum;
	  SubEntries[NumJointSub*2/3 + i].Column = Column;
	  SubEntries[NumJointSub*2/3 + i].Drum = Drum;
     }	  
}

int Find_EqPos ( const int JointID, const int numEntries, const JointInfo_t const* Info )
{
     int i;
     bool found = false;

     i = 0;
     while (i < numEntries && !found ){
	  if (Info[i].JointID == JointID){
	       found = true;
	  } else {
	       i = i + 1;
	  }
     }

     if (!found){
	  printf("Could not find Joint ID %d in the TXE file.\n", JointID);
	  exit( EXIT_FAILURE );
     }

     return i;	 
}


void Reorder_JointEntries ( const int numEntries, JointInfo_t *const JointEntries )
{
     int j,bound = numEntries - 1;
     JointInfo_t temp;
     int i;

     for( i = 0 ; i < numEntries - 2 ; i++ ){
	 for ( j = 0 ; j < bound ; j++ ) {
	      if ( JointEntries[j].JointID > JointEntries[j + 1].JointID ) {
		   temp = JointEntries[j];
		   JointEntries[j] = JointEntries[j + 1];
		   JointEntries[j + 1] = temp;
	      }            
	 }
	 bound = j-1;
     }

     for (i = 0; i < numEntries; i++){
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

void Print_CoupleNodes ( const int NumJointSub, const int NumDrums, const SubEntry_t *const SubEntries, FILE *TheFile )
{
     
     int i;

     fprintf( TheFile, "%d\n", NumJointSub );
     
     for (i = 0; i < NumJointSub; i++ ){
	  switch  (SubEntries[i].Drum){
	  case 1:
	       fprintf(TheFile, "Sim_StoneDrums, 1 %d, 12 %d 1.01 1.0 0.5 0.5 2.0 1.0 0.0 1.0 0.0 1.0 0.0, Column: %d, Drum: %d;\n", SubEntries[i].EqNum, SubEntries[i].PrevEqNum, SubEntries[i].Column, SubEntries[i].Drum);
	       break;
	  case 2:
	       fprintf(TheFile, "Sim_StoneDrums, 1 %d, 12 %d 1.01 1.0 0.5 0.5 2.0 1.0 0.0 1.0 0.0 1.0 0.0, Column: %d, Drum: %d;\n", SubEntries[i].EqNum, SubEntries[i].PrevEqNum, SubEntries[i].Column, SubEntries[i].Drum);
	       break;
	  case 3:
	       fprintf(TheFile, "Sim_StoneDrums, 1 %d, 12 %d 1.01 1.0 0.5 0.5 2.0 1.0 0.0 1.0 0.0 1.0 0.0, Column: %d, Drum: %d;\n", SubEntries[i].EqNum, SubEntries[i].PrevEqNum, SubEntries[i].Column, SubEntries[i].Drum);
	       break;
	  case 4:
	       fprintf(TheFile, "Sim_StoneDrums, 1 %d, 12 %d 1.01 1.0 0.5 0.5 2.0 1.0 0.0 1.0 0.0 1.0 0.0, Column: %d, Drum: %d;\n", SubEntries[i].EqNum, SubEntries[i].PrevEqNum, SubEntries[i].Column, SubEntries[i].Drum);
	       break;
	  case 5:
	       fprintf(TheFile, "Sim_StoneDrums, 1 %d, 12 %d 1.01 1.0 0.5 0.5 2.0 1.0 0.0 1.0 0.0 1.0 0.0, Column: %d, Drum: %d;\n", SubEntries[i].EqNum, SubEntries[i].PrevEqNum, SubEntries[i].Column, SubEntries[i].Drum);
	       break;
	  case 6:
	       fprintf(TheFile, "Sim_StoneDrums, 1 %d, 12 %d 1.01 1.0 0.5 0.5 2.0 1.0 0.0 1.0 0.0 1.0 0.0, Column: %d, Drum: %d;\n", SubEntries[i].EqNum, SubEntries[i].PrevEqNum, SubEntries[i].Column, SubEntries[i].Drum);
	       break;
	  case 7:
	       fprintf(TheFile, "Sim_StoneDrums, 1 %d, 12 %d 1.01 1.0 0.5 0.5 2.0 1.0 0.0 1.0 0.0 1.0 0.0, Column: %d, Drum: %d;\n", SubEntries[i].EqNum, SubEntries[i].PrevEqNum, SubEntries[i].Column, SubEntries[i].Drum);
	       break;
	  default:
	       printf("Drum %d not recognised. It should be between 1 and 7.\n", SubEntries[i].Drum);
	       exit(EXIT_FAILURE);
	  }
     }
};

void print_help( char **argv, const struct option *long_opts )
{
     int i;
     struct option *Current_option = NULL;
     const char *Explanation[5] = {"Input file with joint coordinates",
				   "Output file with substructure entries",
				   "Number of entries",
				   "TXE file from SAP2000",
				   "This help text"};

     Current_option= (struct option *) &long_opts[0];

     printf("Usage: %s [-h] -i <InputFile> -o <OutputFile> -d <desiredDOF> -n <Num_Eqs>", argv[0] );
     printf("\n");
    

     printf("Summary of arguments\n");
     i = 0;
     while( (Current_option != NULL) && (Current_option->name != NULL ) ){
	  printf("  -%c  --%s    \t %s.\n", (char) Current_option->val, Current_option->name, Explanation[i] );
	  Current_option = Current_option + 1;
	  i = i + 1;
     }
}
