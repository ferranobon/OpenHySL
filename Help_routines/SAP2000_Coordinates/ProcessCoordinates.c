#include <stdio.h>
#include <stdlib.h>

int main (int argc, char **argv)
{
     FILE *TheFile;
     char *FileName;
     
     char *line = NULL;
     char dummy[40];
     size_t len = 0;
     ssize_t read;

     int JointID;
     float XCoord, YCoord, ZCoord;
     

     TheFile = fopen("Test.s2k", "r");
     if (TheFile == NULL){
	  exit(EXIT_FAILURE);
     }

     while ((read = getline(&line, &len, TheFile)) != -1){
	  sscanf( line, "   Joint=%d   %s   %s   %s   %s   %s   %s   GlobalX=%f   GlobalY=%f   GlobalZ=%f", &JointID, dummy, dummy, dummy, dummy, dummy, dummy, &XCoord, &YCoord, &ZCoord );
	  printf( "Joint %d is in Coord (%f,%f,%f)\n", JointID, XCoord, YCoord, ZCoord );
     }

     free(line);

     return 0;
}
