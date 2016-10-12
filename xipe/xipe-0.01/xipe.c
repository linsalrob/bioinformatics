/* $Revision: 1.1.1.1 $ $Id: xipe.c,v 1.1.1.1 2007/09/14 22:49:31 linsalrob Exp $ $Date: 2007/09/14 22:49:31 $*/

/* 
#Copyright (C) 2005 Beltran Rodriguez-Brito,
#Pat MacNarnie, and Rober Edwards
#
#This program is free software; you can redistribute it and/or
#modify it under the terms of the GNU General Public License
#as published by the Free Software Foundation; either version 2
#of the License, or (at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program; if not, write to the Free Software
#Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#
#Released as part of the xipe package as supplemental online material
#to "A statistical approach for metagenome comparisons"
# by Beltran Rodriguez-Brito, Forest Rohwer, and Robert A. Edwards.
*/

#include<stdio.h>
#include<stdlib.h>

main(int argc, char** argv)
{
    /* Set defaults for all parameters: */

    /* const int MAXSUBSYS = 350;     this was the old line */
    const int MAXSUBSYS = 32000;
    FILE *infile_1; FILE *infile_2; 
    FILE *outfile_1; FILE *outfile_2;
    int repeats = 100;        int samplesize = 200;
    char* infile1 = NULL;     char* outfile1 = NULL;
    char* infile2 = NULL;     char* outfile2 = NULL;
    float colu1; int colu2; int colu3;
    int i;
    int myerror;
    int numlines1 = 0; int numlines2 = 0;
    int mydata1[MAXSUBSYS]; int myresu1[MAXSUBSYS];
    int mydata2[MAXSUBSYS]; int myresu2[MAXSUBSYS];
    int myresu3[MAXSUBSYS];
    /* clean the arrays :D */
    for (i = 1; i < MAXSUBSYS; i++) {
       mydata1[i] = 0; myresu1[i]=0; 
       mydata2[i] = 0; myresu2[i]=0; }

    /***********************z*********************/
    /* Command line parsing **********************/
    /* Start at i = 1 to skip the command name. */
    for (i = 1; i < argc; i++) {
	/* Check for a switch (leading "-"). */
	if (argv[i][0] == '-') {
	    /* Use the next character to decide what to do. */
	    switch (argv[i][1]) {
		case 'r':	repeats = atoi(argv[++i]); break;
		case 's':	samplesize = atoi(argv[++i]); break;
		case 'f':	infile1 = argv[++i];       break;
		case 'g':	infile2 = argv[++i];       break;
		case 'o':	outfile1 = argv[++i];       break;
		case 'p':	outfile2 = argv[++i];       break;
	    } } }
    /***********************z*********************/

    srand( samplesize + repeats );
    printf("Esteemed user, your parameters are:\n");
    printf("-r = repeats = %d\n", repeats);
    printf("-s = samplesize = %d\n", samplesize);
    if (infile1 != NULL)  printf("-f = input file 1 = \"%s\"\n", infile1);
    if (infile2 != NULL)  printf("-g = input file 2 = \"%s\"\n", infile2);
    if (outfile1 != NULL) printf("-o = medians file = \"%s\"\n", outfile1);
    if (outfile2 != NULL) printf("-p = ranges file  = \"%s\"\n", outfile2);
    if((infile_1 = fopen(infile1, "r")) == NULL) {
       printf(" EEEEEK!!!! \n You are a mean person! \n I can't open that file.\n");
       exit(1);                                }
    if((infile_2 = fopen(infile2, "r")) == NULL) {
       printf(" EEEEEK!!!! \n You are a mean person! \n I can't open that file.\n");
       exit(1);                                }

    numlines1 = 0;
    myerror = fscanf( infile_1, "%f %d %d\n", &colu1, &colu2, &colu3);
/*
    printf( " col 1        2   3   err \n");
    printf( " %f     %d   %d   %d \n", colu1, colu2, colu3, myerror );
*/
    if (myerror ==3) {
       mydata1[colu2]++;
       }
    while ( myerror == 3) {
    numlines1++ ;
    myerror = fscanf( infile_1, "%f %d %d\n", &colu1, &colu2, &colu3);
/*
    printf( " %f     %d   %d   %d \n", colu1, colu2, colu3, myerror );
*/
    if (myerror ==3) {
       mydata1[colu2]++;
       }
      }
    fclose(infile_1);
    printf( "number of lines of file 1 %d \n", numlines1) ;

    numlines2 = 0;
    myerror = fscanf( infile_2, "%f %d %d\n", &colu1, &colu2, &colu3);
/*
    printf( " %f  %d   %d   %d \n", colu1, colu2, colu3, myerror );
*/
    if (myerror ==3) {
       mydata2[colu2]++;
       }
    while ( myerror == 3) {
    numlines2++ ;
    myerror = fscanf( infile_2, "%f %d %d\n", &colu1, &colu2, &colu3);
/*
    printf( " %f  %d   %d   %d \n", colu1, colu2, colu3, myerror );
*/
    if (myerror ==3) {
       mydata2[colu2]++;
       }
      }
    fclose(infile_2);
    printf( "number of lines of file 2 %d \n", numlines2) ;

    /***********************z*********************/
    /* big array here */    {
    /***********************z*********************/
    unsigned int TheWhole[numlines1+1];
    unsigned int TheWhole2[numlines2+1];
/*
    unsigned short int TheWhole3[numlines1+numlines2+1]; */
    int MyCounter; int MyCounter2;
    int j;	int k;	int j1;	int j2;

    if((outfile_1 = fopen(outfile1, "w")) == NULL) {
       printf(" EEEEEK!!!! \n You are a mean person! \n I can't write that file.\n");
       exit(1);                                }
    for (k=1; k < repeats+1; k++) {

    MyCounter = 1;
    for (i= 1; i < MAXSUBSYS; i++) {
        for (j = 1; j < mydata1[i]+1; j++ ) {
        TheWhole[MyCounter] = i;
        MyCounter++;
        }   }
    for (i= 1; i < samplesize+1; i++) {
        j=1+(int) ((1.0*(MyCounter-1))*rand()/(RAND_MAX+1.0));
        myresu1[TheWhole[j]]++;
        }

    MyCounter2 = 1;
    for (i= 1; i < MAXSUBSYS; i++) {
        for (j = 1; j < mydata2[i]+1; j++ ) {
        TheWhole2[MyCounter2] = i;
        MyCounter2++;
        }   }
    for (i= 1; i < samplesize+1; i++) {
        j=1+(int) ((1.0*(MyCounter2-1))*rand()/(RAND_MAX+1.0));
        myresu2[TheWhole2[j]]++;
        }
    for (i= 1; i < MAXSUBSYS; i++) {
        myresu3[i]=myresu2[i]-myresu1[i];
        myresu2[i] = 0;
        myresu1[i] = 0;
        myerror = fprintf( outfile_1, "%d ", myresu3[i]);
        }
    myerror = fprintf( outfile_1, "\n" );

    }
    fclose(outfile_1);


    if((outfile_2 = fopen(outfile2, "w")) == NULL) {
       printf( \
       "EEEK!!!!\n You are a mean person!\n I can't write that file.\n");
       exit(1);                                }

    //printf(" repeats %d samplesize %d j1 %d j2 %d \n", \
    //         repeats, samplesize, j1, j2 );
    //for (k=1; k < MyCounter; k++) {
    //  printf ("i %d, TheWhole[i] %d\n", k, TheWhole[k]); }
    //for (k=1; k < MyCounter2; k++) {
    //  printf ("i %d, TheWhole2[i] %d\n", k, TheWhole2[k]); }
    for (k=1; k < repeats+1; k++) {

      j1 = 1 + (int) ((1.0*samplesize-1)*rand()/(RAND_MAX+1.0));
      j2 = samplesize - j1;
      for (i = 1; i < j1+1; i++) {
        j=1+(int) ((1.0*(MyCounter-1))*rand()/(RAND_MAX+1.0));
        myresu1[TheWhole[j]]++;
        }
      for (i = 1; i < j2+1; i++) {
        j=1+(int) ((1.0*(MyCounter2-1))*rand()/(RAND_MAX+1.0));
        myresu1[TheWhole2[j]]++;
        }

      j1 = 1 + (int) ((1.0*samplesize-1)*rand()/(RAND_MAX+1.0));
      j2 = samplesize - j1;
      //printf("j1 %d j2 %d MyCounter2 %d \n", j1, j2, MyCounter2);
      for (i = 1; i < j1+1; i++) {
        j=1+(int) ((1.0*MyCounter-1)*rand()/(RAND_MAX+1.0));
        //printf(" hola %d %d %d \n", i, j, TheWhole[j]);
        myresu2[TheWhole[j]]++;
        }
      for (i = 1; i < j2+1; i++) {
        j=1+(int) ((1.0*(MyCounter2-1))*rand()/(RAND_MAX+1.0));
        myresu2[TheWhole2[j]]++;
        }

/* replaced by code above
    for (i= 1; i < samplesize+1; i++) {
        j=1+(int) ((1.0*(MyCounter-1))*rand()/(RAND_MAX+1.0));
        myresu1[TheWhole3[j]]++;
        j=1+(int) ((1.0*(MyCounter-1))*rand()/(RAND_MAX+1.0));
        myresu2[TheWhole3[j]]++;
        }
*/
    for (i= 1; i < MAXSUBSYS; i++) {
        myresu3[i] = myresu2[i]-myresu1[i];
        myresu2[i] = 0;
        myresu1[i] = 0;
        myerror = fprintf( outfile_2, "%d ", myresu3[i]);
        }
    myerror = fprintf( outfile_2, "\n" );

    }
    /***********************z*********************/
    /* big array rip here */}
    /***********************z*********************/
}
