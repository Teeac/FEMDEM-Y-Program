/* Copyright (C) 2000, Dr. Antonio Munjiza
 *
 * This code is provided as part of the book entitled "The Combined
 * Finite Discrete Element Method". It is distributed WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
 * PARTICULAR PURPOSE. Inclusion of a part or whole of this code into any other
 * commercial or research or other purpose code is not granted without author's
 * written explicit permission.
 * When results using whole or any part of this code
 * are published, Y code must be mentioned and acknowledgement to Dr Munjiza must be made.
 * Should you modify this source code, the Copyright (C) on the mdified code
 * as a whole belongs to Dr. A. Munjiza regardless of the extent or nature
 * of modifications.
 * Copyright (C) to whole of any code containing any part of this code
 * also belongs to Dr. A.Munjiza. 
 * Any code comprising any part of this source code
 * must be called Y program.
 * If you do not agree with this, you are not allowed to do 
 * any modifications to any part of this source code or included 
 * any part of it in any other program.
 */
#define NeedFunctionPrototypes 1
/****************************************************************************
*                   frame.h
*
*  This header file is included by all C modules. It defines all
*  globally-accessible types and constants.
*
*  from A.Munjiza
*  Copyright 1993 SMILE
*---------------------------------------------------------------------------
*  NOTICE: This source code file is provided so that users may experiment
*  with SMILE under  
*****************************************************************************/
#ifndef FRAMEINCLUDED
#define FRAMEINCLUDED
/* Generic header for all modules */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <time.h>
#include <signal.h>
 
#define  Y_VERSION "1.0"

#define FLT float
#define INS int
#define DBL double
#define INT long int
#define CHR char
#define UINT unsigned int
#define ULON unsigned long

#define SETLINEBUF(fcheck) setvbuf((fcheck), NULL, _IONBF, 0);
/* #define SETLINEBUF(fcheck) setlinebuf((fcheck));		SGI */

#define FILENULL ((FILE*)NULL)
#define CHR2NULL ((CHR **)NULL)
#define DBL1NULL ((DBL*)NULL)
#define DBL2NULL ((DBL**)NULL)
#define DBL3NULL ((DBL***)NULL)
#define INT1NULL ((INT*)NULL)
#define INT2NULL ((INT**)NULL)
#define INT3NULL ((INT***)NULL)

static  DBL DBLTMP=0.0;
#define TRUE 1
#define FALSE 0
#define MYPI 3.1415926535897932384626
#define EPSILON 1.0e-10
#define NEPSILON -1.0e-10
#define BEPSILON 1.0e+15
#define R0 0.0
#define R1 1.0
#define R2 2.0
#define R3 3.0
#define R4 4.0
#define R5 5.0
#define R6 6.0
#define R7 7.0
#define R8 8.0
#define R9 9.0
#define R10 10.0
#define R12 12.0
#define R18 18.0
#define R19 19.0
#define R20 20.0
#define R24 24.0     
#define R36 36.0                 
#define alpha 0.585410196624969  
#define beta 0.138196601125011   
#define R21 21.0
#define RP8 0.8
#define RP2 0.2
#define RP11 1.1
#define RP15 1.5
#define RP25 0.25
#define RP5  0.5
#define RP31  3.1
#define RP9  0.9
#define RP1  0.1
#define RP3  0.3
#define RP7  0.7
#define RP75  0.75
#define RP37  0.37

#define prArbit 10000 /* Arbitrary number to add to the value of pr set assigned to elements */
#define RECL 4

#define ANTEATOI atoi
#define ANTEATOF atof
#define EXP exp
#define ACOS acos
#define ASIN asin
#define ATAN asin
#define SQRT sqrt
#define COS cos
#define SIN sin
#define FSQRT fsqrt
#define FCOS fcos
#define FSIN fsin
#define ABS(x) (((x)<0)?-(x):(x))
#define DABS(x) (((x)<R0)?-(x):(x))
#define MAXIM(x,y) (((x)<(y))?(y):(x))
#define MINIM(x,y) (((x)>(y))?(y):(x))
#define FREE(x) if ((x)!=NULL) free(x)
#define MALLOCED(x) (x)
#define NOTMALLOCED(x) (!x)
#define MALLOC(x) malloc(x)
 
/* =========Small MATRIX======================================== 
****************************************************************************
*  This module contains macros to perform operations on  small matrices
*****************************************************************************/
   /* INVERSE A SMALL MATRIX - also return determinant  */
#define YMATINV3(m,minv,det)\
      {  det=m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])-\
             m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])+\
             m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]);\
         minv[0][0]=(m[1][1]*m[2][2]-m[1][2]*m[2][1])/det;\
         minv[1][0]=(m[1][2]*m[2][0]-m[1][0]*m[2][2])/det;\
         minv[2][0]=(m[1][0]*m[2][1]-m[1][1]*m[2][0])/det;\
         minv[0][1]=(m[0][2]*m[2][1]-m[0][1]*m[2][2])/det;\
         minv[1][1]=(m[0][0]*m[2][2]-m[0][2]*m[2][0])/det;\
         minv[2][1]=(m[0][1]*m[2][0]-m[0][0]*m[2][1])/det;\
         minv[0][2]=(m[0][1]*m[1][2]-m[0][2]*m[1][1])/det;\
         minv[1][2]=(m[0][2]*m[1][0]-m[0][0]*m[1][2])/det;\
         minv[2][2]=(m[0][0]*m[1][1]-m[0][1]*m[1][0])/det;\
    }

#define YMATDET3(m,det)\
      {  det=m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])-\
             m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])+\
             m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]);\
    }
    
#define YMATINV2(m,minv,det)\
      {  det=m[0][0]*m[1][1]-m[1][0]*m[0][1]; \
         minv[0][0]= m[1][1]/det; minv[1][0]=-m[1][0]/det;\
         minv[0][1]=-m[0][1]/det; minv[1][1]= m[0][0]/det;\
    }
/* =========vect3D======================================== 
****************************************************************************
*  This module contains macros to perform operations on 3D vectors.
*****************************************************************************/

/* Misc. Vector Math Macro Definitions */

extern DBL VTemp;
/* Assign; Null */
#define V3DAss(x1,y1,z1,x2,y2,z2) {(x1)=(x2);(y1)=(y2);z1=(z2);}
#define V3DNul(x1,y1,z1) {(x1)=R0;(y1)=R0;(z1)=R0;}
/* Vector Add; Subtract; Scale; Length square; Length; 
   Dot product; Cross product; Volume of triad */
#define V3DAdd(x1,y1,z1,x2,y2,z2,x3,y3,z3) \
              {(x1)=(x2)+(x3);(y1)=(y2)+(y3);(z1)=(z2)+(z3);}
#define V3DSub(x1,y1,z1,x2,y2,z2,x3,y3,z3) \
              {(x1)=(x2)-(x3);(y1)=(y2)-(y3);(z1)=(z2)-(z3);}
#define V3DSca(x1,y1,z1,s){(x1)=(x1)*(s);(y1)=(y1)*(s);(z1)=(z1)*(s);}
#define V3DDiv(x1,y1,z1,s){(x1)=(x1)/(s);(y1)=(y1)/(s);(z1)=(z1)/(s);}
#define V3DLe2(s,x1,y1,z1){(s)=(x1)*(x1)+(y1)*(y1)+(z1)*(z1);}
#define V3DLen(s,x1,y1,z1){(s)=SQRT((x1)*(x1)+(y1)*(y1)+(z1)*(z1));}
#define V3DDot(s,x1,y1,z1,x2,y2,z2) {(s)=((x1)*(x2))+((y1)*(y2))+((z1)*(z2));}
#define V3DCro(x1,y1,z1,x2,y2,z2,x3,y3,z3)\
              {(x1)=(y2)*(z3)-(z2)*(y3); \
               (y1)=(z2)*(x3)-(x2)*(z3); \
               (z1)=(x2)*(y3)-(y2)*(x3); }
#define V3DOmegaRot(x1,y1,z1,x2,y2,z2,x3,y3,z3)\
              {(x1)=(x3)+(y2)*(z3)-(z2)*(y3); \
               (y1)=(y3)+(z2)*(x3)-(x2)*(z3); \
               (z1)=(z3)+(x2)*(y3)-(y2)*(x3); }
#define V3DTranToGl(x1,y1,z1,a,b,c,x2,y2,z2,x3,y3,z3)\
              {(x1)=(a)*(x2)+(b)*(x3)+(c)*((y2)*(z3)-(z2)*(y3)); \
               (y1)=(a)*(y2)+(b)*(y3)+(c)*((z2)*(x3)-(x2)*(z3)); \
               (z1)=(a)*(z2)+(b)*(z3)+(c)*((x2)*(y3)-(y2)*(x3)); }
#define V3DTranToLoc(x1,y1,z1,a,b,c,x2,y2,z2,x3,y3,z3)\
              {(x1)=(a)*(x2)+(b)*(y2)+(c)*(z2);  \
               (y1)=(a)*(x3)+(b)*(y3)+(c)*(z3); \
(z1)=(a)*((y2)*(z3)-(z2)*(y3))+(b)*((x3)*(z2)-(x2)*(z3))+(c)*((x2)*(y3)-(y2)*(x3)); }
#define V3DVol(v,x2,y2,z2,x3,y3,z3,x1,y1,z1) \
              {(v)=((y2)*(z3)-(z2)*(y3))*(x1)+ \
                   ((z2)*(x3)-(x2)*(z3))*(y1)+ \
                   ((x2)*(y3)-(y2)*(x3))*(z1); }
  /*  Normalize a Vector; Halfway Vectors; */
#define V3DNor(s,x1,y1,z1) \
              {(s)=SQRT((x1)*(x1)+(y1)*(y1)+(z1)*(z1));   \
              if((s)>EPSILON)(x1)=(x1)/(s);  \
              if((s)>EPSILON)(y1)=(y1)/(s);  \
              if((s)>EPSILON)(z1)=(z1)/(s);  } 
#define V3DMid(x1,y1,z1,x2,y2,z2,x3,y3,z3) \                                   \
              {(x1)=(RP5)*((x2)+(x3)); \
               (y1)=(RP5)*((y2)+(y3)); \
               (z1)=(RP5)*((z2)+(z3)); }
/* Vector in Between two given vectors  */
#define V3DBet(x1,y1,z1,x2,y2,z2,x3,y3,z3,s2,s3)\
              {(x1)=(x2)*(s2)+(x3)*(s3);     \
               (y1)=(y2)*(s2)+(y3)*(s3);     \
               (z1)=(z2)*(s2)+(z3)*(s3);     }
/* =========vect2D======================================== 
****************************************************************************
*  This module contains macros to perform operations on 2D vectors.
*****************************************************************************/
/* Assign; Null */
#define V2DAss(x1,y1,x2,y2) {(x1)=(x2);(y1)=(y2)}
#define V2DNul(x1,y1) {(x1)=R0;(y1)=R0;} 
/* Vector Add; Subtract; Scale; Length square; Length; Dot product; Cross product;  */
#define V2DAdd(x1,y1,x2,y2,x3,y3) {(x1)=(x2)+(x3);(y1)=(y2)+(y3);}
#define V2DSub(x1,y1,x2,y2,x3,y3) {(x1)=(x2)-(x3);(y1)=(y2)-(y3);}
#define V2DSca(x1,y1,s) {(x1)=(x1)*s;(y1)=(y1)*s;}
#define V2DLe2(s,x1,y1,z1) {(s)=(x1)*(x1)+(y1)*(y1);}
#define V2DLen(s,x1,y1,z1) {(s)=SQRT((x1)*(x1)+(y1)*(y1);}
#define V2DDot(s,x1,y1,x2,y2,z2) {s=(x1*x2)+(y1*y2);}
#define V2DCro(s,x1,y1,x2,y2){s=(x1*y2)-(y1*x2);} 
  
/*  Normalize a Vector; Compute a Vector Halfway Between Two Given Vectors; */ 
#define V2DNor(x1,y1,x2,y2){VTemp=SQRT((x1)*(x1)+(y1)*(y1));   \
              (x1)=(x2)/VTemp;(y1)=(y2)/VTemp;                 } 
#define V2DMid(x1,y1,x2,y2,x3,y3){(x1)=0.5*((x2)+(x3));(y1)=0.5*((y2)+(y3));}
                                               
/* Vector in Between two given vectors  */
#define V2DBet(x1,y1,x2,y2,x3,y3,s2,s3){(x1)=(x2)*(s2)+(x3)*(s3);     \
              (y1)=(y2)*(s2)+(y3)*(s3);     }
/* =========STREAM======================================== 
****************************************************************************
*  This module contains macros to perform operations on STREAMS
*****************************************************************************/
#define FILEND(fptr) (feof(fptr))
static CHR *DBL_S[20]=
{ "%le"     ,"%+1.0le" ,"%+2.0le" ,"%+3.0le" ,"%+4.0le" ,"%+5.0le" ,
  "%+6.0le" ,"%+7.0le" ,"%+8.1le" ,"%+9.2le" ,"%+10.3le",
  "%+11.4le","%+12.5le","%+13.6le","%+14.7le","%+15.8le",
  "%+16.9le","%+17.10le","%+18.11le","%+19.12le" 
};
#define DBL_SR "%le" 
static CHR *INT_S[20]=
{"%ld","%1ld","%2ld","%3ld","%4ld","%5ld","%6ld","%7ld","%8ld","%9ld",
 "%10ld","%11ld","%12ld","%13ld","%14ld",
 "%15ld","%16ld" ,"%17ld","%18ld","%19ld"
};
#define INT_SR   "%ld"

static CHR *CHR_S[30]=
{"%s", "%1s", "%2s", "%3s", "%4s", "%5s", "%6s", "%7s", "%8s", "%9s", "%10s",
 "%11s", "%12s", "%13s", "%14s", "%15s", "%16s", "%17s", "%18s", "%19s", "%20s",
 "%21s", "%22s", "%23s", "%24s", "%25s", "%26s", "%27s", "%28s", "%29s"
}; 
#define CHR_SR "%s"
#define CHRRETURN "\n"
#define CHRSPACE " "
#define CHRTERMINATE 0
#define CHRTAB "\t"
#define INTr(file,x){fscanf((file),INT_SR,((x)));}
#define DBLr(file,x){fscanf((file),DBL_SR,((x)));}
#define CHRr(file,x){fscanf((file),CHR_SR,((x)));}
#define INTw(file,x,ndigit)fprintf((file),(INT_S[(ndigit)]),((x))); 
#define SINTw(file,x,ndigit)sprintf((file),(INT_S[(ndigit)]),((x))); 
#define DBLw(file,x,ndigit) fprintf((file),(DBL_S[((ndigit))]),(x));
#define FLTw(file,x,ndigit) fprintf((file),(FLT_S[((ndigit))]),(x));
 
#define SCHRw(file,x){sprintf((file),CHR_SR,(x));}
#define CHRw(file,x){fprintf((file),CHR_SR,(x));}
#define CHRww(file, x, ndigit){fprintf((file),(CHR_S[(ndigit)]),(x));};
#define CHRend(c,i){c[i]='\0';}
#define CHRwcr(file){fprintf((file),CHR_SR,(CHRRETURN));}
#define CHRwsp(file){fprintf((file),CHR_SR,(CHRSPACE));}
#define CHRwtab(file){fprintf((file),CHR_SR,(CHRTAB));}
#define CHRcmp(c1,c2,n)(strncmp((c1),(c2),(n)))
#define CHRcpy(c1,c2)(strcpy((c1),(c2)))
#define CHRcat(c1,c2)(strcat((c1),(c2)))
void CHRcpynoext(  /* copy no extension */
#if NeedFunctionPrototypes 
  CHR *c1,              /* to copy to    */
  CHR *c2               /* to copy from  */
#endif
); 
void CHRcpynopath( /* copy no extension and path */
#if NeedFunctionPrototypes 
  CHR *c1,              /* to copy to    */
  CHR *c2               /* to copy from  */
#endif
); 

DBL *TalDBL1(  /* alocate array */
#if NeedFunctionPrototypes 
  INT m1               /* size          */
#endif
); 
DBL *TalDBL1I(  /* alocate array */
#if NeedFunctionPrototypes 
  INT m1,              /* size          */
  DBL dini             /* initial value */
#endif
); 

DBL **TalDBL2(  /* alocate array */
#if NeedFunctionPrototypes 
  INT m2,              /* m2 1D arrays size          */
  INT m1               /* of size   m1 DBL           */ 
#endif
); 
DBL **TalDBL2I(  /* alocate array */
#if NeedFunctionPrototypes 
  INT m2,              /* m2 1D arrays size          */
  INT m1,              /* of size   m1 DBL           */ 
  DBL dini             /* initial value              */
#endif
); 
DBL ***TalDBL3(  /* alocate array */
#if NeedFunctionPrototypes 
  INT m3,              /* m3 2D arrays size m2*m1         */
  INT m2,              /* m3*m2 1D arrays                 */
  INT m1               /* of size   m1 DBL                */ 
#endif
);
DBL ***TalDBL3I(  /* alocate array */
#if NeedFunctionPrototypes 
  INT m3,              /* m3 2D arrays size m2*m1         */
  INT m2,              /* m3*m2 1D arrays                 */
  INT m1,              /* of size   m1 DBL                */ 
  DBL dini             /* initial value                   */
#endif
);
INT *TalINT1(  /* alocate array */
#if NeedFunctionPrototypes 
  INT m1               /* size          */
#endif
); 
INT *TalINT1I(  /* alocate array with initial value */
#if NeedFunctionPrototypes 
  INT m1, INT iini               /* size & initial value  */
#endif
); 
INT **TalINT2(  /* alocate array */
#if NeedFunctionPrototypes 
  INT m2,              /* m2 1D arrays size          */
  INT m1               /* of size   m1 INT           */
#endif
); 
INT **TalINT2I(  /* alocate array with initial value  */
#if NeedFunctionPrototypes 
  INT m2,              /* m2 1D arrays size          */
  INT m1,              /* of size   m1 INT           */
  INT iini             /* initial value              */
#endif
); 
INT ***TalINT3(  /* alocate array */
#if NeedFunctionPrototypes 
  INT m3,              /* m3 2D arrays size m2*m1         */
  INT m2,              /* m3*m2 1D arrays                 */
  INT m1               /* of size   m1 INT                */ 
#endif
);

/**********--------------------------------------*******/

INT Getname(  /* return i for argv[i] that follows after "name" like -i file */
#if NeedFunctionPrototypes 
  INT argc,          /* number of argument in argv      */
  CHR **argv,        /* arguments                        */
  CHR *name          /* prefix to be looked for like -i  */
#endif
);
void TformDBL1(  /* read array */
#if NeedFunctionPrototypes 
  FILE *fptr,          /* file          */
  DBL dinit,           /* initialise to */
  INT n1,              /* size          */
  DBL **d1aray          /* array pointer */
#endif
);

void TformDBL2(  /* read array */
#if NeedFunctionPrototypes 
  FILE *fptr,          /* file          */
  DBL dinit,           /* initialise to */
  INT n1,              /* size          */
  INT n2,              /* size          */
  DBL ***d2aray         /* array pointer */
#endif
);
void TformDBL3(  /* read array */
#if NeedFunctionPrototypes 
  FILE *fptr,          /* file          */
  DBL dinit,           /* initialise to */
  INT n1,              /* size          */
  INT n2,              /* size          */
  INT n3,              /* size          */
  DBL ****d3aray         /* array pointer */
#endif
);
void TformINT1(  /* read array */
#if NeedFunctionPrototypes 
  FILE *fptr,          /* file          */
  INT iinit,           /* initialise to */
  INT n1,              /* size          */
  INT **i1aray          /* array pointer */
#endif
);

void TformINT2(  /* read array */
#if NeedFunctionPrototypes 
  FILE *fptr,          /* file          */
  INT iinit,           /* initialise to */
  INT n1,              /* size          */
  INT n2,              /* size          */
  INT ***i2aray         /* array pointer */
#endif
);
void TformINT3(  /* read array */
#if NeedFunctionPrototypes 
  FILE *fptr,          /* file          */
  INT iinit,           /* initialise to */
  INT n1,              /* size          */
  INT n2,              /* size          */
  INT n3,              /* size          */
  INT ****i3aray         /* array pointer */
#endif
);

void TreadDBL1(  /* read array */
#if NeedFunctionPrototypes 
  FILE *fptr,          /* file          */
  INT n1,              /* size          */
  DBL *d1aray          /* array pointer */
#endif
);

void TreadDBL2(  /* read array */
#if NeedFunctionPrototypes 
  FILE *fptr,          /* file          */
  INT n1,              /* size          */
  INT n2,              /* size          */
  DBL **d2aray         /* array pointer */
#endif
);

void TreadINT1(  /* read array */
#if NeedFunctionPrototypes 
  FILE *fptr,          /* file          */
  INT n1,              /* size          */
  INT *i1aray          /* array pointer */
#endif
);

void TreadINT2(  /* read array */
#if NeedFunctionPrototypes 
  FILE *fptr,          /* file          */
  INT n1,              /* size          */
  INT n2,              /* size          */
  INT **i2aray         /* array pointer */
#endif
);
/**********--------------------------------------*******/
FILE *TwriteCHR2(  /* write array */
#if NeedFunctionPrototypes 
  CHR *namep,          /* file name                 */
  FILE *fpt,           /* file                      */
  INT ndigit,          /*number of digits  ignored  */
  INT nperli,          /* numbers per line          */
  INT n1,              /* size                      */
  CHR **c2aray,        /* array pointer             */
  CHR *name            /* array name                */
#endif
);

FILE *TwriteDBL1(  /* write array */
#if NeedFunctionPrototypes 
  CHR *namep,          /* file name       */
  FILE *fpt,           /* file             */
  INT ndigit,          /*number of digits  */
  INT nperli,          /* numbers per line */
  INT n1,              /* size             */
  DBL *d1aray,         /* array pointer    */
  CHR *name            /* array name       */
#endif
);

FILE *TwriteDBL2(  /* write array */
#if NeedFunctionPrototypes 
  CHR  *namep,         /* file name     */
  FILE *fpt,           /* file             */
  INT ndigit,          /*number of digits  */
  INT nperli,          /* numbers per line */
  INT ihow,            /* which order      */
  INT n1,              /* size             */
  INT n2,              /* size             */
  DBL **d2aray,        /* array pointer    */
  CHR *name            /* array name       */
#endif
);
FILE *TwriteDBL3(  /* write array */
#if NeedFunctionPrototypes 
  CHR  *namep,         /* file name     */
  FILE *fpt,           /* file             */
  INT ndigit,          /*number of digits  */
  INT nperli,          /* numbers per line */
  INT ihow,            /* which order      */
  INT n1,              /* size             */
  INT n2,              /* size             */
  INT n3,              /* size             */
  DBL ***d3aray,        /* array pointer    */
  CHR *name            /* array name       */
#endif
);
 
FILE *TwriteINT1(  /* write array */
#if NeedFunctionPrototypes 
  CHR  *namep,         /* file name     */
  FILE *fpt,           /* file             */
  INT ndigit,          /*number of digits  */
  INT nperli,          /* numbers per line */
  INT n1,              /* size             */
  INT *i1aray,         /* array pointer    */
  CHR *name            /* array name       */
#endif
);

FILE *TwriteINT2(  /* write array */
#if NeedFunctionPrototypes 
  CHR  *namep,         /* file name     */
  FILE *fpt,           /* file             */
  INT ndigit,          /*number of digits  */
  INT nperli,          /* numbers per line */
  INT ihow,            /* numbers per line */
  INT n1,              /* size             */
  INT n2,              /* size             */
  INT **i2aray,        /* array pointer    */
  CHR *name            /* array name       */
#endif
);
FILE *TwriteINT3(  /* write array */
#if NeedFunctionPrototypes 
  CHR  *namep,         /* file name     */
  FILE *fpt,           /* file             */
  INT ndigit,          /*number of digits  */
  INT nperli,          /* numbers per line */
  INT ihow,            /* numbers per line */
  INT n1,              /* size             */
  INT n2,              /* size             */
  INT n3,              /* size             */
  INT ***i3aray,        /* array pointer    */
  CHR *name            /* array name       */
#endif
);
/*********SORTING ARRAYS**********************/
void TsortDBL1(  /* sort array = smallest ... largest*/
#if NeedFunctionPrototypes 
  INT n,               /* size                                  */
  DBL *d1x,            /* array pointer (array to be sorted)    */
  INT *i1y             /* array pointer (array to be rearanged) */
#endif
);
void TsortINT(   /* sort array = smallest ... largest*/
#if NeedFunctionPrototypes 
  INT n,               /* size                                  */
  INT nsort,           /* total number of arrays to be sorted   */
  INT nrear,           /* total number of arrays                */
  INT **i2             /* array pointers (sorted or rearanged)  */
#endif
); 
/*************SPACESAVING by 90 FORMAT ************/
void readcode(  /* read coded array */
#if NeedFunctionPrototypes 
  CHR *c1code,         /* coded array [0]=ndigit; [1]=nnum     */
  INT *i1num           /* INT array [0]=ndigit; [1]=nnum
                          the rest between 0 and base^ndigit    */
#endif
);
void codeCHRtoINT(  /* write array as coded */
#if NeedFunctionPrototypes 
  CHR *c1code,         /* coded array [0]=ndigit; [1]=nnum     */
  INT *i1num           /* INT array [0]=ndigit; [1]=nnum
                          the rest between 0 and base^ndigit    */
#endif
);
void codeINTtoCHR(  /* write array as coded */
#if NeedFunctionPrototypes 
  CHR *c1code,         /* coded array [0]=ndigit; [1]=nnum     */
  INT *i1num           /* INT array [0]=ndigit; [1]=nnum
                          the rest between 0 and base^ndigit    */
#endif
);

void codeINTtoDBL(  /* translate INT to DBL array */
#if NeedFunctionPrototypes 
  DBL *d1num,          /* DBL array [0]=ndigit; [1]=nnum
                          all between -1 and +1                */  
  INT *i1num           /* INT array [0]=ndigit; [1]=nnum
                          the rest between 0 and base^ndigit    */
#endif
);

void codeDBLtoINT(  /* translate DBL to INT array */
#if NeedFunctionPrototypes 
  DBL *d1num,          /* DBL array [0]=ndigit; [1]=nnum
                          all between -1 and +1                */  
  INT *i1num           /* INT array [0]=ndigit; [1]=nnum
                          the rest between 0 and base^ndigit    */
#endif
);
 
void UnVec(         /* calc. unit vector from vector r  */
#if NeedFunctionPrototypes 
  DBL *e1x,         /* x member of unit vector          */
  DBL *e1y,         /* y member of unit vector          */
  DBL rx,           /* x member of vector r             */
  DBL ry            /* y member of vector r             */
#endif
);

void Rot90(         /* rotate by -90 degrees            */
#if NeedFunctionPrototypes
  DBL *e1x,         /* x member of rotated unit vector  */
  DBL *e1y,         /* y member of rotated unit vector  */
  DBL e2x,          /* x member of unit vector to rotate*/
  DBL e2y           /* y member of unit vector to rotate*/
#endif
);

void V2DTranToLoc(  /* transform node to local coord.   */
#if NeedFunctionPrototypes
  DBL *u,           /* 1st local coordinate (=>x)       */
  DBL *v,           /* 2nd local coordinate (=>y)       */
  DBL rx,           /* x coord. of the node             */
  DBL ry,           /* y coord. of the node             */
  DBL e1x,          /* x member of 1st unit vector      */
  DBL e1y,          /* y member of 1st unit vector      */
  DBL e2x,          /* x member of 2nd unit vector      */
  DBL e2y           /* y member of 2nd unit vector      */
#endif
);

/*********************************************************/
void Push_into_INT1( /* push an integer record into an array */
#if NeedFunctionPrototypes 
  INT *arr,          /* INT array to be processed            */
  INT ni,            /* position to be pushed                */
  INT nj             /* INT record to be pushed              */
#endif
);

INT Get_from_INT1(  /* get records from the array one by one */
#if NeedFunctionPrototypes 
  INT *arr,         /* INT array to be processed             */
  INT ni            /* header position                       */
#endif
);

/********************* TIME FUNCTION *********************/
/* Returns character data on the real-time clock and date */
void date_and_time( 
#if NeedFunctionPrototypes 
char *date_time
#endif
);

/* Prints the current YMDHMS date as a time stamp */
void timestamp(void);

/* Returns the number of seconds elapsed since the program was launched. */
double run_time(void);

/********************* CTRL-C Pause *********************/
int enablc(
#if NeedFunctionPrototypes
DBL dctime, INT ncstep, INT mcstep
#endif
);

/*********************** FILE POSITION ******************/
INT nline_current_position(
#if NeedFunctionPrototypes
FILE *fptr
#endif
);
#endif
