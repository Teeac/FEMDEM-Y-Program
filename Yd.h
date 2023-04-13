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
 * Should you modify this source code, the Copyright (C) on the modified code
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
/* file Yd.h  Y data base description */
#include "Ytypes.h"
#ifndef FRAMEINCL
#include "frame.h"
#define FRAMEINCL
#endif
#ifndef YTYPESINCL
#include "Ytypes.h"
#define YTYPESINCL
#endif

typedef struct YDC_struct *YDC;
struct YDC_struct
{ INT  mcstep, ncstep;       /* maximum/current number of time steps                  */  
  CHR  cruntime[256];        /* run time                                              */
  CHR  cfiname[256];         /* input file name                                       */
  FILE *finp, *fcheck;

  FILE *files[10];           /* input files' pointer                                  */
  INT nfile;                 /* ID of the current input file                          */
  
  DBL  dcgrax;               /* dcgray gravity x (deft: 0.0)                          */
  DBL  dcgray;               /* dcgray gravity y (deft: 0.0)                          */
  DBL  dcgraz;               /* dcgray gravity z (deft: 0.0)                          */
  DBL  dcsizc;               /* size coord.      (deft: 1.0)                          */
  DBL  dcsizf;               /* size force       (deft: 1.0)                          */
  DBL  dcsizs;               /* size stress      (deft: 1.0)                          */
  DBL  dcsizv;               /* size velocity    (deft: 1.0)                          */
  DBL  dcstec;               /* current time step size  (deft: 0.0)                   */
  DBL  dctime;               /* current time            (deft: 0.0)                   */
  DBL  dcurelx;              /* under relaxation for mass matrix (deft: 1.0)          */
  INT  initer;               /* number of interation for multi-pass algorithm(deft: 2)*/
  INT  icoutt;               /* output file type. 0-ym(deft.), 1-binary               */
  INT  icoutf;               /* write output frequency  (deft: 10)                    */
  INT  icouti;               /* current write output No (deft: 0)                     */
  INT  icoutp;               /* output precision - digits per number (deft: 4)        */
  INT  icshtf;               /* show time frequency. (deft: 100)                      */
  INT  iwfast;               /* fast mode. on: 0, off: >0 (deft: off)                 */
  INT  ietype;               /* element type - 4-nodes or 10-nodes (deft: -1)         */  
};

typedef struct YDE_struct *YDE;
struct YDE_struct
{ INT melem, nelem;          /* maximum (actual) number of elements                   */
  INT melst, nelst;          /* maximum (actual) number of elemen. states var.        */
  INT melno, nelno;          /* maximum (actual) number of elemen. nodes              */

  INT   *i1elcf;             /*[melem]    contacting couple first                     */
  INT   *i1elpr;             /*[melem]    element property                            */
  INT   *i1elprtmp;          /*[melem]    element property  used for meshing          */
  DBL   *d1elfs;             /*[melem]    shear strength at joint (Mohr-Coulomb)      */
  DBL   **d2elst;            /*[melst][melem]    - element state                      */
  INT   **i2elto;            /*[melno][melem]    - element topology                   */
  DBL   ***d3tcs;            /*[ndime][ndime][melem] -Cauchy stress                   */
  DBL   *d1emct;             /*[melem]  total elemental mass                          */
  INT   *i1elbe;             /*[melem] boundary element                               */
  INT   **i2elcn;            /* Z.LEI (06/05/2010 11:39:00)                           */
                             /* [2][melem] connective elements of each joint element  */  
};

typedef struct YDI_struct *YDI;
struct YDI_struct
{ INT micoup, nicoup;        /* maximum possible number of contacting couples         */
  INT    iiecff;             /* interaction element contact. couple free first        */

  DBL    diedi;              /* travel since last detection                           */
  DBL    diezon;             /* buffer zone size                                      */
  DBL   *d1iesl;             /*[mcoup] contact sliding distance between couples       */
  INT   *i1iecn;             /*[mcoup] next contacting couple                         */
  INT   *i1iect;             /*[mcoup] couple target                                  */
  INT   mistate;             /* number of states for d2sldis                          */
  DBL    **d2sldis;          /*[mistate][mcoup] sliding distance btw  couples         */
};

typedef struct YDN_struct *YDN;
struct YDN_struct
{ INT mnodim, nnodim;        /* max(actual) nodal dimensions number                   */
  INT mnopo, nnopo;          /* maximum (actual) number of nodal points               */

  DBL   *d1nmct;             /* [mnopo] nodal mass current translation                */
  DBL  **d2ncc;              /* [mnodim][mnopo] nodal coordinate current              */
  DBL  **d2nci;              /* [mnodim][mnopo] nodal coordinate initial              */
  DBL  **d2nfc;              /* [mnodim][mnopo] nodal contact force current           */
  DBL  **d2nft;              /* [mnodim][mnopo] nodal total current force             */  
  DBL  **d2nvc;              /* [mnodim][mnopo] nodal velocity current                */
  INT  *i1nobf;              /* [mnopo] nodal boundary >0 is boundary                 */
  INT  *i1nopr;              /* [mnopo] nodal boundary condition                      */
  INT  *i1niid;              /* [mnopo] nodal initial ID - Z.LEI 06/05/2010 11:16:26  */
  INT  *i1ncnp;              /* [mnopo] nodal initial connective relationship - Z.LEI */
  
};

typedef struct YDB_struct *YDB;
struct YDB_struct
{ INT mborh, nborh;          /* maximum (actual) number of boreholes                  */
  INT mbdim, nbdim;          /* maximum (actual) number of dimensions of boreholes    */
  INT nbpaf;                 /* number of amplitude factors                           */

  DBL   **d2bca;             /* [mbdim][mborh] coordinates of point A in borehole     */
  DBL   **d2bcb;             /* [mbdim][mborh] coordinates of point B in borehole     */
  DBL   *d1brad;             /* [mborh] radii of boreholes                            */
  DBL   *d1bpaf;             /* [nbpf] amplitude factor of pressure amplitude         */
  DBL   *d1bpts;             /* [mborh] start time of pressure load on boreholes      */
  DBL   *d1bpte;             /* [mborh] end time of pressure load on boreholes        */
  DBL   *d1bvdt;             /* [mborh] velocity of detonation on boreholes           */
  DBL   *d1bprs;             /* [mborh] amplitudes of pressure for each borehole      */
  DBL   dblmax;              /* max length of all boreholes                           */
  DBL   dbbuf;               /* buffer (max. size of element)                         */
};

typedef struct YDS_struct *YDS;
struct YDS_struct
{ INT msour, nsour;          /* maximum (actual) number of sources                    */
  INT msdim, nsdim;          /* maximum (actual) number of dimensions of sources      */
  INT nspaf;                 /* number of pressure amplitude factors (p=p(t) )        */
  INT nssaf;                 /* number of pressure amplitude factors (s=s(r) )        */

  DBL   **d2scs;             /* [msdim][msour] coordinates of sources                 */
  DBL   *d1spaf;             /* [nspaf] pressure amplitude factors (p=p(t) )          */
  DBL   *d1ssaf;             /* [nssaf] pressure amplitude factors (s=s(r) )          */
  DBL   *d1spts;             /* [msour] start time of pressure load                   */
  DBL   *d1spte;             /* [msour] end time of pressure load                     */
  DBL   *d1svpr;             /* [msour] velocity of pressure propagation              */
  DBL   *d1sprs;             /* [msour] amplitudes of pressures (p=p(t) )             */
  DBL   *d1ssir;             /* [msour] source initial radius                         */
  DBL   dsbuf;               /* buffer (max. size of element)                         */
};

typedef struct YDO_struct *YDO;
struct YDO_struct
{ INT mohys, nohys;          /* maximum (actual) number of hystory variables          */

  DBL     dohyp;             /* output hystory accuracy                               */

  DBL   *d1ohyf;             /*[mohys] output hystory factor to scale state           */
  DBL   *d1ohyc;             /*[mohys] output hystory factor to scale time            */
  DBL   *d1ohys;             /*[mohys] output hystory state                           */
  DBL   *d1ohyt;             /*[mohys] output hystory time                            */
  DBL   *d1ohyx;             /*[mohys] output history x coordinate of the point       */
  DBL   *d1ohyy;             /*[mohys] output history y coordinate of the point       */
  DBL   *d1ohyz;             /*[mohys] output history z coordinate of the point       */
  FILE  **f2ohyf;            /*[mohys] output history files                           */

  INT   *i1ohyt;             /*[mohys] output hystory type, i.e. which variable       */
};

/* Y Database Properties - Elements */
typedef struct YDPE_struct *YDPE;
struct YDPE_struct
{ INT mprop, nprop;          /* maximum (actual) number of properties                 */

  DBL   *d1peca;             /*[mprop] child age - procreation                        */
  DBL   *d1pecl;             /*[mprop] child life - interval for procreation          */
  DBL   *d1peks;             /*[mprop] dpeks=2hbeta*sqrt(E*ro) in 2D or 3D,0<beta<1   */
  DBL   *d1pela;             /*[mprop] property lamda - Lame elastic constant         */
  DBL   *d1pemu;             /*[mprop] property mu    - Lame elastic constant         */
  DBL   *d1pepe;             /*[mprop] contact penalty parameter                      */
  DBL   *d1pept;             /*[mprop] tangential penalty param., (0.1*ydpj->d1pjpe)  */
  DBL   *d1pefr;             /*[mprop] Coloumb friction                               */
  DBL   *d1pera;             /*[mprop] property radius of sphere                      */
  DBL   *d1pero;             /*[mprop] property ro    - density                       */
  DBL   *d1pevi;             /*[mprop] viscosity for  granular flow                   */

  DBL   *d1psem;             /*[mprop] maximum tensile stretch                        */

  INT   *i1pecn;             /*[mprop] property No to be assigned to child's nodes    */
  INT   *i1pecp;             /*[mprop] permanent property to be assigned to child     */
  INT   *i1pect;             /*[mprop] temporary property to be assigned to child     */
  INT   *i1pefr;             /*[mprop] if >, fracture                                 */
  INT   *i1pejp;             /*[mprop] joint property; if<0, no joints                */
  INT   *i1pemb;             /*[mprop] mark boundary nodes 1 yes 0 no                 */
  INT   *i1pemn;             /*[mprop] number of mesh refinements                     */

  INT   *i1pnib;             /*[mprop] 1 if borhole regardles boundary, else 0        */
  INT   *i1ptyp;             /*[mprop] property type                                  */ 
  INT   *i1psde;             /*[mprop] state damage elastic id                        */
  
  DBL    dkeneg;             /* kinetic energy                                        */
};

/* Y Database Properties for Joints */
typedef struct YDPJ_struct *YDPJ;
struct YDPJ_struct
{ INT mpjset, npjset;        /* maximum (actual) number of Joint Property sets       */

  INT  nintpt;               /* nb. integration points for joint elements            */        
  DBL *d1pjfs;               /* [mpjset] ultimate shear strength at joint            */
  DBL *d1pjft;               /* [mpjset] ultimate tensile strength at joint          */
  DBL *d1pjgf;               /* [mpjset] ultimate fracture energy at joint           */
  DBL *d1pjco;               /* [mpjset] cohesion at joint                           */
  DBL *d1pjfr;               /* [mpjset] friction at joint                           */

  DBL *d1pjpe;               /* [mpjset] fracture penalty parameter                  */

  INT *i1psde;               /* [mpjset] state damage elastic id                     */
  INT *i1ptyp;               /* [mpjset] Property Joint Type                         */
};

/* Y Database Properties - Nodes, Nodal Properties or Boundary Conditions for nodes  */
typedef struct YDPN_struct *YDPN;
struct YDPN_struct
{ INT mpnset, npnset;        /* maximum (actual) number of node property sets        */
  INT mpnfact, npnfact;      /* maximum (actual) number of factors                   */
  DBL ***d3pnfac;            /*[2][mbc][mfact] time and amplitude factor             */

  INT   *i1pnfx;             /*[mpnset] fixity x direction 1 force; 2 acc. 3 vel.    */
  INT   *i1pnfy;             /*[mpnset] fixity y direction 1 force; 2 acc. 3 vel.    */
  INT   *i1pnfz;             /*[mpnset] fixity z direction 1 force; 2 acc. 3 vel.    */

  DBL   *d1pnaf;             /*[mpnset] amplitude factor all ampltd multp by it      */
  DBL   *d1pnap;             /*[mpnset] amplitude of element surface pressure        */
  DBL   *d1pnat;             /*[mpnset] amplitude of element surface traction        */
  DBL   *d1pnax;             /*[mpnset] amplitude of force/velocity x                */
  DBL   *d1pnay;             /*[mpnset] amplitude of force/velocity y                */
  DBL   *d1pnaz;             /*[mpnset] amplitude of force/velocity z                */

  DBL   *d1pnxx;             /*[mpnset] direction of local x                         */
  DBL   *d1pnxy;             /*[mpnset] direction of local x                         */
  DBL   *d1pnxz;             /*[mpnset] direction of local x                         */
  DBL   *d1pnyx;             /*[mpnset] direction of local y                         */
  DBL   *d1pnyy;             /*[mpnset] direction of local y                         */
  DBL   *d1pnyz;             /*[mpnset] direction of local y                         */
  DBL   *d1pnzx;             /*[mpnset] direction of local z                         */
  DBL   *d1pnzy;             /*[mpnset] direction of local z                         */
  DBL   *d1pnzz;             /*[mpnset] direction of local z                         */
};

/* Y Database Properties - Meshing */
typedef struct YDPM_struct *YDPM;
struct YDPM_struct
{ INT mpmcom;                /* combination of pr. sets (rows in I2PMSET)            */
  INT mpmcol;                /* Number of columns  in I2PMSET                        */
  INT **i2pmset;             /*[mpmcol][mpmcom] mesh pr.sets i & j together          */
  INT mpmrow;                /* Number of rows in I2PMIJ                             */
  INT **i2pmij;              /* Meshing combinations defining joints btw 2 pr. sets  */
  INT iautbn;                /* identify boundary nodes automatically (zhou LEI) 
                                default > 0: identify boundary nodes automatically   */
};



typedef struct YD_struct *YD;
struct YD_struct
{ struct YDC_struct ydc;     /* control structure                                    */
  struct YDE_struct yde;     /* element description structure                        */
  struct YDI_struct ydi;     /* interaction structure                                */
  struct YDN_struct ydn;     /* node description structure                           */
  struct YDB_struct ydb;     /* borehole description structure                       */
  struct YDS_struct yds;     /* inter. fluid (source) description structure          */
  struct YDO_struct ydo;     /* output description structure                         */
  struct YDPE_struct ydpe;   /* property description structure for elements          */
  struct YDPN_struct ydpn;   /* property description structure for nodes             */
  struct YDPJ_struct ydpj;   /* property description structure for joints            */
  struct YDPM_struct ydpm;   /* property description structure for meshing           */
};

