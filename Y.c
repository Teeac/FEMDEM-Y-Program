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
/* File   Y.c */
#include "Yproto.h"
void show_time_info(YDC ydc);
INT argc;
main(argc,argv)
char **argv;
{ CHR c1name[300];         /* name of the problem i.e. input file */
  struct YD_struct yd;     /* Y database                          */
  YDC ydc=&(yd.ydc);       /* Y control database                  */
  YDE yde=&(yd.yde);       /* Y element database                  */
  YDI ydi=&(yd.ydi);       /* Y interaction database              */
  YDN ydn=&(yd.ydn);       /* Y node database                     */
  YDB ydb=&(yd.ydb);       /* Y borehole database                 */
  YDS yds=&(yd.yds);       /* Y source (inter. fluid) database    */
  YDO ydo=&(yd.ydo);       /* Y output database                   */
  YDPE ydpe=&(yd.ydpe);    /* Y property database  for elements   */
  YDPN ydpn=&(yd.ydpn);    /* Y property database  for nodes (BC) */
  YDPJ ydpj=&(yd.ydpj);    /* Y property database  for joints     */
  YDPM ydpm=&(yd.ydpm);    /* Y property database  for meshing    */
  INT Tctrlc, itimes=0;
  
  /* process data */
  CHRw(stdout,"Copyright (C) 2000, Dr. Antonio Munjiza. \n");
  CHRw(stdout,"This code is provided as part of the book entitled The Combined \n");
  CHRw(stdout,"Finite Discrete Element Method. It is distributed WITHOUT ANY WARRANTY; \n");
  CHRw(stdout,"without even the implied warranty of MERCHANTABILITY or FITNESS FOR A \n");
  CHRw(stdout,"PARTICULAR PURPOSE. Inclusion of a part or whole of this code into any other \n");
  CHRw(stdout,"commercial or research or other purpose code is not granted without author's \n");
  CHRw(stdout,"written explicit permission. \n");
  CHRw(stdout,"When results using whole or any part of this code \n");
  CHRw(stdout,"are published, Y code must be mentioned and acknowledgement to \n");
  CHRw(stdout,"Dr Munjiza must be made. \n");
  CHRw(stdout,"Should you modify this source code, the Copyright (C) on the modified code \n");
  CHRw(stdout,"as a whole belongs to Dr. A. Munjiza regardless of the extent or nature \n");
  CHRw(stdout,"of modifications. \n");
  CHRw(stdout,"Copyright (C) to whole of any code containing any part of this code \n");
  CHRw(stdout,"also belongs to Dr. A.Munjiza. \n");
  CHRw(stdout,"Any code comprising any part of this source code \n");
  CHRw(stdout,"must be called Y program.\n");
  CHRw(stdout,"If you do not agree with this, you are not allowed to do \n");
  CHRw(stdout,"any modifications to any part of this source code or include \n");
  CHRw(stdout,"any part of it in any other program. \n\n");
    
  /* get name of the problem */
  if(argv[1]!=NULL)
  { CHRcpy(c1name,argv[1]);
  }
  else
  { CHRwcr(stdout);
    CHRw(stdout,"  please define input file names: "); CHRwcr(stdout);
    CHRw(stdout," >");
    gets(c1name);
  }
  strcpy(ydc->cfiname, c1name);   ydc->cfiname[255]='\0';
  ydc->finp=FILENULL; ydc->fcheck=FILENULL;
  
  /* Process while any input */
  while(Yrd(c1name,&yd)>0)
  { itimes=itimes+1;
    CHRw(stdout,"NEW INPUT: "); CHRw(stdout, c1name); CHRwcr(stdout);
    if(Ycheck(&yd)<0) break;    date_and_time(ydc->cruntime); timestamp();
    CHRw(stdout, "Start calculating ...\n");
    for(ydc->ncstep=ydc->ncstep;ydc->ncstep<ydc->mcstep;ydc->ncstep++)
    { show_time_info(ydc);                            /* show time information    */
      Ymd(ydc,yde,ydi,ydn,ydpe, ydpn,ydpm);           /* mesh elements            */
      Yfd(ydc,yde,ydn,ydi,ydo,ydpe,ydpn,ydpj);        /* nodal forces             */
      Ybor(ydc,yde,ydn,ydb,yds,ydpe,ydpj,ydpn);       /* borholes, inter. fluid   */
      Ycd(ydc,yde,ydi,ydn,ydpe,ydpn);                 /* contact detection        */
      Yid(ydc,yde,ydi,ydn,ydo,ydpe,ydpn, ydpj,ydpm);  /* interaction              */
      Yod(c1name,&yd);                                /* output results           */
      Ysd(ydc,yde,ydn,ydo,ydpe,ydpn );                /* solve equations          */
      Yfrd(ydc,yde,ydi,ydn,ydpe,ydpn,ydpj,ydpm);      /* fracture                 */
      ydc->dctime=ydc->dctime+ydc->dcstec;            /* update time              */
      /* CTRL-C Interruption */
      Tctrlc = enablc(ydc->dctime, ydc->ncstep, ydc->mcstep);    if(Tctrlc!=1) break;
    }
  }
  
  /* Termination */
  CHRw(stderr,"   ***** Y HAS ORDERLY FINISHED *****");  CHRwcr(stderr);
  CHRw(stderr,"Press a key to continue");  CHRwcr(stderr);
  getchar();
}

void show_time_info(YDC ydc)
{ INT static itimes = 0;
  FILE *fptr=FILENULL, *fcheck=FILENULL;
  CHR  datetime[300];
  DBL runtime, ctime;
  INT ncstep, icshtf, il;
 
  fcheck = ydc->fcheck;
  ncstep = ydc->ncstep;
  ctime =  ydc->dctime;
  icshtf = ydc->icshtf;
  
  if(itimes==0)
  { itimes = 1;
    for(il=0; il<2; il++)
    { fptr=stdout; if(il==1) fptr=fcheck;
      CHRww(fptr, "Ncycle", 8); CHRww(fptr, "Physi. time", 17); 
      CHRww(fptr, "Date and Time", 22);
      CHRww(fptr, "Elapsed time", 15);
      CHRwcr(fptr);
    }
  }
 
  if(ncstep%icshtf==0)
  { for(il=0; il<2; il++)
    { fptr=stdout; if(il==1) fptr=fcheck;
      date_and_time(datetime);
      runtime = run_time();
      INTw(fptr, ydc->ncstep, 8); fprintf(fptr, "  %15.7e",ctime); CHRwsp(fptr);
      CHRwsp(fptr);  CHRwsp(fptr);  CHRw(fptr, datetime);  
      fprintf(fptr, "%15.3f", runtime);
      CHRwcr(fptr);
    }
  }
}
