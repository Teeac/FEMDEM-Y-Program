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
/* File  Yproto.h */
#ifndef YDINCL
#include "Yd.h"
#define YDINCL
#endif

void Ycd(         /*   contact detection                   */
#if NeedFunctionPrototypes
 YDC ydc,         /* = control database                    */
 YDE yde,         /* = element database                    */
 YDI ydi,         /* = interaction database                */
 YDN ydn,         /* = nodal database                      */
 YDPE ydpe,       /* - property database for elements      */
 YDPN ydpn        /* - boundary conditions database        */
 #endif
);

void Yfd(         /*   nodal forces                        */
#if NeedFunctionPrototypes
 YDC ydc,         /* - control database                    */
 YDE yde,         /* = element database                    */
 YDN ydn,         /* = nodal database                      */
 YDI ydi,         /* = interaction database                */
 YDO ydo,         /* - output database                     */
 YDPE ydpe,       /* - property database for elements      */
 YDPN ydpn,       /* - property database for nodes (BC)    */
 YDPJ ydpj        /* - property database for elements      */
#endif
);

void Ybor(        /* pressure on borholes                  */
#if NeedFunctionPrototypes
 YDC ydc,         /* = control database                    */
 YDE yde,         /* = element database                    */
 YDN ydn,         /* = nodal database                      */
 YDB ydb,         /* = borhole database                    */
 YDS yds,         /* = source (inter. fluid) database      */
 YDPE ydpe,       /* - property database for elements      */
 YDPJ ydpj,       /* - property database for joints        */
 YDPN ydpn        /* - boundary conditions database        */
#endif
);

void Yfrd(        /*  fracture elements                    */
#if NeedFunctionPrototypes
 YDC ydc,         /* = control database                    */
 YDE yde,         /* = element database                    */
 YDI ydi,         /* = interaction database                */
 YDN ydn,         /* = nodal database                      */
 YDPE ydpe,       /* - property database for elements      */
 YDPN ydpn,       /* - boundary conditions database        */
 YDPJ ydpj,       /* - property database for elements      */
 YDPM ydpm        /* - property database for meshing       */
#endif
);

void Ymd(         /*  mesh i.e. subdivide elements         */
#if NeedFunctionPrototypes
 YDC ydc,         /* = control database                    */
 YDE yde,         /* = element database                    */
 YDI ydi,         /* = interaction database                */
 YDN ydn,         /* = nodal database                      */
 YDPE ydpe,       /* - property database for elements      */
 YDPN ydpn,       /* - property database for nodes (BC)    */
 YDPM ydpm        /* - property database for meshing       */
#endif
);

void Yid(         /*   procreate                           */
#if NeedFunctionPrototypes
 YDC ydc,         /* = control database                    */
 YDE yde,         /* = element database                    */
 YDI ydi,         /* = interaction database                */
 YDN ydn,         /* = nodal database                      */
 YDO ydo,         /* - output database                     */
 YDPE ydpe,       /* - property database for elements      */
 YDPN ydpn,       /* - boundary conditions database        */
 YDPJ ydpj,       /* - property database for elements      */
 YDPM ydpm        /* - property database for meshing       */
#endif
);

void Ysd(         /*  solve equations                      */
#if NeedFunctionPrototypes
 YDC ydc,         /* - control database                    */
 YDE yde,         /* = element database                    */
 YDN ydn,         /* = nodal database                      */
 YDO ydo,         /* = output database                     */
 YDPE ydpe,       /* - property database for elements      */
 YDPN ydpn        /* - boundary conditions database        */
#endif
);

INT Yrd(          /*  read input                           */
#if NeedFunctionPrototypes
 CHR *namep,      /* - name of the problem i.e. input file */
 YD  yd           /* = database                            */
#endif
);

INT Ycheck(         /* check model data before calculation            */
#if NeedFunctionPrototypes
 YD  yd           /* = database                                       */
#endif
);

void Yod(         /* output results in space-save format   */
#if NeedFunctionPrototypes
 CHR *namep,      /* - name of the problem i.e. input file */
 YD  yd           /* = database                            */
#endif
);
