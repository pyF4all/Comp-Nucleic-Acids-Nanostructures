#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "nabcode.h"
extern char NAB_rsbuf[];
static int mytaskid, numtasks;



static MOLECULE_T *m,  *m1,  *m_final,  *m_temp;

static INT_T i, j, k, t;

static REAL_T trise, ttwist;

static MATRIX_T mat_x, mat_drz;

static STRING_T *n = NULL,  *s[100];

static REAL_T re_arr, disp, re_arr2, re_arr3, repla;


int main( argc, argv )
	int	argc;
	char	*argv[];
{
	nabout = stdout; /*default*/

	mytaskid=0; numtasks=1;
static REAL_T __ft0001__;
static REAL_T __ft0002__;
static REAL_T __ft0003__;
static REAL_T __ft0004__;
static STRING_T *__st0001__ = NULL;
static STRING_T *__st0002__ = NULL;
static STRING_T *__st0003__ = NULL;
static STRING_T *__st0004__ = NULL;
static STRING_T *__st0005__ = NULL;
re_arr = 0.000000E+00;
re_arr2 = 4.000000E+00;
re_arr3 = 1.000000E+00;
repla =  - 2.000000E+00;
disp = 2.200000E+01;


m_final = newmolecule(  );
addstrand( m_final, "S" );
addstrand( m_final, "B" );
addstrand( m_final, "C" );
addstrand( m_final, "D" );
addstrand( m_final, "E" );
addstrand( m_final, "F" );
addstrand( m_final, "G" );
addstrand( m_final, "H" );
addstrand( m_final, "I" );
addstrand( m_final, "J" );
addstrand( m_final, "K" );
addstrand( m_final, "L" );
addstrand( m_final, "M" );
addstrand( m_final, "N" );
addstrand( m_final, "O" );
addstrand( m_final, "P" );
addstrand( m_final, "Q" );
addstrand( m_final, "R" );
addstrand( m_final, "My" );
addstrand( m_final, "Ny" );
addstrand( m_final, "Oy" );
addstrand( m_final, "Py" );
addstrand( m_final, "Qy" );
addstrand( m_final, "Ry" );



m1 = newmolecule(  );
addstrand( m1, "A" );


for( i = 5;i <= 54;i = i + 1 ){
trise = ( i ) * 2.810000E+00;
ttwist = ( i ) * 3.270000E+01;

NAB_strcpy(  &s[5 - 1], "u" );NAB_strcpy(  &s[6 - 1], "u" );NAB_strcpy(  &s[7 - 1], "c" );NAB_strcpy(  &s[8 - 1], "u" );NAB_strcpy(  &s[9 - 1], "a" );NAB_strcpy(  &s[10 - 1], "u" );NAB_strcpy(  &s[11 - 1], "c" );NAB_strcpy(  &s[12 - 1], "u" );NAB_strcpy(  &s[13 - 1], "u" );NAB_strcpy(  &s[14 - 1], "a" );
NAB_strcpy(  &s[15 - 1], "c" );NAB_strcpy(  &s[16 - 1], "a" );NAB_strcpy(  &s[17 - 1], "u" );NAB_strcpy(  &s[18 - 1], "u" );NAB_strcpy(  &s[19 - 1], "u" );NAB_strcpy(  &s[20 - 1], "c" );NAB_strcpy(  &s[21 - 1], "u" );NAB_strcpy(  &s[22 - 1], "u" );NAB_strcpy(  &s[23 - 1], "c" );NAB_strcpy(  &s[24 - 1], "g" );
NAB_strcpy(  &s[25 - 1], "c" );NAB_strcpy(  &s[26 - 1], "c" );NAB_strcpy(  &s[27 - 1], "g" );NAB_strcpy(  &s[28 - 1], "u" );NAB_strcpy(  &s[29 - 1], "c" );NAB_strcpy(  &s[30 - 1], "u" );NAB_strcpy(  &s[31 - 1], "u" );NAB_strcpy(  &s[32 - 1], "u" );NAB_strcpy(  &s[33 - 1], "a" );NAB_strcpy(  &s[34 - 1], "a" );
NAB_strcpy(  &s[35 - 1], "u" );NAB_strcpy(  &s[36 - 1], "u" );NAB_strcpy(  &s[37 - 1], "g" );NAB_strcpy(  &s[38 - 1], "u" );NAB_strcpy(  &s[39 - 1], "c" );NAB_strcpy(  &s[40 - 1], "a" );NAB_strcpy(  &s[41 - 1], "a" );NAB_strcpy(  &s[42 - 1], "c" );NAB_strcpy(  &s[43 - 1], "u" );NAB_strcpy(  &s[44 - 1], "u" );
NAB_strcpy(  &s[45 - 1], "c" );NAB_strcpy(  &s[46 - 1], "a" );NAB_strcpy(  &s[47 - 1], "u" );NAB_strcpy(  &s[48 - 1], "a" );NAB_strcpy(  &s[49 - 1], "u" );NAB_strcpy(  &s[50 - 1], "c" );NAB_strcpy(  &s[51 - 1], "c" );NAB_strcpy(  &s[52 - 1], "u" );NAB_strcpy(  &s[53 - 1], "g" );NAB_strcpy(  &s[54 - 1], "g" );
n = NULL;

NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );

printf( "%s %s\n", s[i - 1], n );
if( i == 5 ){
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "rna" ),  &s[i - 1], STEMP( __st0003__, "" ), STEMP( __st0004__, "rna" ), FTEMP( __ft0001__, 6.860000E+00 ), FTEMP( __ft0002__, 1.671000E+01 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "a3" ) );
}
else if( i == 54 ){
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "rna" ),  &s[i - 1], STEMP( __st0003__, "" ), STEMP( __st0004__, "rna" ), FTEMP( __ft0001__, 6.860000E+00 ), FTEMP( __ft0002__, 1.671000E+01 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "a5" ) );
}
else{
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "rna" ),  &s[i - 1], STEMP( __st0003__, "" ), STEMP( __st0004__, "rna" ), FTEMP( __ft0001__, 6.860000E+00 ), FTEMP( __ft0002__, 1.671000E+01 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}

transformmol( mat_drz, m, NULL );
mergestr( m1, "A", "first", m, "anti", "last" );
freemolecule( m );
n = NULL;
}

m_temp = copymolecule( m1 );
NAB_matcpy( mat_x, newtransform( disp, 0.000000E+00, repla, 0.000000E+00, 0.000000E+00, re_arr3 * 3.270000E+01 ) );
transformmol( mat_x, m_temp, NULL );
putpdb( "rna1.pdb", m_temp, NULL );
mergestr( m_final, "S", "last", m_temp, "A", "first" );


freemolecule( m_temp );
freemolecule( m1 );
m1 = newmolecule(  );
addstrand( m1, "A" );
n = NULL;


for( i = 1;i <= 50;i = i + 1 ){
trise = ( i - 1 ) * 2.810000E+00;
ttwist = ( i - 1 ) * 3.270000E+01;
NAB_strcpy(  &s[1 - 1], "g" );NAB_strcpy(  &s[2 - 1], "g" );NAB_strcpy(  &s[3 - 1], "u" );NAB_strcpy(  &s[4 - 1], "c" );NAB_strcpy(  &s[5 - 1], "c" );NAB_strcpy(  &s[6 - 1], "u" );NAB_strcpy(  &s[7 - 1], "a" );NAB_strcpy(  &s[8 - 1], "u" );NAB_strcpy(  &s[9 - 1], "a" );NAB_strcpy(  &s[10 - 1], "c" );
NAB_strcpy(  &s[11 - 1], "u" );NAB_strcpy(  &s[12 - 1], "u" );NAB_strcpy(  &s[13 - 1], "c" );NAB_strcpy(  &s[14 - 1], "a" );NAB_strcpy(  &s[15 - 1], "a" );NAB_strcpy(  &s[16 - 1], "c" );NAB_strcpy(  &s[17 - 1], "u" );NAB_strcpy(  &s[18 - 1], "g" );NAB_strcpy(  &s[19 - 1], "u" );NAB_strcpy(  &s[20 - 1], "u" );
NAB_strcpy(  &s[21 - 1], "a" );NAB_strcpy(  &s[22 - 1], "a" );NAB_strcpy(  &s[23 - 1], "u" );NAB_strcpy(  &s[24 - 1], "u" );NAB_strcpy(  &s[25 - 1], "u" );NAB_strcpy(  &s[26 - 1], "c" );NAB_strcpy(  &s[27 - 1], "u" );NAB_strcpy(  &s[28 - 1], "g" );NAB_strcpy(  &s[29 - 1], "c" );NAB_strcpy(  &s[30 - 1], "c" );
NAB_strcpy(  &s[31 - 1], "g" );NAB_strcpy(  &s[32 - 1], "c" );NAB_strcpy(  &s[33 - 1], "u" );NAB_strcpy(  &s[34 - 1], "u" );NAB_strcpy(  &s[35 - 1], "c" );NAB_strcpy(  &s[36 - 1], "u" );NAB_strcpy(  &s[37 - 1], "u" );NAB_strcpy(  &s[38 - 1], "u" );NAB_strcpy(  &s[39 - 1], "a" );NAB_strcpy(  &s[40 - 1], "c" );
NAB_strcpy(  &s[41 - 1], "a" );NAB_strcpy(  &s[42 - 1], "u" );NAB_strcpy(  &s[43 - 1], "u" );NAB_strcpy(  &s[44 - 1], "c" );NAB_strcpy(  &s[45 - 1], "u" );NAB_strcpy(  &s[46 - 1], "a" );NAB_strcpy(  &s[47 - 1], "u" );NAB_strcpy(  &s[48 - 1], "c" );NAB_strcpy(  &s[49 - 1], "u" );NAB_strcpy(  &s[50 - 1], "u" );
n = NULL;
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );


printf( "%s %s\n", s[i - 1], n );
if( i == 1 ){
m = wc_helix(  &s[i - 1], STEMP( __st0001__, "" ), STEMP( __st0002__, "rna" ),  &n, STEMP( __st0003__, "" ), STEMP( __st0004__, "rna" ), FTEMP( __ft0001__, 6.860000E+00 ), FTEMP( __ft0002__, 1.671000E+01 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "s5" ) );
}
else if( i == 50 ){
m = wc_helix(  &s[i - 1], STEMP( __st0001__, "" ), STEMP( __st0002__, "rna" ),  &n, STEMP( __st0003__, "" ), STEMP( __st0004__, "rna" ), FTEMP( __ft0001__, 6.860000E+00 ), FTEMP( __ft0002__, 1.671000E+01 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "s3" ) );
}
else{
m = wc_helix(  &s[i - 1], STEMP( __st0001__, "" ), STEMP( __st0002__, "rna" ),  &n, STEMP( __st0003__, "" ), STEMP( __st0004__, "rna" ), FTEMP( __ft0001__, 6.860000E+00 ), FTEMP( __ft0002__, 1.671000E+01 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}

transformmol( mat_drz, m, NULL );
mergestr( m1, "A", "last", m, "sense", "first" );
freemolecule( m );
n = NULL;
}
m_temp = copymolecule( m1 );
NAB_matcpy( mat_x, newtransform( 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, re_arr2 * 3.270000E+01 ) );
transformmol( mat_x, m_temp, NULL );
putpdb( "rna2.pdb", m_temp, NULL );
mergestr( m_final, "B", "last", m_temp, "A", "first" );




freemolecule( m1 );
m1 = newmolecule(  );
addstrand( m1, "A" );
for( i = 1;i <= 28;i = i + 1 ){

n = NULL;
NAB_strcpy(  &s[1 - 1], "g" );NAB_strcpy(  &s[2 - 1], "g" );NAB_strcpy(  &s[3 - 1], "a" );NAB_strcpy(  &s[4 - 1], "c" );NAB_strcpy(  &s[5 - 1], "c" );NAB_strcpy(  &s[6 - 1], "a" );NAB_strcpy(  &s[7 - 1], "a" );NAB_strcpy(  &s[8 - 1], "g" );NAB_strcpy(  &s[9 - 1], "a" );NAB_strcpy(  &s[10 - 1], "u" );
NAB_strcpy(  &s[11 - 1], "a" );NAB_strcpy(  &s[12 - 1], "g" );NAB_strcpy(  &s[13 - 1], "a" );NAB_strcpy(  &s[14 - 1], "a" );NAB_strcpy(  &s[15 - 1], "u" );NAB_strcpy(  &s[16 - 1], "g" );NAB_strcpy(  &s[17 - 1], "a" );NAB_strcpy(  &s[18 - 1], "g" );NAB_strcpy(  &s[19 - 1], "u" );NAB_strcpy(  &s[20 - 1], "u" );
NAB_strcpy(  &s[21 - 1], "g" );NAB_strcpy(  &s[22 - 1], "a" );NAB_strcpy(  &s[23 - 1], "a" );NAB_strcpy(  &s[24 - 1], "g" );NAB_strcpy(  &s[25 - 1], "u" );NAB_strcpy(  &s[26 - 1], "a" );NAB_strcpy(  &s[27 - 1], "u" );NAB_strcpy(  &s[28 - 1], "a" );
if( i <= 16 ){
trise = ( 55 - i ) * 2.810000E+00;
t = ( 55 - i );
ttwist = t * 3.270000E+01;

NAB_matcpy( mat_x, newtransform( 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, re_arr2 * 3.270000E+01 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
if( i == 1 ){
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "rna" ),  &s[i - 1], STEMP( __st0003__, "" ), STEMP( __st0004__, "rna" ), FTEMP( __ft0001__, 6.860000E+00 ), FTEMP( __ft0002__, 1.671000E+01 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "a5" ) );
}
else{
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "rna" ),  &s[i - 1], STEMP( __st0003__, "" ), STEMP( __st0004__, "rna" ), FTEMP( __ft0001__, 6.860000E+00 ), FTEMP( __ft0002__, 1.671000E+01 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}
transformmol( mat_drz, m, NULL );
transformmol( mat_x, m, NULL );
mergestr( m1, "A", "last", m, "anti", "first" );
}
else{
trise = ( 21 + i ) * 2.810000E+00;
t = ( 21 + i );
ttwist = t * 3.270000E+01;
NAB_matcpy( mat_x, newtransform( disp, 0.000000E+00, repla, 0.000000E+00, 0.000000E+00, re_arr3 * 3.270000E+01 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
if( i == 28 ){
m = wc_helix(  &s[i - 1], STEMP( __st0001__, "" ), STEMP( __st0002__, "rna" ),  &n, STEMP( __st0003__, "" ), STEMP( __st0004__, "rna" ), FTEMP( __ft0001__, 6.860000E+00 ), FTEMP( __ft0002__, 1.671000E+01 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "s3" ) );
}
else{
m = wc_helix(  &s[i - 1], STEMP( __st0001__, "" ), STEMP( __st0002__, "rna" ),  &n, STEMP( __st0003__, "" ), STEMP( __st0004__, "rna" ), FTEMP( __ft0001__, 6.860000E+00 ), FTEMP( __ft0002__, 1.671000E+01 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}
transformmol( mat_drz, m, NULL );
transformmol( mat_x, m, NULL );
mergestr( m1, "A", "last", m, "sense", "first" );

}
freemolecule( m );
n = NULL;
printf( "%s %s %d %d %lf\n", s[i - 1], n, i, t, ttwist );
}
putpdb( "ssd1.pdb", m1, NULL );
mergestr( m_final, "C", "last", m1, "A", "first" );
putpdb( "6h1.pdb", m_final, NULL );

freemolecule( m1 );
m1 = newmolecule(  );
addstrand( m1, "A" );
for( i = 1;i <= 44;i = i + 1 ){

n = NULL;
NAB_strcpy(  &s[1 - 1], "c" );NAB_strcpy(  &s[2 - 1], "a" );NAB_strcpy(  &s[3 - 1], "g" );NAB_strcpy(  &s[4 - 1], "a" );NAB_strcpy(  &s[5 - 1], "a" );NAB_strcpy(  &s[6 - 1], "a" );NAB_strcpy(  &s[7 - 1], "u" );NAB_strcpy(  &s[8 - 1], "u" );NAB_strcpy(  &s[9 - 1], "a" );NAB_strcpy(  &s[10 - 1], "a" );
NAB_strcpy(  &s[11 - 1], "c" );NAB_strcpy(  &s[12 - 1], "u" );NAB_strcpy(  &s[13 - 1], "a" );NAB_strcpy(  &s[14 - 1], "a" );NAB_strcpy(  &s[15 - 1], "a" );NAB_strcpy(  &s[16 - 1], "g" );NAB_strcpy(  &s[17 - 1], "a" );NAB_strcpy(  &s[18 - 1], "a" );NAB_strcpy(  &s[19 - 1], "g" );NAB_strcpy(  &s[20 - 1], "c" );
NAB_strcpy(  &s[21 - 1], "g" );NAB_strcpy(  &s[22 - 1], "g" );NAB_strcpy(  &s[23 - 1], "c" );NAB_strcpy(  &s[24 - 1], "a" );NAB_strcpy(  &s[25 - 1], "g" );NAB_strcpy(  &s[26 - 1], "a" );NAB_strcpy(  &s[27 - 1], "a" );NAB_strcpy(  &s[28 - 1], "a" );NAB_strcpy(  &s[29 - 1], "u" );NAB_strcpy(  &s[30 - 1], "u" );
NAB_strcpy(  &s[31 - 1], "a" );NAB_strcpy(  &s[32 - 1], "a" );NAB_strcpy(  &s[33 - 1], "c" );NAB_strcpy(  &s[34 - 1], "u" );NAB_strcpy(  &s[35 - 1], "a" );NAB_strcpy(  &s[36 - 1], "a" );NAB_strcpy(  &s[37 - 1], "a" );NAB_strcpy(  &s[38 - 1], "g" );NAB_strcpy(  &s[39 - 1], "a" );NAB_strcpy(  &s[40 - 1], "a" );
NAB_strcpy(  &s[41 - 1], "g" );NAB_strcpy(  &s[42 - 1], "c" );NAB_strcpy(  &s[43 - 1], "g" );NAB_strcpy(  &s[44 - 1], "g" );
if( i <= 11 ){
trise = ( 28 - i ) * 2.810000E+00;
t = ( 28 - i );
ttwist = t * 3.270000E+01;

NAB_matcpy( mat_x, newtransform( 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, re_arr2 * 3.270000E+01 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
if( i == 1 ){
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "rna" ),  &s[i - 1], STEMP( __st0003__, "" ), STEMP( __st0004__, "rna" ), FTEMP( __ft0001__, 6.860000E+00 ), FTEMP( __ft0002__, 1.671000E+01 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "a5" ) );
}
else{
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "rna" ),  &s[i - 1], STEMP( __st0003__, "" ), STEMP( __st0004__, "rna" ), FTEMP( __ft0001__, 6.860000E+00 ), FTEMP( __ft0002__, 1.671000E+01 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}
transformmol( mat_drz, m, NULL );
transformmol( mat_x, m, NULL );
mergestr( m1, "A", "last", m, "anti", "first" );
}
else if( i > 11 && i <= 33 ){
trise = ( 4 + i ) * 2.810000E+00;
t = ( 4 + i );
ttwist = t * 3.270000E+01;
NAB_matcpy( mat_x, newtransform( disp, 0.000000E+00, repla, 0.000000E+00, 0.000000E+00, re_arr3 * 3.270000E+01 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
if( i == 20 ){
m = wc_helix(  &s[i - 1], STEMP( __st0001__, "" ), STEMP( __st0002__, "rna" ),  &n, STEMP( __st0003__, "" ), STEMP( __st0004__, "rna" ), FTEMP( __ft0001__, 6.860000E+00 ), FTEMP( __ft0002__, 1.671000E+01 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}
else{
m = wc_helix(  &s[i - 1], STEMP( __st0001__, "" ), STEMP( __st0002__, "rna" ),  &n, STEMP( __st0003__, "" ), STEMP( __st0004__, "rna" ), FTEMP( __ft0001__, 6.860000E+00 ), FTEMP( __ft0002__, 1.671000E+01 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}
transformmol( mat_drz, m, NULL );
transformmol( mat_x, m, NULL );
mergestr( m1, "A", "last", m, "sense", "first" );

}
else{
trise = ( 72 - i ) * 2.810000E+00;
t = ( 72 - i );
ttwist = t * 3.270000E+01;

NAB_matcpy( mat_x, newtransform( 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, re_arr2 * 3.270000E+01 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
if( i == 44 ){
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "rna" ),  &s[i - 1], STEMP( __st0003__, "" ), STEMP( __st0004__, "rna" ), FTEMP( __ft0001__, 6.860000E+00 ), FTEMP( __ft0002__, 1.671000E+01 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "a3" ) );
}
else{
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "rna" ),  &s[i - 1], STEMP( __st0003__, "" ), STEMP( __st0004__, "rna" ), FTEMP( __ft0001__, 6.860000E+00 ), FTEMP( __ft0002__, 1.671000E+01 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}
transformmol( mat_drz, m, NULL );
transformmol( mat_x, m, NULL );
mergestr( m1, "A", "last", m, "anti", "first" );
}
freemolecule( m );
n = NULL;
printf( "%s %s %d %d %lf\n", s[i - 1], n, i, t, ttwist );
}
putpdb( "ssd2.pdb", m1, NULL );
mergestr( m_final, "D", "last", m1, "A", "first" );
putpdb( "6h2.pdb", m_final, NULL );





freemolecule( m1 );
m1 = newmolecule(  );
addstrand( m1, "A" );
for( i = 1;i <= 28;i = i + 1 ){

n = NULL;
NAB_strcpy(  &s[1 - 1], "g" );NAB_strcpy(  &s[2 - 1], "g" );NAB_strcpy(  &s[3 - 1], "a" );NAB_strcpy(  &s[4 - 1], "c" );NAB_strcpy(  &s[5 - 1], "c" );NAB_strcpy(  &s[6 - 1], "a" );NAB_strcpy(  &s[7 - 1], "a" );NAB_strcpy(  &s[8 - 1], "g" );NAB_strcpy(  &s[9 - 1], "a" );NAB_strcpy(  &s[10 - 1], "u" );
NAB_strcpy(  &s[11 - 1], "a" );NAB_strcpy(  &s[12 - 1], "g" );NAB_strcpy(  &s[13 - 1], "a" );NAB_strcpy(  &s[14 - 1], "a" );NAB_strcpy(  &s[15 - 1], "u" );NAB_strcpy(  &s[16 - 1], "g" );NAB_strcpy(  &s[17 - 1], "a" );NAB_strcpy(  &s[18 - 1], "g" );NAB_strcpy(  &s[19 - 1], "u" );NAB_strcpy(  &s[20 - 1], "u" );
NAB_strcpy(  &s[21 - 1], "g" );NAB_strcpy(  &s[22 - 1], "a" );NAB_strcpy(  &s[23 - 1], "a" );NAB_strcpy(  &s[24 - 1], "g" );NAB_strcpy(  &s[25 - 1], "u" );NAB_strcpy(  &s[26 - 1], "a" );NAB_strcpy(  &s[27 - 1], "u" );NAB_strcpy(  &s[28 - 1], "a" );
if( i <= 16 ){
trise = ( i - 1 ) * 2.810000E+00;
t = ( i - 1 );
ttwist = t * 3.270000E+01;
NAB_matcpy( mat_x, newtransform( disp, 0.000000E+00, repla, 0.000000E+00, 0.000000E+00, re_arr3 * 3.270000E+01 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
if( i == 1 ){
m = wc_helix(  &s[i - 1], STEMP( __st0001__, "" ), STEMP( __st0002__, "rna" ),  &n, STEMP( __st0003__, "" ), STEMP( __st0004__, "rna" ), FTEMP( __ft0001__, 6.860000E+00 ), FTEMP( __ft0002__, 1.671000E+01 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "s5" ) );
}
else{
m = wc_helix(  &s[i - 1], STEMP( __st0001__, "" ), STEMP( __st0002__, "rna" ),  &n, STEMP( __st0003__, "" ), STEMP( __st0004__, "rna" ), FTEMP( __ft0001__, 6.860000E+00 ), FTEMP( __ft0002__, 1.671000E+01 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}
transformmol( mat_drz, m, NULL );
transformmol( mat_x, m, NULL );
mergestr( m1, "A", "last", m, "sense", "first" );
}
else{
trise = ( 33 - i ) * 2.810000E+00;
t = ( 33 - i );
ttwist = t * 3.270000E+01;
NAB_matcpy( mat_x, newtransform( 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, re_arr2 * 3.270000E+01 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
if( i == 28 ){
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "rna" ),  &s[i - 1], STEMP( __st0003__, "" ), STEMP( __st0004__, "rna" ), FTEMP( __ft0001__, 6.860000E+00 ), FTEMP( __ft0002__, 1.671000E+01 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "a3" ) );
}
else{
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "rna" ),  &s[i - 1], STEMP( __st0003__, "" ), STEMP( __st0004__, "rna" ), FTEMP( __ft0001__, 6.860000E+00 ), FTEMP( __ft0002__, 1.671000E+01 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}
transformmol( mat_drz, m, NULL );
transformmol( mat_x, m, NULL );
mergestr( m1, "A", "last", m, "anti", "first" );

}
freemolecule( m );
n = NULL;
printf( "%s %s %d %d %lf\n", s[i - 1], n, i, t, ttwist );
}
putpdb( "ssd3.pdb", m1, NULL );
mergestr( m_final, "E", "last", m1, "A", "first" );
putpdb( "6h3.pdb", m_final, NULL );


	exit( 0 );
}
