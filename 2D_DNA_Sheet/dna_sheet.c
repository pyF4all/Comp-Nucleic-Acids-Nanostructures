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

static STRING_T *n = NULL,  *s = NULL;

static REAL_T re_arr, disp;


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
disp = 2.000000E+01;


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


for( i = 1;i <= 56;i = i + 1 ){
trise = ( i ) * 3.380000E+00;
ttwist = ( i ) * 3.428000E+01;

NAB_strcpy(  &s, "g" );
n = NULL;

NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );

printf( "%s %s\n", s, n );
if( i == 1 ){
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &s, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "a3" ) );
}
else if( i == 56 ){
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &s, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "a5" ) );
}
else{
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &s, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}

transformmol( mat_drz, m, NULL );
mergestr( m1, "A", "first", m, "anti", "last" );
freemolecule( m );
n = NULL;
NAB_strcpy(  &s, "g" );
}

m_temp = copymolecule( m1 );
NAB_matcpy( mat_x, newtransform( disp, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, re_arr * 3.428000E+01 ) );
transformmol( mat_x, m_temp, NULL );

mergestr( m_final, "S", "last", m_temp, "A", "first" );

m_temp = copymolecule( m1 );
NAB_matcpy( mat_x, newtransform( 3 * disp, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, re_arr * 3.428000E+01 ) );
transformmol( mat_x, m_temp, NULL );

mergestr( m_final, "B", "last", m_temp, "A", "first" );

m_temp = copymolecule( m1 );
NAB_matcpy( mat_x, newtransform( 5 * disp, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, re_arr * 3.428000E+01 ) );
transformmol( mat_x, m_temp, NULL );
mergestr( m_final, "C", "last", m_temp, "A", "first" );


freemolecule( m_temp );
freemolecule( m1 );
m1 = newmolecule(  );
addstrand( m1, "A" );
n = NULL;
NAB_strcpy(  &s, "g" );

for( i = 1;i <= 56;i = i + 1 ){
trise = ( i ) * 3.380000E+00;
ttwist = ( i ) * 3.428000E+01;

NAB_strcpy(  &s, "g" );
n = NULL;
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );


printf( "%s %s\n", s, n );
if( i == 1 ){
m = wc_helix(  &s, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &n, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "s5" ) );
}
else if( i == 56 ){
m = wc_helix(  &s, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &n, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "s3" ) );
}
else{
m = wc_helix(  &s, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &n, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}

transformmol( mat_drz, m, NULL );
mergestr( m1, "A", "last", m, "sense", "first" );
freemolecule( m );
n = NULL;
NAB_strcpy(  &s, "g" );
}
m_temp = copymolecule( m1 );
NAB_matcpy( mat_x, newtransform( 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, re_arr * 3.428000E+01 ) );
transformmol( mat_x, m_temp, NULL );

mergestr( m_final, "D", "last", m_temp, "A", "first" );

m_temp = copymolecule( m1 );
NAB_matcpy( mat_x, newtransform( 2 * disp, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, re_arr * 3.428000E+01 ) );
transformmol( mat_x, m_temp, NULL );

mergestr( m_final, "E", "last", m_temp, "A", "first" );

m_temp = copymolecule( m1 );
NAB_matcpy( mat_x, newtransform( 4 * disp, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, re_arr * 3.428000E+01 ) );
transformmol( mat_x, m_temp, NULL );

mergestr( m_final, "F", "last", m_temp, "A", "first" );





freemolecule( m1 );
m1 = newmolecule(  );
addstrand( m1, "A" );
for( i = 1;i <= 14;i = i + 1 ){

n = NULL;
NAB_strcpy(  &s, "c" );
if( i <= 7 ){
trise = ( 57 - i ) * 3.380000E+00;
t = ( 57 - i );
ttwist = t * 3.428000E+01;

NAB_matcpy( mat_x, newtransform( 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, re_arr * 3.428000E+01 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
if( i == 1 ){
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &s, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "a5" ) );
}
else{
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &s, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}
transformmol( mat_drz, m, NULL );
transformmol( mat_x, m, NULL );
mergestr( m1, "A", "last", m, "anti", "first" );
}
else{
trise = ( 42 + i ) * 3.380000E+00;
t = ( 42 + i );
ttwist = t * 3.428000E+01;
NAB_matcpy( mat_x, newtransform( disp, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, re_arr * 3.428000E+01 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
if( i == 14 ){
m = wc_helix(  &s, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &n, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "s3" ) );
}
else{
m = wc_helix(  &s, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &n, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}
transformmol( mat_drz, m, NULL );
transformmol( mat_x, m, NULL );
mergestr( m1, "A", "last", m, "sense", "first" );

}
freemolecule( m );
n = NULL;
NAB_strcpy(  &s, "c" );
printf( "%s %s %d %d %lf\n", s, n, i, t, ttwist );
}

mergestr( m_final, "G", "last", m1, "A", "first" );





freemolecule( m1 );
m1 = newmolecule(  );
addstrand( m1, "A" );
for( i = 1;i <= 14;i = i + 1 ){

n = NULL;
NAB_strcpy(  &s, "c" );
if( i <= 7 ){
trise = ( 57 - i ) * 3.380000E+00;
t = ( 57 - i );
ttwist = t * 3.428000E+01;

NAB_matcpy( mat_x, newtransform( 2 * disp, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, re_arr * 3.428000E+01 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
if( i == 1 ){
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &s, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "a5" ) );
}
else{
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &s, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}
transformmol( mat_drz, m, NULL );
transformmol( mat_x, m, NULL );
mergestr( m1, "A", "last", m, "anti", "first" );
}
else{
trise = ( 42 + i ) * 3.380000E+00;
t = ( 42 + i );
ttwist = t * 3.428000E+01;
NAB_matcpy( mat_x, newtransform( 3 * disp, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, re_arr * 3.428000E+01 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
if( i == 14 ){
m = wc_helix(  &s, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &n, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "s3" ) );
}
else{
m = wc_helix(  &s, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &n, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}
transformmol( mat_drz, m, NULL );
transformmol( mat_x, m, NULL );
mergestr( m1, "A", "last", m, "sense", "first" );

}
freemolecule( m );
n = NULL;
NAB_strcpy(  &s, "c" );
printf( "%s %s %d %d %lf\n", s, n, i, t, ttwist );
}

mergestr( m_final, "H", "last", m1, "A", "first" );




freemolecule( m1 );
m1 = newmolecule(  );
addstrand( m1, "A" );
for( i = 1;i <= 14;i = i + 1 ){

n = NULL;
NAB_strcpy(  &s, "c" );
if( i <= 7 ){
trise = ( 57 - i ) * 3.380000E+00;
t = ( 57 - i );
ttwist = t * 3.428000E+01;

NAB_matcpy( mat_x, newtransform( 4 * disp, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, re_arr * 3.428000E+01 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
if( i == 1 ){
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &s, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "a5" ) );
}
else{
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &s, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}
transformmol( mat_drz, m, NULL );
transformmol( mat_x, m, NULL );
mergestr( m1, "A", "last", m, "anti", "first" );
}
else{
trise = ( 42 + i ) * 3.380000E+00;
t = ( 42 + i );
ttwist = t * 3.428000E+01;
NAB_matcpy( mat_x, newtransform( 5 * disp, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, re_arr * 3.428000E+01 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
if( i == 14 ){
m = wc_helix(  &s, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &n, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "s3" ) );
}
else{
m = wc_helix(  &s, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &n, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}
transformmol( mat_drz, m, NULL );
transformmol( mat_x, m, NULL );
mergestr( m1, "A", "last", m, "sense", "first" );

}
freemolecule( m );
n = NULL;
NAB_strcpy(  &s, "c" );
printf( "%s %s %d %d %lf\n", s, n, i, t, ttwist );
}

mergestr( m_final, "I", "last", m1, "A", "first" );





freemolecule( m1 );
m1 = newmolecule(  );
addstrand( m1, "A" );
for( i = 1;i <= 14;i = i + 1 ){

n = NULL;
NAB_strcpy(  &s, "c" );
if( i <= 7 ){
trise = ( i ) * 3.380000E+00;
t = ( i );
ttwist = t * 3.428000E+01;

NAB_matcpy( mat_x, newtransform( disp, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, re_arr * 3.428000E+01 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
if( i == 1 ){
m = wc_helix(  &s, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &n, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "s5" ) );
}
else{
m = wc_helix(  &s, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &n, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}
transformmol( mat_drz, m, NULL );
transformmol( mat_x, m, NULL );
mergestr( m1, "A", "last", m, "sense", "first" );
}
else{
trise = ( 15 - i ) * 3.380000E+00;
t = ( 15 - i );
ttwist = t * 3.428000E+01;
NAB_matcpy( mat_x, newtransform( 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, re_arr * 3.428000E+01 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
if( i == 14 ){
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &s, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "a3" ) );
}
else{
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &s, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}
transformmol( mat_drz, m, NULL );
transformmol( mat_x, m, NULL );
mergestr( m1, "A", "last", m, "anti", "first" );

}
freemolecule( m );
n = NULL;
NAB_strcpy(  &s, "c" );
printf( "%s %s %d %d %lf\n", s, n, i, t, ttwist );
}

mergestr( m_final, "J", "last", m1, "A", "first" );





freemolecule( m1 );
m1 = newmolecule(  );
addstrand( m1, "A" );
for( i = 1;i <= 14;i = i + 1 ){

n = NULL;
NAB_strcpy(  &s, "c" );
if( i <= 7 ){
trise = ( i ) * 3.380000E+00;
t = ( i );
ttwist = t * 3.428000E+01;

NAB_matcpy( mat_x, newtransform( 3 * disp, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, re_arr * 3.428000E+01 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
if( i == 1 ){
m = wc_helix(  &s, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &n, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "s5" ) );
}
else{
m = wc_helix(  &s, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &n, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}
transformmol( mat_drz, m, NULL );
transformmol( mat_x, m, NULL );
mergestr( m1, "A", "last", m, "sense", "first" );
}
else{
trise = ( 15 - i ) * 3.380000E+00;
t = ( 15 - i );
ttwist = t * 3.428000E+01;
NAB_matcpy( mat_x, newtransform( 2 * disp, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, re_arr * 3.428000E+01 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
if( i == 14 ){
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &s, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "a3" ) );
}
else{
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &s, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}
transformmol( mat_drz, m, NULL );
transformmol( mat_x, m, NULL );
mergestr( m1, "A", "last", m, "anti", "first" );

}
freemolecule( m );
n = NULL;
NAB_strcpy(  &s, "c" );
printf( "%s %s %d %d %lf\n", s, n, i, t, ttwist );
}

mergestr( m_final, "K", "last", m1, "A", "first" );





freemolecule( m1 );
m1 = newmolecule(  );
addstrand( m1, "A" );
for( i = 1;i <= 14;i = i + 1 ){

n = NULL;
NAB_strcpy(  &s, "c" );
if( i <= 7 ){
trise = ( i ) * 3.380000E+00;
t = ( i );
ttwist = t * 3.428000E+01;

NAB_matcpy( mat_x, newtransform( 5 * disp, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, re_arr * 3.428000E+01 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
if( i == 1 ){
m = wc_helix(  &s, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &n, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "s5" ) );
}
else{
m = wc_helix(  &s, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &n, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}
transformmol( mat_drz, m, NULL );
transformmol( mat_x, m, NULL );
mergestr( m1, "A", "last", m, "sense", "first" );
}
else{
trise = ( 15 - i ) * 3.380000E+00;
t = ( 15 - i );
ttwist = t * 3.428000E+01;
NAB_matcpy( mat_x, newtransform( 4 * disp, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, re_arr * 3.428000E+01 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
if( i == 14 ){
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &s, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "a3" ) );
}
else{
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &s, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}
transformmol( mat_drz, m, NULL );
transformmol( mat_x, m, NULL );
mergestr( m1, "A", "last", m, "anti", "first" );

}
freemolecule( m );
n = NULL;
NAB_strcpy(  &s, "c" );
printf( "%s %s %d %d %lf\n", s, n, i, t, ttwist );
}

mergestr( m_final, "L", "last", m1, "A", "first" );




freemolecule( m1 );
m1 = newmolecule(  );
addstrand( m1, "A" );
for( i = 1;i <= 25;i = i + 1 ){

n = NULL;
NAB_strcpy(  &s, "c" );
if( i <= 14 ){
trise = ( 22 - i ) * 3.380000E+00;
t = ( 22 - i );
ttwist = t * 3.428000E+01;

NAB_matcpy( mat_x, newtransform( 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, re_arr * 3.428000E+01 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
if( i == 1 ){
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &s, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "a5" ) );
}
else{
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &s, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}
transformmol( mat_drz, m, NULL );
transformmol( mat_x, m, NULL );
mergestr( m1, "A", "last", m, "anti", "first" );
}

else if( i > 14 && i <= 21 ){
trise = ( i - 7 ) * 3.380000E+00;
t = ( i - 7 );
ttwist = t * 3.428000E+01;

NAB_matcpy( mat_x, newtransform( disp, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, re_arr * 3.428000E+01 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
if( i == 15 ){
m = wc_helix(  &s, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &n, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}
else{
m = wc_helix(  &s, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &n, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}
transformmol( mat_drz, m, NULL );
transformmol( mat_x, m, NULL );
mergestr( m1, "A", "last", m, "sense", "first" );
}

else{
trise = ( 36 - i ) * 3.380000E+00;
t = ( 36 - i );
ttwist = t * 3.428000E+01;
NAB_matcpy( mat_x, newtransform( 2 * disp, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, re_arr * 3.428000E+01 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
if( i == 25 ){
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &s, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "a3" ) );
}
else{
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &s, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}
transformmol( mat_drz, m, NULL );
transformmol( mat_x, m, NULL );
mergestr( m1, "A", "last", m, "anti", "first" );

}
freemolecule( m );
n = NULL;
NAB_strcpy(  &s, "c" );
printf( "%s %s %d %d %lf\n", s, n, i, t, ttwist );
}

mergestr( m_final, "M", "last", m1, "A", "first" );






freemolecule( m1 );
m1 = newmolecule(  );
addstrand( m1, "A" );
for( i = 1;i <= 25;i = i + 1 ){

n = NULL;
NAB_strcpy(  &s, "c" );
if( i <= 14 ){
trise = ( 43 - i ) * 3.380000E+00;
t = ( 43 - i );
ttwist = t * 3.428000E+01;

NAB_matcpy( mat_x, newtransform( 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, re_arr * 3.428000E+01 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
if( i == 1 ){
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &s, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "a5" ) );
}
else{
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &s, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}
transformmol( mat_drz, m, NULL );
transformmol( mat_x, m, NULL );
mergestr( m1, "A", "last", m, "anti", "first" );
}

else if( i > 14 && i <= 21 ){
trise = ( i + 14 ) * 3.380000E+00;
t = ( i + 14 );
ttwist = t * 3.428000E+01;

NAB_matcpy( mat_x, newtransform( disp, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, re_arr * 3.428000E+01 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
if( i == 15 ){
m = wc_helix(  &s, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &n, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}
else{
m = wc_helix(  &s, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &n, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}
transformmol( mat_drz, m, NULL );
transformmol( mat_x, m, NULL );
mergestr( m1, "A", "last", m, "sense", "first" );
}

else{
trise = ( 57 - i ) * 3.380000E+00;
t = ( 57 - i );
ttwist = t * 3.428000E+01;
NAB_matcpy( mat_x, newtransform( 2 * disp, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, re_arr * 3.428000E+01 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
if( i == 25 ){
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &s, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "a3" ) );
}
else{
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &s, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}
transformmol( mat_drz, m, NULL );
transformmol( mat_x, m, NULL );
mergestr( m1, "A", "last", m, "anti", "first" );

}
freemolecule( m );
n = NULL;
NAB_strcpy(  &s, "c" );
printf( "%s %s %d %d %lf\n", s, n, i, t, ttwist );
}

mergestr( m_final, "N", "last", m1, "A", "first" );





freemolecule( m1 );
m1 = newmolecule(  );
addstrand( m1, "A" );
for( i = 1;i <= 28;i = i + 1 ){

n = NULL;
NAB_strcpy(  &s, "c" );
if( i <= 7 ){
trise = ( 22 - i ) * 3.380000E+00;
t = ( 22 - i );
ttwist = t * 3.428000E+01;

NAB_matcpy( mat_x, newtransform( 2 * disp, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, re_arr * 3.428000E+01 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
if( i == 1 ){
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &s, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "a5" ) );
}
else{
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &s, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}
transformmol( mat_drz, m, NULL );
transformmol( mat_x, m, NULL );
mergestr( m1, "A", "last", m, "anti", "first" );
}

else if( i > 7 && i <= 21 ){
trise = ( i + 7 ) * 3.380000E+00;
t = ( i + 7 );
ttwist = t * 3.428000E+01;

NAB_matcpy( mat_x, newtransform( disp, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, re_arr * 3.428000E+01 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
if( i == 15 ){
m = wc_helix(  &s, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &n, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}
else{
m = wc_helix(  &s, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &n, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}
transformmol( mat_drz, m, NULL );
transformmol( mat_x, m, NULL );
mergestr( m1, "A", "last", m, "sense", "first" );
}

else{
trise = ( 50 - i ) * 3.380000E+00;
t = ( 50 - i );
ttwist = t * 3.428000E+01;
NAB_matcpy( mat_x, newtransform( 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, re_arr * 3.428000E+01 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
if( i == 28 ){
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &s, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "a3" ) );
}
else{
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &s, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}
transformmol( mat_drz, m, NULL );
transformmol( mat_x, m, NULL );
mergestr( m1, "A", "last", m, "anti", "first" );

}
freemolecule( m );
n = NULL;
NAB_strcpy(  &s, "c" );
printf( "%s %s %d %d %lf\n", s, n, i, t, ttwist );
}

mergestr( m_final, "O", "last", m1, "A", "first" );






freemolecule( m1 );
m1 = newmolecule(  );
addstrand( m1, "A" );
for( i = 1;i <= 28;i = i + 1 ){

n = NULL;
NAB_strcpy(  &s, "c" );
if( i <= 7 ){
trise = ( 22 - i ) * 3.380000E+00;
t = ( 22 - i );
ttwist = t * 3.428000E+01;

NAB_matcpy( mat_x, newtransform( 4 * disp, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, re_arr * 3.428000E+01 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
if( i == 1 ){
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &s, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "a5" ) );
}
else{
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &s, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}
transformmol( mat_drz, m, NULL );
transformmol( mat_x, m, NULL );
mergestr( m1, "A", "last", m, "anti", "first" );
}

else if( i > 7 && i <= 21 ){
trise = ( i + 7 ) * 3.380000E+00;
t = ( i + 7 );
ttwist = t * 3.428000E+01;

NAB_matcpy( mat_x, newtransform( 3 * disp, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, re_arr * 3.428000E+01 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
if( i == 15 ){
m = wc_helix(  &s, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &n, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}
else{
m = wc_helix(  &s, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &n, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}
transformmol( mat_drz, m, NULL );
transformmol( mat_x, m, NULL );
mergestr( m1, "A", "last", m, "sense", "first" );
}

else{
trise = ( 50 - i ) * 3.380000E+00;
t = ( 50 - i );
ttwist = t * 3.428000E+01;
NAB_matcpy( mat_x, newtransform( 2 * disp, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, re_arr * 3.428000E+01 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
if( i == 28 ){
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &s, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "a3" ) );
}
else{
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &s, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}
transformmol( mat_drz, m, NULL );
transformmol( mat_x, m, NULL );
mergestr( m1, "A", "last", m, "anti", "first" );

}
freemolecule( m );
n = NULL;
NAB_strcpy(  &s, "c" );
printf( "%s %s %d %d %lf\n", s, n, i, t, ttwist );
}

mergestr( m_final, "P", "last", m1, "A", "first" );




freemolecule( m1 );
m1 = newmolecule(  );
addstrand( m1, "A" );
for( i = 1;i <= 28;i = i + 1 ){

n = NULL;
NAB_strcpy(  &s, "c" );
if( i <= 7 ){
trise = ( 43 - i ) * 3.380000E+00;
t = ( 43 - i );
ttwist = t * 3.428000E+01;

NAB_matcpy( mat_x, newtransform( 2 * disp, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, re_arr * 3.428000E+01 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
if( i == 1 ){
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &s, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "a5" ) );
}
else{
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &s, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}
transformmol( mat_drz, m, NULL );
transformmol( mat_x, m, NULL );
mergestr( m1, "A", "last", m, "anti", "first" );
}

else if( i > 7 && i <= 21 ){
trise = ( i + 28 ) * 3.380000E+00;
t = ( i + 28 );
ttwist = t * 3.428000E+01;

NAB_matcpy( mat_x, newtransform( disp, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, re_arr * 3.428000E+01 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
if( i == 15 ){
m = wc_helix(  &s, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &n, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}
else{
m = wc_helix(  &s, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &n, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}
transformmol( mat_drz, m, NULL );
transformmol( mat_x, m, NULL );
mergestr( m1, "A", "last", m, "sense", "first" );
}

else{
trise = ( 71 - i ) * 3.380000E+00;
t = ( 71 - i );
ttwist = t * 3.428000E+01;
NAB_matcpy( mat_x, newtransform( 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, re_arr * 3.428000E+01 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
if( i == 28 ){
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &s, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "a3" ) );
}
else{
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &s, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}
transformmol( mat_drz, m, NULL );
transformmol( mat_x, m, NULL );
mergestr( m1, "A", "last", m, "anti", "first" );

}
freemolecule( m );
n = NULL;
NAB_strcpy(  &s, "c" );
printf( "%s %s %d %d %lf\n", s, n, i, t, ttwist );
}

mergestr( m_final, "R", "last", m1, "A", "first" );





freemolecule( m1 );
m1 = newmolecule(  );
addstrand( m1, "A" );
for( i = 1;i <= 28;i = i + 1 ){

n = NULL;
NAB_strcpy(  &s, "c" );
if( i <= 7 ){
trise = ( 43 - i ) * 3.380000E+00;
t = ( 43 - i );
ttwist = t * 3.428000E+01;

NAB_matcpy( mat_x, newtransform( 4 * disp, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, re_arr * 3.428000E+01 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
if( i == 1 ){
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &s, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "a5" ) );
}
else{
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &s, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}
transformmol( mat_drz, m, NULL );
transformmol( mat_x, m, NULL );
mergestr( m1, "A", "last", m, "anti", "first" );
}

else if( i > 7 && i <= 21 ){
trise = ( i + 28 ) * 3.380000E+00;
t = ( i + 28 );
ttwist = t * 3.428000E+01;

NAB_matcpy( mat_x, newtransform( 3 * disp, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, re_arr * 3.428000E+01 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
if( i == 15 ){
m = wc_helix(  &s, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &n, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}
else{
m = wc_helix(  &s, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &n, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}
transformmol( mat_drz, m, NULL );
transformmol( mat_x, m, NULL );
mergestr( m1, "A", "last", m, "sense", "first" );
}

else{
trise = ( 71 - i ) * 3.380000E+00;
t = ( 71 - i );
ttwist = t * 3.428000E+01;
NAB_matcpy( mat_x, newtransform( 2 * disp, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, re_arr * 3.428000E+01 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
if( i == 28 ){
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &s, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "a3" ) );
}
else{
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &s, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}
transformmol( mat_drz, m, NULL );
transformmol( mat_x, m, NULL );
mergestr( m1, "A", "last", m, "anti", "first" );

}
freemolecule( m );
n = NULL;
NAB_strcpy(  &s, "c" );
printf( "%s %s %d %d %lf\n", s, n, i, t, ttwist );
}

mergestr( m_final, "My", "last", m1, "A", "first" );







freemolecule( m1 );
m1 = newmolecule(  );
addstrand( m1, "A" );
for( i = 1;i <= 24;i = i + 1 ){

n = NULL;
NAB_strcpy(  &s, "c" );
if( i <= 3 ){
trise = ( 11 - i ) * 3.380000E+00;
t = ( 11 - i );
ttwist = t * 3.428000E+01;

NAB_matcpy( mat_x, newtransform( 2 * disp, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, re_arr * 3.428000E+01 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
if( i == 1 ){
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &s, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "a5" ) );
}
else{
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &s, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}
transformmol( mat_drz, m, NULL );
transformmol( mat_x, m, NULL );
mergestr( m1, "A", "last", m, "anti", "first" );
}

else if( i > 3 && i <= 10 ){
trise = ( i + 4 ) * 3.380000E+00;
t = ( i + 4 );
ttwist = t * 3.428000E+01;

NAB_matcpy( mat_x, newtransform( 3 * disp, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, re_arr * 3.428000E+01 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
if( i == 8 ){
m = wc_helix(  &s, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &n, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}
else{
m = wc_helix(  &s, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &n, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}
transformmol( mat_drz, m, NULL );
transformmol( mat_x, m, NULL );
mergestr( m1, "A", "last", m, "sense", "first" );
}
else if( i > 10 && i <= 17 ){
trise = ( 25 - i ) * 3.380000E+00;
t = ( 25 - i );
ttwist = t * 3.428000E+01;

NAB_matcpy( mat_x, newtransform( 4 * disp, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, re_arr * 3.428000E+01 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
if( i == 15 ){
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &s, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}
else{
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &s, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}
transformmol( mat_drz, m, NULL );
transformmol( mat_x, m, NULL );
mergestr( m1, "A", "last", m, "anti", "first" );
}
else{
trise = ( i - 10 ) * 3.380000E+00;
t = ( i - 10 );
ttwist = t * 3.428000E+01;
NAB_matcpy( mat_x, newtransform( 5 * disp, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, re_arr * 3.428000E+01 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
if( i == 24 ){
m = wc_helix(  &s, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &n, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "s3" ) );
}
else{
m = wc_helix(  &s, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &n, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}
transformmol( mat_drz, m, NULL );
transformmol( mat_x, m, NULL );
mergestr( m1, "A", "last", m, "sense", "first" );

}
freemolecule( m );
n = NULL;
NAB_strcpy(  &s, "c" );
printf( "%s %s %d %d %lf\n", s, n, i, t, ttwist );
}

mergestr( m_final, "Ny", "last", m1, "A", "first" );




freemolecule( m1 );
m1 = newmolecule(  );
addstrand( m1, "A" );
for( i = 1;i <= 24;i = i + 1 ){

n = NULL;
NAB_strcpy(  &s, "c" );
if( i <= 3 ){
trise = ( 32 - i ) * 3.380000E+00;
t = ( 32 - i );
ttwist = t * 3.428000E+01;

NAB_matcpy( mat_x, newtransform( 2 * disp, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, re_arr * 3.428000E+01 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
if( i == 1 ){
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &s, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "a5" ) );
}
else{
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &s, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}
transformmol( mat_drz, m, NULL );
transformmol( mat_x, m, NULL );
mergestr( m1, "A", "last", m, "anti", "first" );
}

else if( i > 3 && i <= 10 ){
trise = ( i + 25 ) * 3.380000E+00;
t = ( i + 25 );
ttwist = t * 3.428000E+01;

NAB_matcpy( mat_x, newtransform( 3 * disp, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, re_arr * 3.428000E+01 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
if( i == 8 ){
m = wc_helix(  &s, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &n, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}
else{
m = wc_helix(  &s, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &n, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}
transformmol( mat_drz, m, NULL );
transformmol( mat_x, m, NULL );
mergestr( m1, "A", "last", m, "sense", "first" );
}
else if( i > 10 && i <= 17 ){
trise = ( 46 - i ) * 3.380000E+00;
t = ( 46 - i );
ttwist = t * 3.428000E+01;

NAB_matcpy( mat_x, newtransform( 4 * disp, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, re_arr * 3.428000E+01 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
if( i == 15 ){
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &s, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}
else{
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &s, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}
transformmol( mat_drz, m, NULL );
transformmol( mat_x, m, NULL );
mergestr( m1, "A", "last", m, "anti", "first" );
}
else{
trise = ( i + 11 ) * 3.380000E+00;
t = ( i + 11 );
ttwist = t * 3.428000E+01;
NAB_matcpy( mat_x, newtransform( 5 * disp, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, re_arr * 3.428000E+01 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
if( i == 24 ){
m = wc_helix(  &s, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &n, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "s3" ) );
}
else{
m = wc_helix(  &s, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &n, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}
transformmol( mat_drz, m, NULL );
transformmol( mat_x, m, NULL );
mergestr( m1, "A", "last", m, "sense", "first" );

}
freemolecule( m );
n = NULL;
NAB_strcpy(  &s, "c" );
printf( "%s %s %d %d %lf\n", s, n, i, t, ttwist );
}

mergestr( m_final, "Oy", "last", m1, "A", "first" );






freemolecule( m1 );
m1 = newmolecule(  );
addstrand( m1, "A" );
for( i = 1;i <= 21;i = i + 1 ){

n = NULL;
NAB_strcpy(  &s, "c" );
if( i <= 14 ){
trise = ( i + 14 ) * 3.380000E+00;
t = ( i + 14 );
ttwist = t * 3.428000E+01;

NAB_matcpy( mat_x, newtransform( 5 * disp, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, re_arr * 3.428000E+01 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
if( i == 1 ){
m = wc_helix(  &s, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &n, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "s5" ) );
}
else{
m = wc_helix(  &s, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &n, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}
transformmol( mat_drz, m, NULL );
transformmol( mat_x, m, NULL );
mergestr( m1, "A", "last", m, "sense", "first" );
}
else{
trise = ( 43 - i ) * 3.380000E+00;
t = ( 43 - i );
ttwist = t * 3.428000E+01;
NAB_matcpy( mat_x, newtransform( 4 * disp, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, re_arr * 3.428000E+01 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
if( i == 21 ){
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &s, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "a3" ) );
}
else{
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &s, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}
transformmol( mat_drz, m, NULL );
transformmol( mat_x, m, NULL );
mergestr( m1, "A", "last", m, "anti", "first" );

}
freemolecule( m );
n = NULL;
NAB_strcpy(  &s, "c" );
printf( "%s %s %d %d %lf\n", s, n, i, t, ttwist );
}

mergestr( m_final, "Py", "last", m1, "A", "first" );









freemolecule( m1 );
m1 = newmolecule(  );
addstrand( m1, "A" );
for( i = 1;i <= 21;i = i + 1 ){

n = NULL;
NAB_strcpy(  &s, "c" );
if( i <= 14 ){
trise = ( i + 35 ) * 3.380000E+00;
t = ( i + 35 );
ttwist = t * 3.428000E+01;

NAB_matcpy( mat_x, newtransform( 5 * disp, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, re_arr * 3.428000E+01 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
if( i == 1 ){
m = wc_helix(  &s, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &n, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "s5" ) );
}
else{
m = wc_helix(  &s, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &n, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}
transformmol( mat_drz, m, NULL );
transformmol( mat_x, m, NULL );
mergestr( m1, "A", "last", m, "sense", "first" );
}
else{
trise = ( 64 - i ) * 3.380000E+00;
t = ( 64 - i );
ttwist = t * 3.428000E+01;
NAB_matcpy( mat_x, newtransform( 4 * disp, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, re_arr * 3.428000E+01 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
if( i == 21 ){
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &s, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "a3" ) );
}
else{
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &s, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}
transformmol( mat_drz, m, NULL );
transformmol( mat_x, m, NULL );
mergestr( m1, "A", "last", m, "anti", "first" );

}
freemolecule( m );
n = NULL;
NAB_strcpy(  &s, "c" );
printf( "%s %s %d %d %lf\n", s, n, i, t, ttwist );
}

mergestr( m_final, "Qy", "last", m1, "A", "first" );


putpdb( "dna_sheet_exp.pdb", m_final, NULL );


	exit( 0 );
}
