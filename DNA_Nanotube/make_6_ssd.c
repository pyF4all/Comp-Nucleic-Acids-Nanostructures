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

m1 = newmolecule(  );
addstrand( m1, "A" );


for( i = 1;i <= 59;i = i + 1 ){
trise = ( i ) * 3.380000E+00;
ttwist = ( i ) * 3.428000E+01;

NAB_strcpy(  &s, "a" );
n = NULL;

NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );

printf( "%s %s\n", s, n );
if( i == 1 ){
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &s, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "a3" ) );
}
else if( i == 59 ){
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &s, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "a5" ) );
}
else{
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &s, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}

transformmol( mat_drz, m, NULL );
mergestr( m1, "A", "first", m, "anti", "last" );
freemolecule( m );
n = NULL;
NAB_strcpy(  &s, "a" );
}

m_temp = copymolecule( m1 );
NAB_matcpy( mat_x, newtransform( 0.000000E+00, 2.000000E+01, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00 ) );
transformmol( mat_x, m_temp, NULL );
putpdb( "dna1.pdb", m_temp, NULL );
mergestr( m_final, "S", "last", m_temp, "A", "first" );

m_temp = copymolecule( m1 );
NAB_matcpy( mat_x, newtransform( 1.732000E+01,  - 1.000000E+01, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00 ) );
transformmol( mat_x, m_temp, NULL );
putpdb( "dna3.pdb", m_temp, NULL );
mergestr( m_final, "B", "last", m_temp, "A", "first" );

m_temp = copymolecule( m1 );
NAB_matcpy( mat_x, newtransform(  - 1.732000E+01,  - 1.000000E+01, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00 ) );
transformmol( mat_x, m_temp, NULL );
mergestr( m_final, "C", "last", m_temp, "A", "first" );
putpdb( "dna5.pdb", m_temp, NULL );

freemolecule( m_temp );
freemolecule( m1 );
m1 = newmolecule(  );
addstrand( m1, "A" );
n = NULL;
NAB_strcpy(  &s, "a" );

for( i = 1;i <= 59;i = i + 1 ){
trise = ( i ) * 3.380000E+00;
ttwist = ( i ) * 3.428000E+01;

NAB_strcpy(  &s, "a" );
n = NULL;
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );


printf( "%s %s\n", s, n );
if( i == 1 ){
m = wc_helix(  &s, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &n, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "s5" ) );
}
else if( i == 59 ){
m = wc_helix(  &s, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &n, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "s3" ) );
}
else{
m = wc_helix(  &s, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &n, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
}

transformmol( mat_drz, m, NULL );
mergestr( m1, "A", "last", m, "sense", "first" );
freemolecule( m );
n = NULL;
NAB_strcpy(  &s, "a" );
}
m_temp = copymolecule( m1 );
NAB_matcpy( mat_x, newtransform( 1.732000E+01, 1.000000E+01, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00 ) );
transformmol( mat_x, m_temp, NULL );
putpdb( "dna2.pdb", m_temp, NULL );
mergestr( m_final, "D", "last", m_temp, "A", "first" );

m_temp = copymolecule( m1 );
NAB_matcpy( mat_x, newtransform( 0.000000E+00,  - 2.000000E+01, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00 ) );
transformmol( mat_x, m_temp, NULL );
putpdb( "dna4.pdb", m_temp, NULL );
mergestr( m_final, "E", "last", m_temp, "A", "first" );

m_temp = copymolecule( m1 );
NAB_matcpy( mat_x, newtransform(  - 1.732000E+01, 1.000000E+01, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00 ) );
transformmol( mat_x, m_temp, NULL );
putpdb( "dna6.pdb", m_temp, NULL );
mergestr( m_final, "F", "last", m_temp, "A", "first" );
putpdb( "6h.pdb", m_final, NULL );



freemolecule( m1 );
m1 = newmolecule(  );
addstrand( m1, "A" );
for( i = 1;i <= 24;i = i + 1 ){

n = NULL;
NAB_strcpy(  &s, "t" );
if( i <= 12 ){
trise = ( i ) * 3.380000E+00;
t = ( i );
ttwist = t * 3.428000E+01;

NAB_matcpy( mat_x, newtransform( 0.000000E+00, 2.000000E+01, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00 ) );
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
trise = ( 25 - i ) * 3.380000E+00;
t = ( 25 - i );
ttwist = t * 3.428000E+01;
NAB_matcpy( mat_x, newtransform( 1.732000E+01, 1.000000E+01, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
if( i == 24 ){
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
NAB_strcpy(  &s, "t" );
printf( "%s %s %d %d %lf\n", s, n, i, t, ttwist );
}
putpdb( "ssd1.pdb", m1, NULL );
mergestr( m_final, "G", "last", m1, "A", "first" );
putpdb( "6h1.pdb", m_final, NULL );


freemolecule( m1 );
m1 = newmolecule(  );
addstrand( m1, "A" );
for( i = 1;i <= 24;i = i + 1 ){
n = NULL;
NAB_strcpy(  &s, "t" );
if( i <= 12 ){
trise = ( i ) * 3.380000E+00;
t = ( i );
ttwist = t * 3.428000E+01;
NAB_matcpy( mat_x, newtransform(  - 1.732000E+01,  - 1.000000E+01, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00 ) );
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
trise = ( 25 - i ) * 3.380000E+00;
t = ( 25 - i );
ttwist = t * 3.428000E+01;
NAB_matcpy( mat_x, newtransform( 0.000000E+00,  - 2.000000E+01, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
if( i == 24 ){
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
NAB_strcpy(  &s, "t" );
printf( "%s %s %d %d %lf\n", s, n, i, t, ttwist );
}
putpdb( "ssd9.pdb", m1, NULL );
mergestr( m_final, "H", "last", m1, "A", "first" );
putpdb( "6h2.pdb", m_final, NULL );

freemolecule( m1 );
m1 = newmolecule(  );
addstrand( m1, "A" );
for( i = 1;i <= 24;i = i + 1 ){

n = NULL;
NAB_strcpy(  &s, "t" );
if( i <= 12 ){
t = ( 60 - i );
trise = t * 3.380000E+00;
ttwist = t * 3.428000E+01;
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
NAB_matcpy( mat_x, newtransform( 1.732000E+01, 1.000000E+01, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00 ) );
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
t = ( 35 + i );
trise = t * 3.380000E+00;
ttwist = t * 3.428000E+01;
NAB_matcpy( mat_x, newtransform( 1.732000E+01,  - 1.000000E+01, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00 ) );
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
NAB_strcpy(  &s, "t" );
printf( "%s %s %d %d %lf\n", s, n, i, t, ttwist );
}
putpdb( "ssd6.pdb", m1, NULL );
mergestr( m_final, "I", "last", m1, "A", "first" );
putpdb( "6h3.pdb", m_final, NULL );

freemolecule( m1 );
m1 = newmolecule(  );
addstrand( m1, "A" );
for( i = 1;i <= 24;i = i + 1 ){
n = NULL;
NAB_strcpy(  &s, "t" );
if( i <= 12 ){
t = ( 60 - i );
trise = t * 3.380000E+00;
ttwist = t * 3.428000E+01;
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
NAB_matcpy( mat_x, newtransform(  - 1.732000E+01, 1.000000E+01, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00 ) );
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
t = ( 35 + i );
trise = t * 3.380000E+00;
ttwist = t * 3.428000E+01;
NAB_matcpy( mat_x, newtransform(  - 1.732000E+01,  - 1.000000E+01, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
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
NAB_strcpy(  &s, "t" );
printf( "%s %s %d %d %lf\n", s, n, i, t, ttwist );
}
putpdb( "ssd12.pdb", m1, NULL );
mergestr( m_final, "J", "last", m1, "A", "first" );
putpdb( "6h4.pdb", m_final, NULL );

freemolecule( m1 );
m1 = newmolecule(  );
addstrand( m1, "A" );
for( i = 1;i <= 33;i = i + 1 ){
n = NULL;
NAB_strcpy(  &s, "t" );
if( i <= 19 ){
t = ( i );
trise = t * 3.380000E+00;
ttwist = t * 3.428000E+01;
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
NAB_matcpy( mat_x, newtransform( 1.732000E+01,  - 1.000000E+01, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00 ) );
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
else if( i > 19 && i <= 26 ){
t = ( 39 - i );
trise = t * 3.380000E+00;
ttwist = t * 3.428000E+01;
NAB_matcpy( mat_x, newtransform( 0.000000E+00,  - 2.000000E+01, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &s, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
transformmol( mat_drz, m, NULL );
transformmol( mat_x, m, NULL );
mergestr( m1, "A", "last", m, "anti", "first" );
}
else{
t = ( i - 14 );
trise = t * 3.380000E+00;
ttwist = t * 3.428000E+01;
NAB_matcpy( mat_x, newtransform(  - 1.732000E+01,  - 1.000000E+01, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
if( i == 33 ){
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
NAB_strcpy(  &s, "t" );
printf( "%s %s %d %d %lf\n", s, n, i, t, ttwist );
}
putpdb( "ssd4A.pdb", m1, NULL );
mergestr( m_final, "K", "last", m1, "A", "first" );
putpdb( "6h5.pdb", m_final, NULL );

freemolecule( m1 );
m1 = newmolecule(  );
addstrand( m1, "A" );
for( i = 1;i <= 33;i = i + 1 ){
n = NULL;
NAB_strcpy(  &s, "t" );
if( i <= 19 ){
t = ( 60 - i );
trise = t * 3.380000E+00;
ttwist = t * 3.428000E+01;
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
NAB_matcpy( mat_x, newtransform( 0.000000E+00,  - 2.000000E+01, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00 ) );
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
else if( i > 19 && i <= 26 ){
t = ( 21 + i );
trise = t * 3.380000E+00;
ttwist = t * 3.428000E+01;
NAB_matcpy( mat_x, newtransform( 1.732000E+01,  - 1.000000E+01, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );

m = wc_helix(  &s, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &n, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );

transformmol( mat_drz, m, NULL );
transformmol( mat_x, m, NULL );
mergestr( m1, "A", "last", m, "sense", "first" );
}
else{
t = ( 74 - i );
trise = t * 3.380000E+00;
ttwist = t * 3.428000E+01;
NAB_matcpy( mat_x, newtransform( 1.732000E+01, 1.000000E+01, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
if( i == 33 ){
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
NAB_strcpy(  &s, "t" );
printf( "%s %s %d %d %lf\n", s, n, i, t, ttwist );
}
putpdb( "ssd3A.pdb", m1, NULL );
mergestr( m_final, "L", "last", m1, "A", "first" );
putpdb( "6h6.pdb", m_final, NULL );

freemolecule( m1 );
m1 = newmolecule(  );
addstrand( m1, "A" );
for( i = 1;i <= 42;i = i + 1 ){
n = NULL;
NAB_strcpy(  &s, "t" );
if( i <= 7 ){
t = ( 19 + i );
trise = t * 3.380000E+00;
ttwist = t * 3.428000E+01;
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
NAB_matcpy( mat_x, newtransform(  - 1.732000E+01,  - 1.000000E+01, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00 ) );
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
else if( i > 7 && i <= 14 ){
t = ( 34 - i );
trise = t * 3.380000E+00;
ttwist = t * 3.428000E+01;
NAB_matcpy( mat_x, newtransform(  - 1.732000E+01, 1.000000E+01, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &s, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );

transformmol( mat_drz, m, NULL );
transformmol( mat_x, m, NULL );
mergestr( m1, "A", "last", m, "anti", "first" );
}
else if( i > 14 && i <= 28 ){
t = ( i + 5 );
trise = t * 3.380000E+00;
ttwist = t * 3.428000E+01;
NAB_matcpy( mat_x, newtransform( 0.000000E+00, 2.000000E+01, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
m = wc_helix(  &s, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &n, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );

transformmol( mat_drz, m, NULL );
transformmol( mat_x, m, NULL );
mergestr( m1, "A", "last", m, "sense", "first" );
}
else if( i > 28 && i <= 35 ){
t = ( 62 - i );
trise = t * 3.380000E+00;
ttwist = t * 3.428000E+01;
NAB_matcpy( mat_x, newtransform( 1.732000E+01, 1.000000E+01, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &s, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );

transformmol( mat_drz, m, NULL );
transformmol( mat_x, m, NULL );
mergestr( m1, "A", "last", m, "anti", "first" );
}

else{
t = ( i - 9 );
trise = t * 3.380000E+00;
ttwist = t * 3.428000E+01;
NAB_matcpy( mat_x, newtransform( 1.732000E+01,  - 1.000000E+01, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
if( i == 42 ){
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
NAB_strcpy(  &s, "t" );
printf( "%s %s %d %d %lf\n", s, n, i, t, ttwist );
}
putpdb( "ssd4B.pdb", m1, NULL );
mergestr( m_final, "M", "last", m1, "A", "first" );
putpdb( "6h7.pdb", m_final, NULL );

freemolecule( m1 );
m1 = newmolecule(  );
addstrand( m1, "A" );
for( i = 1;i <= 42;i = i + 1 ){
n = NULL;
NAB_strcpy(  &s, "t" );
if( i <= 7 ){
t = ( 34 - i );
trise = t * 3.380000E+00;
ttwist = t * 3.428000E+01;
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
NAB_matcpy( mat_x, newtransform(  - 1.732000E+01, 1.000000E+01, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00 ) );
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
else if( i > 7 && i <= 14 ){
t = ( 19 + i );
trise = t * 3.380000E+00;
ttwist = t * 3.428000E+01;
NAB_matcpy( mat_x, newtransform(  - 1.732000E+01,  - 1.000000E+01, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
m = wc_helix(  &s, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &n, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );

transformmol( mat_drz, m, NULL );
transformmol( mat_x, m, NULL );
mergestr( m1, "A", "last", m, "sense", "first" );
}
else if( i > 14 && i <= 28 ){
t = ( 48 - i );
trise = t * 3.380000E+00;
ttwist = t * 3.428000E+01;
NAB_matcpy( mat_x, newtransform( 0.000000E+00,  - 2.000000E+01, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &s, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );

transformmol( mat_drz, m, NULL );
transformmol( mat_x, m, NULL );
mergestr( m1, "A", "last", m, "anti", "first" );
}
else if( i > 28 && i <= 35 ){
t = ( i - 9 );
trise = t * 3.380000E+00;
ttwist = t * 3.428000E+01;
NAB_matcpy( mat_x, newtransform( 1.732000E+01,  - 1.000000E+01, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
m = wc_helix(  &s, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &n, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );

transformmol( mat_drz, m, NULL );
transformmol( mat_x, m, NULL );
mergestr( m1, "A", "last", m, "sense", "first" );
}

else{
t = ( 62 - i );
trise = t * 3.380000E+00;
ttwist = t * 3.428000E+01;
NAB_matcpy( mat_x, newtransform( 1.732000E+01, 1.000000E+01, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
if( i == 42 ){
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
NAB_strcpy(  &s, "t" );
printf( "%s %s %d %d %lf\n", s, n, i, t, ttwist );
}
putpdb( "ssd3C.pdb", m1, NULL );
mergestr( m_final, "N", "last", m1, "A", "first" );
putpdb( "6h8.pdb", m_final, NULL );

freemolecule( m1 );
m1 = newmolecule(  );
addstrand( m1, "A" );
for( i = 1;i <= 21;i = i + 1 ){
n = NULL;
NAB_strcpy(  &s, "t" );
if( i <= 7 ){
t = ( 33 + i );
trise = t * 3.380000E+00;
ttwist = t * 3.428000E+01;
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
NAB_matcpy( mat_x, newtransform( 1.732000E+01,  - 1.000000E+01, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00 ) );
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
else if( i > 7 && i <= 14 ){
t = ( 48 - i );
trise = t * 3.380000E+00;
ttwist = t * 3.428000E+01;
NAB_matcpy( mat_x, newtransform( 0.000000E+00,  - 2.000000E+01, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &s, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );

transformmol( mat_drz, m, NULL );
transformmol( mat_x, m, NULL );
mergestr( m1, "A", "last", m, "anti", "first" );
}
else{
t = ( i + 19 );
trise = t * 3.380000E+00;
ttwist = t * 3.428000E+01;
NAB_matcpy( mat_x, newtransform(  - 1.732000E+01,  - 1.000000E+01, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
if( i == 21 ){
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
NAB_strcpy(  &s, "t" );
printf( "%s %s %d %d %lf\n", s, n, i, t, ttwist );
}
putpdb( "ssd4C.pdb", m1, NULL );
mergestr( m_final, "O", "last", m1, "A", "first" );
putpdb( "6h9.pdb", m_final, NULL );

freemolecule( m1 );
m1 = newmolecule(  );
addstrand( m1, "A" );
for( i = 1;i <= 21;i = i + 1 ){
n = NULL;
NAB_strcpy(  &s, "t" );
if( i <= 7 ){
t = ( 41 - i );
trise = t * 3.380000E+00;
ttwist = t * 3.428000E+01;
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
NAB_matcpy( mat_x, newtransform( 1.732000E+01, 1.000000E+01, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00 ) );
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
else if( i > 7 && i <= 14 ){
t = ( i + 26 );
trise = t * 3.380000E+00;
ttwist = t * 3.428000E+01;
NAB_matcpy( mat_x, newtransform( 0.000000E+00, 2.000000E+01, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
m = wc_helix(  &s, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &n, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );

transformmol( mat_drz, m, NULL );
transformmol( mat_x, m, NULL );
mergestr( m1, "A", "last", m, "sense", "first" );
}
else{
t = ( 55 - i );
trise = t * 3.380000E+00;
ttwist = t * 3.428000E+01;
NAB_matcpy( mat_x, newtransform(  - 1.732000E+01, 1.000000E+01, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00 ) );
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
NAB_strcpy(  &s, "t" );
printf( "%s %s %d %d %lf\n", s, n, i, t, ttwist );
}
putpdb( "ssd3B.pdb", m1, NULL );
mergestr( m_final, "P", "last", m1, "A", "first" );
putpdb( "6h10.pdb", m_final, NULL );

freemolecule( m1 );
m1 = newmolecule(  );
addstrand( m1, "A" );
for( i = 1;i <= 33;i = i + 1 ){
n = NULL;
NAB_strcpy(  &s, "t" );
if( i <= 7 ){
t = ( 40 + i );
trise = t * 3.380000E+00;
ttwist = t * 3.428000E+01;
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
NAB_matcpy( mat_x, newtransform(  - 1.732000E+01,  - 1.000000E+01, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00 ) );
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
else if( i > 7 && i <= 14 ){
t = ( 55 - i );
trise = t * 3.380000E+00;
ttwist = t * 3.428000E+01;
NAB_matcpy( mat_x, newtransform(  - 1.732000E+01, 1.000000E+01, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
m = wc_helix(  &n, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &s, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );
transformmol( mat_drz, m, NULL );
transformmol( mat_x, m, NULL );
mergestr( m1, "A", "last", m, "anti", "first" );
}
else{
t = ( i + 26 );
trise = t * 3.380000E+00;
ttwist = t * 3.428000E+01;
NAB_matcpy( mat_x, newtransform( 0.000000E+00, 2.000000E+01, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
if( i == 33 ){
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
NAB_strcpy(  &s, "t" );
printf( "%s %s %d %d %lf\n", s, n, i, t, ttwist );
}
putpdb( "ssd4D.pdb", m1, NULL );
mergestr( m_final, "Q", "last", m1, "A", "first" );
putpdb( "6h11.pdb", m_final, NULL );

freemolecule( m1 );
m1 = newmolecule(  );
addstrand( m1, "A" );
for( i = 1;i <= 33;i = i + 1 ){
n = NULL;
NAB_strcpy(  &s, "t" );
if( i <= 7 ){
t = ( 20 - i );
trise = t * 3.380000E+00;
ttwist = t * 3.428000E+01;
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
NAB_matcpy( mat_x, newtransform( 1.732000E+01, 1.000000E+01, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00 ) );
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
else if( i > 7 && i <= 14 ){
t = ( i + 5 );
trise = t * 3.380000E+00;
ttwist = t * 3.428000E+01;
NAB_matcpy( mat_x, newtransform( 0.000000E+00, 2.000000E+01, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
m = wc_helix(  &s, STEMP( __st0001__, "" ), STEMP( __st0002__, "dna" ),  &n, STEMP( __st0003__, "" ), STEMP( __st0004__, "dna" ), FTEMP( __ft0001__, 2.250000E+00 ), FTEMP( __ft0002__,  - 4.960000E+00 ), FTEMP( __ft0003__, 0.000000E+00 ), FTEMP( __ft0004__, 0.000000E+00 ), STEMP( __st0005__, "" ) );

transformmol( mat_drz, m, NULL );
transformmol( mat_x, m, NULL );
mergestr( m1, "A", "last", m, "sense", "first" );
}
else{
t = ( 34 - i );
trise = t * 3.380000E+00;
ttwist = t * 3.428000E+01;
NAB_matcpy( mat_x, newtransform(  - 1.732000E+01, 1.000000E+01, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00 ) );
NAB_matcpy( mat_drz, newtransform( 0.000000E+00, 0.000000E+00, trise, 0.000000E+00, 0.000000E+00, ttwist ) );
if( i == 33 ){
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
NAB_strcpy(  &s, "t" );
printf( "%s %s %d %d %lf\n", s, n, i, t, ttwist );
}
putpdb( "ssd3D.pdb", m1, NULL );
mergestr( m_final, "R", "last", m1, "A", "first" );
putpdb( "6h12.pdb", m_final, NULL );

putpdb( "dna_nanotube.pdb", m_final, NULL );



	exit( 0 );
}
