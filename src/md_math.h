/***************************************************************************
* Copyright (C) 2013-2014 Alexander V. Popov.
* 
* This file is part of Amateur Modeling of Biopolymers: Restoration, 
* Optimization, Solvation & Initial Analysis (AMBROSIA) source code.
* 
* AMBROSIA source code is free software; you can redistribute it and/or 
* modify it under the terms of the GNU General Public License as 
* published by the Free Software Foundation; either version 2 of 
* the License, or (at your option) any later version.
* 
* AMBROSIA source code is distributed in the hope that it will be 
* useful, but WITHOUT ANY WARRANTY; without even the implied 
* warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
* See the GNU General Public License for more details.
* 
* You should have received a copy of the GNU General Public License
* along with this program; if not, write to the Free Software 
* Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
***************************************************************************/
#ifndef MD_MATH_H
#define MD_MATH_H
/**
* @file
* @brief	Math functions.
*/

/** @brief	The Pi is always the Pi. */
#ifndef M_PI
#define M_PI				vec_t(3.14159265358979323846)
#endif

/** @brief	Epsilon value for convergence. */
#define ON_EPSILON			vec_t( 1e-5 )

/** @brief	Golden ratio constant. */
#define GOLDEN_RATIO		vec_t( 1.6180339887 )

/** @brief	Macro to convert degrees to radians. */
#define DEG2RAD( a )		( ( (a) * M_PI ) / 180 )
/** @brief	Macro to convert radians to degrees. */
#define RAD2DEG( a )		( ( (a) * 180 ) / M_PI )

/** @brief	Macro to set proper sign of a depending on b. */
#define CHECK_SIGN( a, b )	( ( (b) > 0 ) ? fabs( a ) : -fabs( a ) )

/* Macros for successive squaring */
/** @brief	Calculate a square of a number. */
#define SQR( a )			( ( a ) * ( a ) )
/** @brief	Calculate a cube of a number. */
#define CUBE( a )			( ( a ) * ( a ) * ( a ) )
/** @brief	Calculate a sixth power of a number. */
#define SIXTH( a )			SQR( ( a ) * ( a ) * ( a ) )
/** @brief	Calculate an eighth power of a number. */
#define EIGHTH( a )			SQR( ( a ) * ( a ) * ( a ) * ( a ) )
/** @brief	Calculate a twelfth power of a number. */
#define TWELFTH( a )		SIXTH( ( a ) * ( a ) )

/** @brief	Macro for clamping. */
#define CLAMP( f, a, b )	{ if ( ( f ) < ( a ) ) { ( f ) = ( a ); } else if ( ( f ) > ( b ) ) { ( f ) = ( b ); } }

/* Debug functions and macros for NAN checks */
/** @brief	Check a single-precision IEEE-754 floating-point number for a NaN value. */
inline bool checkNAN( float f ) {
	return (((*(int*)&f)&(0xff<<23))==(0xff<<23));
}
/** @brief	Check a double-precision IEEE-754 floating-point number for a NaN value. */
inline bool checkNAN( double f ) {
	return (((*((int*)&f+1))&(0x7ff<<20))==(0x7ff<<20));
}
#if defined( _DEBUG )
/** @brief	Debug macro for checking for a NaN value (empty in release). */
#define CHECK_NAN(x)	if ( checkNAN(x) ) __asm int 3
/** @brief	Debug macro for checking for a NaN vector (empty in release). */
#define CHECK_NAN3(x)	if ( checkNAN(x[0])||checkNAN(x[1])||checkNAN(x[2]) ) __asm int 3
#else
/** @brief	Debug macro for checking for a NaN value (empty in release). */
#define CHECK_NAN(x)
/** @brief	Debug macro for checking for a NaN vector (empty in release). */
#define CHECK_NAN3(x)
#endif

/** 
* @brief	Calculate difference between angles.
* @param	a1 : angle 1, in radians.
* @param	a2 : angle 2, in radians.
* @return	Difference between angles, in radians, modulo 2pi.
*/
inline vec_t AngleDiff( vec_t a1, vec_t a2 )
{
	vec_t delta = a1 - a2;
	if ( a1 > a2 ) {
		if ( delta >= M_PI ) delta -= 2*M_PI;
	} else {
		if ( delta <= -M_PI ) delta += 2*M_PI;
	}
	return delta;
}

/** 
* @brief	Clear a 3-component floating-point vector to zero.
* @param	a : source vector, will be cleared.
*/
inline void Vec3Clear( vec_t *a )
{
	a[0] = 0;
	a[1] = 0;
	a[2] = 0;
}


/** 
* @brief	Copy a 3-component floating-point vector.
* @param	a : source vector.
* @param	c : destination vector.
*/
inline void Vec3Copy( const vec_t *a, vec_t *c )
{
	c[0] = a[0];
	c[1] = a[1];
	c[2] = a[2];
}

/** 
* @brief	Add two 3-component floating-point vectors and store result.
* @details	Performs vector addition: C = A + B.
* @param	a : vector 1.
* @param	b : vector 2.
* @param	c : result.
*/
inline void Vec3Add( const vec_t *a, const vec_t *b, vec_t *c )
{
	c[0] = a[0] + b[0];
	c[1] = a[1] + b[1];
	c[2] = a[2] + b[2];
}
/** 
* @brief	Subtract two 3-component floating-point vectors and store result.
* @details	Performs vector addition: C = A - B.
* @param	a : menuend.
* @param	b : subtrahend.
* @param	c : result.
*/
inline void Vec3Sub( const vec_t *a, const vec_t *b, vec_t *c )
{
	c[0] = a[0] - b[0];
	c[1] = a[1] - b[1];
	c[2] = a[2] - b[2];
}
/** 
* @brief	Negate a 3-component floating-point.
* @details	Performs negation: C = -A.
* @param	a : source vector.
* @param	c : result.
*/
inline void Vec3Negate( const vec_t *a, vec_t *c )
{
	c[0] = -a[0];
	c[1] = -a[1];
	c[2] = -a[2];
}
/** 
* @brief	Scale a 3-component floating-point vector by a scalar value.
* @details	Performs scaling by scalar: C = A * S.
* @param	a : source vector.
* @param	s : scale scalar.
* @param	c : result.
*/
inline void Vec3Scale( const vec_t *a, const vec_t s, vec_t *c )
{
	c[0] = a[0] * s;
	c[1] = a[1] * s;
	c[2] = a[2] * s;
}
/** 
* @brief	Multication/addition of 3-component floating-point vectors.
* @details	Performs multiply-add operation: C = A + S * B.
* @param	a : vector 1.
* @param	s : scale scalar.
* @param	b : vector 2.
* @param	c : result.
*/
inline void Vec3MA( const vec_t *a, const vec_t s, const vec_t *b, vec_t *c )
{
	c[0] = a[0] + b[0] * s;
	c[1] = a[1] + b[1] * s;
	c[2] = a[2] + b[2] * s;
}
/** 
* @brief	Dot product of two 3-component floating-point vectors.
* @details	Calculates dot product: A . B.
* @param	a : vector 1.
* @param	b : vector 2.
* @return	Scalar dot product.
*/
inline vec_t Vec3Dot( const vec_t *a, const vec_t *b )
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
/** 
* @brief	Cross product of two 3-component floating-point vectors.
* @details	Calculates cross product: C = A x B.
* @param	a : vector 1.
* @param	b : vector 2.
* @param	c : result.
*/
inline void Vec3Cross( const vec_t *a, const vec_t *b, vec_t *c )
{
	c[0] = a[1] * b[2] - a[2] * b[1];
	c[1] = a[2] * b[0] - a[0] * b[2];
	c[2] = a[0] * b[1] - a[1] * b[0];
}
/** 
* @brief	Scalar triple product two 3-component floating-point vectors.
* @details	Calculates scalar triple product: A . (B x C).
* @param	a : vector 1.
* @param	b : vector 2.
* @param	c : vector 3.
* @return	Scalar triple product.
*/
inline vec_t Vec3Stp( const vec_t *a, const vec_t *b, const vec_t *c )
{
	return a[0] * ( b[1] * c[2] - b[2] * c[1] ) + a[1] * ( b[2] * c[0] - b[0] * c[2] ) + a[2] * ( b[0] * c[1] - b[1] * c[0] );
}
/** 
* @brief	Returns length of 3-component floating-point vector.
* @param	v : source vector.
* @return	Length of the source vector.
*/
inline vec_t Vec3Len( vec_t *v )
{
	return (vec_t)sqrt( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] );
}
/** 
* @brief	Returns square of length of 3-component floating-point vector.
* @param	v : source vector.
* @return	Square of length of the source vector.
*/
inline vec_t Vec3LenSq( const vec_t *v )
{
	return ( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] );
}
/** 
* @brief	Normalizes 3-component floating-point vector in place.
* @param	v : vector to normalize, will be overwritten with result.
* @return	Length of the original (non-normalized) vector.
*/
inline vec_t Vec3Nrm( vec_t *v )
{
	vec_t length, ilength;

	length = (vec_t)sqrt( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] );

	if ( length ) {
		ilength = vec_t( 1.0 ) / length;
		v[0] *= ilength;
		v[1] *= ilength;
		v[2] *= ilength;
	}
		
	return length;
}
/** 
* @brief	Normalizes 3-component floating-point vector in place.
* @param	v : vector to normalize, will be overwritten with result.
* @param	di : inverse of the length of the original (non-normalized) vector.
* @return	Length of the original (non-normalized) vector.
*/
inline vec_t Vec3Nrm2( vec_t *v, vec_t &di )
{
	vec_t length, ilength = 0;

	length = (vec_t)sqrt( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] );

	if ( length ) {
		ilength = vec_t( 1.0 ) / length;
		v[0] *= ilength;
		v[1] *= ilength;
		v[2] *= ilength;
	}
		
	di = ilength;
	return length;
}
/** 
* @brief	Normalizes 3-component floating-point vector in place.
* @details	This can a bit faster than #Vec3Nrm, since no checks are performed and nothing is returned.
* @param	v : vector to normalize, will be overwritten with result.
*/
inline void Vec3FastNrm( vec3_t v )
{
	vec_t ilength = vec_t( 1.0 ) / sqrt( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] );

	v[0] *= ilength;
	v[1] *= ilength;
	v[2] *= ilength;
}


/** 
* @brief	Line minimization of a multidimensional function.
* @details	This is used in conjugate gradients energy optimization algorithm.
* @param	dimension : number of arguments of a multidimensional function.
* @param	arguments : array of arguments (will be modified).
* @param	gradient : array of derivatives (will be modified).
* @param	F : interface to multidimensional function.
* @return	Value of a function at the minimum.
*/
extern vec_t LineMinimization( int dimension, vec4_t *arguments, vec4_t *gradient, vec_t (*F)(vec_t) );

/** 
* @brief	Replaces an n-by-n matrix, a, with the LU decomposition of a row-wise permutation of itself.
* @param	a : source matrix.
* @param	n : dimension of the matrix (matrix is square).
* @param	indx : the vector which records the row permutation effected by the partial pivoting.
* @param	vv : temporary matrix (dimension is the same as for a).
* @return	True if LU decomposition is successful, false if there is a singularity.
*/
extern bool ludcmp( vec_t *a, const int n, int *indx, vec_t *vv );
/** 
* @brief	Solves the set of n linear equations Ax = b.
* @details	LUBKSB must be used with the procedure LUDCMP to do this.
* @param	a : LU matrix.
* @param	n : dimension of the LU matrix (matrix is square).
* @param	indx : the vector which holds the row permutation effected by the partial pivoting (from LUDCMP).
* @param	b : input values of b; output results.
*/
extern void lubksb( const vec_t *a, const int n, const int *indx, vec_t *b );

/** 
* @brief	Sets random seed for the generator.
* @param	seed : seed value ( 0 = use current time).
*/
extern void SetRandomSeed( int seed );
/** 
* @brief	Generates a random float in a range [A,B].
* @param	flLow : minimum allowed value.
* @param	flHigh : maximum allowed value.
* @return	A random float value in a range.
*/
extern vec_t RandomFloat( vec_t flLow, vec_t flHigh );
/** 
* @brief	Generates a random integer in a range [A,B].
* @param	lLow : minimum allowed value.
* @param	lHigh : maximum allowed value.
* @return	A random integer value in a range.
*/
extern int RandomInt( int lLow, int lHigh );

#endif /*MD_MATH_H*/
