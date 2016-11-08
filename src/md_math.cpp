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
#include "md_main.h"
#include "md_math.h"

/**
* Adapted function mnbrak from Numeric Recipes.
*/
static void BracketMinimization( vec_t *ax, vec_t *bx, vec_t *cx, vec_t *fa, vec_t *fb, vec_t *fc, vec_t (*F)(vec_t) )
{
	const vec_t tiny = vec_t( 1.0e-20 );
	const vec_t glimit = 100;
	vec_t fu, temp;
	
	*fa = F( *ax );
	*fb = F( *bx );

	if ( *fb > *fa ) {
		temp = *ax; *ax = *bx; *bx = temp;
		temp = *fb; *fb = *fa; *fa = temp;
	}

	*cx = (*bx) + GOLDEN_RATIO * ( *bx - *ax );
	*fc = F( *cx );

	while ( *fb > *fc ) {
		vec_t r = ( *bx - *ax ) * ( *fb - *fc );
		vec_t q = ( *bx - *cx ) * ( *fb - *fa );
		vec_t u = ( *bx ) - ( ( *bx - *cx ) * q - ( *bx - *ax ) * r ) / ( 2.0f * CHECK_SIGN( std::max( fabs( q - r ), tiny ), q - r ) );
		vec_t ulim = ( *bx ) + glimit * ( *cx - *bx );

		if ( ( *bx - u ) * ( u - *cx ) > 0 ) {
			fu = F( u );
			if ( fu < *fc ) {
				*ax = *bx;
				*bx = u;
				*fa = *fb;
				*fb = fu;
				return;
			} else if ( fu > *fb ) {
				*cx = u;
				*fc = fu;
				return;
			}
			u = ( *cx ) + GOLDEN_RATIO * ( *cx - *bx );
			fu = F( u );
		} else if ( ( *cx - u ) * ( u - ulim ) > 0 ) {
			fu = F( u );
			if ( fu < *fc ) {
				*bx = *cx; *cx = u; u = *cx + GOLDEN_RATIO * ( *cx - *bx );
				*fb = *fc; *fc = fu; fu = F( u );
			}
		} else if ( ( u - ulim ) * ( ulim - *cx ) >= 0 ) {
			u = ulim;
			fu = F( u );
		} else {
			u = ( *cx ) + GOLDEN_RATIO * ( *cx - *bx );
			fu = F( u );
		}

		*ax = *bx; *bx = *cx; *cx = u;
		*fa = *fb; *fb = *fc; *fc = fu;
	}
}

/**
* Adapted function brent from Numeric Recipes.
*/
static vec_t BrentMinimization( vec_t ax, vec_t bx, vec_t cx, vec_t tolerance, vec_t (*F)(vec_t), vec_t *xmin )
{
	const int maxIterations = 100;
	const vec_t zeps = vec_t( 1.0e-10 );
	const vec_t cGold = vec_t( 0.3819660 );

	vec_t x, w, v, u, p, q, r;
	vec_t fx, fw, fv, fu;
	vec_t etemp, e = 0;
	vec_t d = 0;

	vec_t a = std::min( ax, cx );
	vec_t b = std::max( ax, cx );

	x = w = v = bx;
	fx = fw = fv = F( x );

	for ( int iteration = 0; iteration < maxIterations; ++iteration ) {
		vec_t xm = vec_t(0.5) * ( a + b );
		vec_t tol1 = tolerance * fabs(x) + zeps;
		vec_t tol2 = 2 * tol1;

		if ( fabs( x - xm ) <= ( tol2 - 0.5f * ( b - a ) ) ) {
			*xmin = x;
			return fx;
		}

		if ( fabs( e ) > tol1 ) {
			r = ( x - w ) * ( fx - fv );
			q = ( x - v ) * ( fx - fw );
			p = ( x - v ) * q - ( x - w ) * r;
			q = 2.0f * ( q - r );
			if ( q > 0 ) p = -p;
			q = fabs( q );
			etemp = e;
			e = d;
			if ( fabs( p ) >= fabs( 0.5 * q * etemp ) || p <= q * ( a - x ) || p >= q * ( b - x ) ) {
				d = cGold * ( e = ( x >= xm ? a - x : b - x ) );
			} else {
				d = p / q;
				u = x + d;
				if ( u - a < tol2 || b - u < tol2 )
					d = CHECK_SIGN( tol1, xm - x );
			}
		} else {
			d = cGold * ( e = ( x >= xm ? a - x : b - x ) );
		}
		u = ( fabs( d ) >= tol1 ? ( x + d ) : ( x + CHECK_SIGN( tol1, d ) ) );
		fu = F( u );
		if ( fu <= fx ) {
			if ( u >= x ) a = x; else b = x;
			v = w; w = x; x = u;
			fv = fw; fw = fx; fx = fu;
		} else {
			if ( u < x ) a = u; else b = u;
			if ( fu <= fw || w == x ) {
				v = w; w = u; fv = fw; fw = fu;
			} else if ( fu <= fv || v == x || v == w ) {
				v = u; fv = fu;
			}
		}
	}

	Log().Fatal( "BrentMinimization: too many iterations (limit = %i)\n", maxIterations );
	*xmin = x;
	return fx;
}

/**
* Adapted function linmin from Numeric Recipes.
*/
vec_t LineMinimization( int dimension, vec4_t *arguments, vec4_t *gradient, vec_t (*F)(vec_t) )
{
	const vec_t tolerance = vec_t( 0.01 );
	vec_t fa, fb, fc, xmin, fret;
	vec_t lambda_init = vec_t( 0.1 );
	vec_t gradSquare = 0;

	for ( int i = 0; i < dimension; ++i )
		gradSquare += Vec3LenSq( gradient[i] );

	if ( gradSquare > 1 ) 
		lambda_init = lambda_init / sqrt( gradSquare );

	vec_t ax = 0.0f;
	vec_t bx = gradSquare;
	vec_t cx = 0.0f;

	BracketMinimization( &ax, &bx, &cx, &fa, &fb, &fc, F );
	fret = BrentMinimization( ax, bx, cx, tolerance, F, &xmin );

	for ( int i = 0; i < dimension; ++i ) {
		Vec3Scale( gradient[i], xmin, gradient[i] );
		Vec3Add( arguments[i], gradient[i], arguments[i] );
	}

	return fret;
}

/**
* LUDCMP
* Replaces an n-by-n matrix, a, with the LU decomposition of a row-wise permutation of itself
*/
bool ludcmp( vec_t *a, const int n, int *indx, vec_t *vv )
{
	int imax = 0;
	bool bsng = false;

    for ( int i = 0; i < n; i++) {
		vec_t big = 0;
		vec_t temp;
		for ( int j = 0; j < n; j++ ) {
			if ((temp = fabs(a[i*n+j])) > big) big = temp;
		}
        if (big == 0) bsng = true;
        vv[i] = vec_t( 1.0 ) / big;
    }

    for ( int j = 0; j < n; j++ ) {
		for ( int i = 0; i < j; i++ ) {
			vec_t sum = a[i*n+j];
			for ( int k = 0; k < i; k++ ) sum -= a[i*n+k]*a[k*n+j];
			a[i*n+j] = sum;
        }
        vec_t big = 0;
		vec_t temp;
        for ( int i = j; i < n; i++ ) {
			vec_t sum = a[i*n+j];
            for ( int k = 0; k < j; k++ ) sum -= a[i*n+k]*a[k*n+j];
            a[i*n+j] = sum;
            if ( (temp = vv[i]*fabs(sum)) >= big) {
				big = temp;
				imax = i;
            }
        }
        if (j != imax) {
			for ( int k = 0; k < n; k++ ) {
				temp = a[imax*n+k];
                a[imax*n+k] = a[j*n+k];
                a[j*n+k] = temp;
            }
            vv[imax] = vv[j];
        }
        indx[j] = imax;
        if (a[j*n+j] == 0) 
			a[j*n+j] = vec_t(1.0e-20);

        if (j != n) {
			temp = vec_t(1.0)/(a[j*n+j]);
            for ( int i = j+1; i < n; i++) a[i*n+j] *= temp;
        }
    }
	return bsng;
}
/**
* LUBKSB
* Solves the set of n linear equations Ax = b (LUBKSB must be used with the procedure LUDCMP to do this)
*/
void lubksb( const vec_t *a, const int n, const int *indx, vec_t *b )
{
	int ii = -1;
    for ( int i = 0; i < n; i++ ) {
		int ip = indx[i];
        vec_t sum = b[ip];
        b[ip] = b[i];
        if (ii >= 0) {
			for ( int j = ii; j <= i-1; j++) 
				sum -= a[i*n+j]*b[j];
		} else if (sum) {
			ii = i;
		}
        b[i] = sum;
    }

    for ( int i = n-1; i >= 0; i-- ) {
		vec_t sum = b[i];
        for ( int j = i+1; j < n; j++ ) 
			sum -= a[i*n+j]*b[j];
        b[i] = sum / a[i*n+i];
    }
}

/*
* Random number generation
*/
#define MAX_RANDOM_RANGE	0x7FFFFFFFUL
#define MAX_RANDOM_RANGE	0x7FFFFFFFUL
#define IA					16807
#define IM					2147483647
#define IQ					127773
#define IR					2836
#define NTAB				32
#define NDIV				(1+(IM-1)/NTAB)
#define AM					(1.0/IM)
#define EPS					1.2e-7
#define RNMX				(1.0 - EPS)

/**
* SetRandomSeed
* Sets random seed for the generator.
*/
static int idum = 0;

void SetRandomSeed( int seed )
{
	if ( seed ) idum = seed;
	else idum = -(int)time( NULL );

	if ( 1000 < idum ) idum = -idum;
	else if ( -1000 < idum ) idum -= 22261048;
}

/**
* IRAN1
* Generates random integer value with uniform distribution.
*/
static int iran1( void )
{
	int	j, k;

	static int iy = 0;
	static int iv[NTAB];
	
	if ( idum <= 0 || !iy ) {
		if ( -(idum) < 1 ) idum = 1;
		else idum = -(idum);
		for ( j = NTAB + 7; j >= 0; j-- ) {
			k = (idum) / IQ;
			idum = IA * (idum - k * IQ) - IR * k;
			if ( idum < 0 ) idum += IM;
			if ( j < NTAB ) iv[j] = idum;
		}
		iy = iv[0];
	}
	
	k = (idum)/IQ;
	idum = IA * (idum - k * IQ) - IR * k;
	if ( idum < 0 ) idum += IM;
	j = iy / NDIV;
	iy = iv[j];
	iv[j] = idum;
	return iy;
}

/**
* IRAN1
* Generates random float value with uniform distribution.
*/
static vec_t fran1( void )
{
	vec_t temp = (vec_t)AM * iran1();
	if ( temp > RNMX ) return (vec_t)RNMX;
	else return temp;
}

vec_t RandomFloat( vec_t flLow, vec_t flHigh )
{
	if ( idum == 0 ) SetRandomSeed( 0 );
	return (fran1() * (flHigh - flLow)) + flLow;
}

int RandomInt( int lLow, int lHigh )
{
	unsigned int maxAcceptable;
	unsigned int n, x = lHigh-lLow + 1; 	

	if ( idum == 0 ) SetRandomSeed( 0 );

	if ( x <= 0 || MAX_RANDOM_RANGE < x-1 )
		return lLow;

	maxAcceptable = MAX_RANDOM_RANGE - ( (MAX_RANDOM_RANGE+1) % x );
	do {
		n = iran1();
	} while ( n > maxAcceptable );
	
	return lLow + (n % x);
}
