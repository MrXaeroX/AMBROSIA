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
#ifndef MD_TYPES_H
#define MD_TYPES_H

/**
* @file
* @brief	Custom data type definitions.
* @details	The file defines custom data types.
*/

/** @brief	Unsigned 8-bit integer. */
#if defined(_MSC_VER)
typedef unsigned __int8 byte;
#else
typedef unsigned char byte;
#endif

/** @brief	Unsigned 16-bit integer. */
#if defined(_MSC_VER)
typedef unsigned __int16 word;
#else
typedef unsigned short word;
#endif

/** @brief	Unsigned 32-bit integer. */
#if defined(_MSC_VER)
typedef unsigned __int32 dword;
#else
typedef unsigned int dword;
#endif

/** @brief	Unsigned 64-bit integer. */
#if defined(_MSC_VER)
typedef unsigned __int64 qword;
#else
#if defined(__GNUC__)
#pragma GCC diagnostic ignored "-Wlong-long"
#endif
typedef unsigned long long qword;
#endif

/** @brief	Double-precision floating-point value.
* @details	May be changed to single (float), if such precision is a must, or you are RAM-limited (NOT recommended!).
*/
#if !defined(SINGLEPRECISION)
typedef double vec_t;
#else
typedef float vec_t;
#endif

/** @brief	2-component floating-point vector. */
typedef vec_t vec2_t[2];

/** @brief	3-component floating-point vector. */
typedef vec_t vec3_t[3];

/** @brief	4-component floating-point vector. */
typedef vec_t vec4_t[4];

#endif /*MD_TYPES_H*/
