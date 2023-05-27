// stdafx.h : include file for standard system include files,
//  or project specific include files that are used frequently, but
//      are changed infrequently
//

#if !defined(AFX_STDAFX_H__43BFE2DD_E831_47E1_B38A_23A4ABAA0610__INCLUDED_)
#define AFX_STDAFX_H__43BFE2DD_E831_47E1_B38A_23A4ABAA0610__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#define VC_EXTRALEAN		// Exclude rarely-used stuff from Windows headers

#include <afxwin.h>         // MFC core and standard components
//#include <afxext.h>         // MFC extensions
//#include <afxdisp.h>        // MFC Automation classes
//#include <afxdtctl.h>		// MFC support for Internet Explorer 4 Common Controls
//#ifndef _AFX_NO_AFXCMN_SUPPORT
//#include <afxcmn.h>			// MFC support for Windows Common Controls
//#endif // _AFX_NO_AFXCMN_SUPPORT
//#include "./h/tchar.h"
#include "Point2_T.h"
#include "shxUtils.h"
// TODO: reference additional headers your program requires here


//  MBCS and Unicode Translation Flags.
//  Please use Unicode, either UTF-16 (WCHAR) or UTF-8 (CP_UTF8)
//
// MB_PRECOMPOSED and MB_COMPOSITE are deprecated, not recommended, and
// provide out-of-date behavior.
// Windows typically uses Unicode Normalization Form C type sequences,
// If explicit normalization forms are required, please use NormalizeString.
#define MB_PRECOMPOSED            0x00000001  // DEPRECATED: use single precomposed characters when possible.
#define MB_COMPOSITE              0x00000002  // DEPRECATED: use multiple discrete characters when possible.
#define MB_USEGLYPHCHARS          0x00000004  // DEPRECATED: use glyph chars, not ctrl chars
#define MB_ERR_INVALID_CHARS      0x00000008  // error for invalid chars

// WC_COMPOSITECHECK, WC_DISCARDNS and WC_SEPCHARS are deprecated, not recommended,
// and provide out-of-date behavior.
// Windows typically uses Unicode Normalization Form C type sequences,
// If explicit normalization forms are required, please use NormalizeString.
#define WC_COMPOSITECHECK         0x00000200  // convert composite to precomposed
#define WC_DISCARDNS              0x00000010  // discard non-spacing chars          // Used with WC_COMPOSITECHECK
#define WC_SEPCHARS               0x00000020  // generate separate chars            // Used with WC_COMPOSITECHECK
#define WC_DEFAULTCHAR            0x00000040  // replace w/ default char            // Used with WC_COMPOSITECHECK
#if (WINVER >= 0x0600)
#define WC_ERR_INVALID_CHARS      0x00000080  // error for invalid chars
#endif

#if(WINVER >= 0x0500)
#define WC_NO_BEST_FIT_CHARS      0x00000400  // do not use best fit chars
#endif /* WINVER >= 0x0500 */

#define GDI_ERROR (0xFFFFFFFFL)

const double	IC_PI			= 3.14159265358979323846;
const double	IC_TWOPI		= 2 * IC_PI;
#define IC_ZRO           1.0E-10

#define _TEOF       EOF
#define __T(x)      x

#define _T(x)       __T(x)
#define _TEXT(x)    __T(x)


#ifndef _TCHAR_DEFINED
typedef char TCHAR, * PTCHAR;
typedef unsigned char TBYTE, * PTBYTE;
#define _TCHAR_DEFINED
#endif /* !_TCHAR_DEFINED */

#ifndef __TCHAR_DEFINED
typedef char            _TCHAR;
typedef signed char     _TSCHAR;
typedef unsigned char   _TUCHAR;
typedef unsigned char   _TXCHAR;
typedef unsigned int    _TINT;
#define __TCHAR_DEFINED
#endif  /* __TCHAR_DEFINED */


typedef unsigned long       DWORD;
typedef int                 BOOL;
typedef unsigned char       BYTE;
typedef unsigned short      WORD;
typedef float               FLOAT;

typedef int                 INT;
typedef unsigned int        UINT;
typedef unsigned int* PUINT;

//typedef struct POINT
//{
//	long  x;
//	long  y;
//} POINT;
//
//typedef struct GLYPHMETRICS {
//	UINT    gmBlackBoxX;
//	UINT    gmBlackBoxY;
//	POINT   gmptGlyphOrigin;
//	short   gmCellIncX;
//	short   gmCellIncY;
//} GLYPHMETRICS;
//
//#define DECLARE_HANDLE(name) struct name##__{int unused;}; typedef struct name##__ *name
//DECLARE_HANDLE(HDC);

#define strnsame(a,b,c) (!strncmp((a),(b),(c)))

#ifndef FABS
//on XP we observed that the fabs() function
//performs 10+ times slower than the macro version underhere.
//For some drawings 15+ % of total processing time was spent in fabs()!!!
//We did not yet profile fabs() on older operating systems
#define FABS(a)  ((a) >= 0.0 ? (a) : -(a))
//#define FABS(a)  fabs(a)
//#define FABS(a)  fastAbs(a)

#endif
const double IcadFuzzyDistanceX = 1.0e-11;
inline bool icadPointEqual(const double p1[3], const double p2[3], double tol= IcadFuzzyDistanceX)
{
	if (FABS(p2[0])<=tol) {
		if (FABS(p1[0])>tol)
			return false;
	}
	else if (FABS((p1[0]/p2[0])-1.0) >= tol)
		return false;

	if (FABS(p2[1])<=tol) {
		if (FABS(p1[1])>tol)
			return false;
	}
	else if (FABS((p1[1]/p2[1])-1.0) >= tol)
		return false;

	if (FABS(p2[2])<=tol) {
		if (FABS(p1[2])>tol)
			return false;
	}
	else if (FABS((p1[2]/p2[2])-1.0) >= tol)
		return false;

	return true;
}

// If you don't know the scale of the underlying measurements,
// using the test "abs(a/b - 1) < epsilon"
// is likely to be more robust than simply comparing the difference

inline bool icadRealEqual(double r1, double r2, double tol = IcadFuzzyDistanceX)
{
	if (FABS(r2)<=tol) {
		if (FABS(r1)>tol)
			return false;
	}
	else {
		if (FABS((r1/r2)-1.0) > tol)
			return false;
	}
	return true;
}
//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_STDAFX_H__43BFE2DD_E831_47E1_B38A_23A4ABAA0610__INCLUDED_)
