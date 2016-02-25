#ifndef __CONFIG_H
#define __CONFIG_H

#ifdef  __cplusplus
# define BEGIN_C_DECLS  extern "C" {
# define END_C_DECLS    }
#else
# define BEGIN_C_DECLS
# define END_C_DECLS
#endif

#define OP_VERSION_MAJOR 1
#define OP_VERSION_MINOR 0
#define OP_VERSION_MICRO 0

#define OP_VERSION_STRING "1.0.0"

#define OP_VERSION_ENCODE(major, minor, micro) (((major)*10000)+((minor)*100)+(micro))

#define OP_VERSION OP_VERSION_ENCODE(OP_VERSION_MAJOR,OP_VERSION_MINOR,OP_VERSION_MICRO)

#define BUILD

#define XP_PC
#define _CRT_SECURE_NO_WARNINGS
#define _CRT_NONSTDC_NO_WARNINGS
#define HAVE_STD_SWAP

#if !defined(__GNUC__)
#define __attribute__(x)
#endif

#define OP_PRINTF(n,m)   __attribute__ ((__format__(__printf__,n,m)))

#if defined(WIN32) || defined(_MSC_VER) || defined(__CYGWIN__)
#ifdef BUILD_DLL
#define OP_EXPORT        __declspec(dllexport) __stdcall
#else
#define OP_EXPORT        __declspec(dllimport) __stdcall
#endif
#define OP_HIDDEN
#else
#if defined(GCC_VISIBILITY) && !defined(__sun)
#define OP_EXPORT        __attribute__ ((visibility("default")))
#define OP_HIDDEN        __attribute__ ((visibility("hidden")))
#elif defined(__SUNPRO_C) && (__SUNPRO_C >= 0x550)
#define OP_EXPORT
#define OP_HIDDEN        __hidden
#else
#define OP_EXPORT
#define OP_HIDDEN
#endif
#endif

#if defined(DEBUG)
#define NO_METAPROGS
#endif

#if defined(_MSC_VER)
#pragma warning(disable: 4227)
#endif

#endif // CONFIG_H
