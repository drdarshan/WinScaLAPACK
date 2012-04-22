#ifndef _IMPORT_H_
#define _IMPORT_H_

#if DYNAMIC
#define DLLIMPORT __declspec(dllimport)
#else
#define DLLIMPORT
#endif
#endif