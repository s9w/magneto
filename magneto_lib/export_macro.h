#pragma once

#ifdef MAGNETOLIB_EXPORTS
#define CLASS_DECLSPEC    __declspec(dllexport)
#else
#define CLASS_DECLSPEC    __declspec(dllimport)
#endif
