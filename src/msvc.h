#pragma once

#define WIN32_LEAN_AND_MEAN

/*
        Disable deprecated CRT functions limitation.
        For further information visit:
        http://msdn.microsoft.com/en-us/library/ms235384.aspx
*/
#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_NON_CONFORMING_SWPRINTFS

#ifndef _DEBUG
        /*
                Disable checked iterators.
                For further information visit:
                http://msdn.microsoft.com/en-us/library/aa985965.aspx
        */
        #define _SECURE_SCL 0
        #define _SECURE_SCL_THROWS 0
        #define _SCL_SECURE_NO_DEPRECATE
        #define _HAS_ITERATOR_DEBUGGING 0
#endif

/*
        Suppress min/max conflicts with STL.
        For further information visit:
        http://support.microsoft.com/kb/143208
*/
#ifndef NOMINMAX
        #define NOMINMAX
#endif

/*
        If the compiler supports C++0x, disable the workarounds
        and call the standard headers.
*/
/*
        Debug help only safety works with Visual C++
*/
#ifdef __USE_MINIDUMP__
        #ifndef __EXCEPTION_TRACER__
                #define __EXCEPTION_TRACER__
        #endif
#endif
