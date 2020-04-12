#ifndef HTOA_MEMORY_H
#define HTOA_MEMORY_H

#include <new>
#include <cstdlib>
#include <stdlib.h>

#if !defined(__APPLE__)
#include <malloc.h>
#endif

#define CACHE_LINE_ALIGNMENT 64
#define SIMD_ALIGNMENT 16;
#ifdef __GNUC__
    #define AS_ALIGN(x)  __attribute__ ((aligned(x)))
#else
    #define AS_ALIGN(x) __declspec(align(x))
#endif

inline void* aligned_malloc(size_t size, size_t alignment)
{
#ifdef __MINGW_GCC_VERSION
    return _aligned_malloc(size, alignment);
#else
    #ifdef __GNUC__
#if _WIN32
    void *pResult = ::memalign(alignment, size);
    if (pResult == NULL)
#else
    void *pResult = NULL;
    if (::posix_memalign(&pResult, alignment, size) != 0)
#endif
    {
        pResult = ::malloc(size);
    }
    return pResult;
#else
    return _aligned_malloc(size, alignment);
#endif
#endif
}

inline void aligned_free(void *pPtr)
{
     if (pPtr == NULL)
        return;

#ifdef __GNUC__
    free(pPtr);
#else
    _aligned_free(pPtr);
#endif
}

template<typename T, size_t kAlignment>
inline T *alignedAlloc(size_t numInstances)
{
    return reinterpret_cast<T*>(aligned_malloc(sizeof(T) * numInstances, kAlignment));
}

template<typename T>
inline void alignedFree(T* data)
{
    aligned_free(data);
}

template<size_t kAlignment>
class AlignBase
{
public:
    void* operator new(size_t size)
    {
        return aligned_malloc(size, kAlignment);
    }

    void* operator new[](size_t size)
    {
        return aligned_malloc(size, kAlignment);
    }

    void operator delete(void *pPtr)
    {
        aligned_free(pPtr);
    }

    void operator delete[](void *pPtr)
    {
        aligned_free(pPtr);
    }
};

#endif // HTOA_MEMORY_H
