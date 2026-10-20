/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  FirstTouch.h : an allocator that places large arrays across NUMA nodes
*
*  Linux places a page of memory on the NUMA node of the thread that first
*  writes it. A std::vector is zero-filled by the thread that resizes it, so
*  every large array of a run ends up on the node of the main thread, and
*  with threads on several nodes every parallel loop then reads across the
*  interconnect from that one node.
*
*  With numa on (util/Opt.h; the default), an allocation of at least kMinBytes is
*  taken from mmap (pages that nobody has written yet) and written once by
*  all threads in a static-schedule loop over its elements, the schedule
*  the parallel loops over the same arrays use: each page lands on the node
*  of the thread that will work on it. The values are written afterwards by
*  the container as before, so nothing a computation sees changes. With numa
*  off, or for small arrays, it is std::allocator.
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/
#ifndef YUZEN_FIRSTTOUCH_H
#define YUZEN_FIRSTTOUCH_H

#include <cstddef>
#include <cstring>
#include <memory>
#include <new>
#include "Opt.h"

#if defined(__unix__) || defined(__APPLE__)
#include <sys/mman.h>
#define YUZEN_HAVE_MMAP 1
#endif

namespace yuzen { namespace util {

const std::size_t kFirstTouchMinBytes = 256 * 1024;

template<class T>
struct FirstTouchAllocator
{
	typedef T value_type;

	FirstTouchAllocator() noexcept {}
	template<class U> FirstTouchAllocator(const FirstTouchAllocator<U>&) noexcept {}

	T* allocate(std::size_t n)
	{
		const std::size_t bytes = n * sizeof(T);
#ifdef YUZEN_HAVE_MMAP
		if(bytes >= kFirstTouchMinBytes && yuzen::opt::On(yuzen::opt::Numa))
		{
			void *p = mmap(0, bytes, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
			if(p == MAP_FAILED) throw std::bad_alloc();
			//  One write per element range, as a parallel loop over the
			//  elements would split it (schedule(static) over n)
			unsigned char *b = static_cast<unsigned char*>(p);
			const long long ne = (long long)n;
#if defined(_OPENMP)
#pragma omp parallel for schedule(static)
#endif
			for(long long i = 0; i < ne; ++i)
				std::memset(b + (std::size_t)i * sizeof(T), 0, sizeof(T));
			return static_cast<T*>(p);
		}
#endif
		return std::allocator<T>().allocate(n);
	}

	void deallocate(T *p, std::size_t n) noexcept
	{
		const std::size_t bytes = n * sizeof(T);
#ifdef YUZEN_HAVE_MMAP
		if(bytes >= kFirstTouchMinBytes && yuzen::opt::On(yuzen::opt::Numa))
		{
			munmap(p, bytes);
			return;
		}
#endif
		std::allocator<T>().deallocate(p, n);
	}

	template<class U> bool operator==(const FirstTouchAllocator<U>&) const noexcept { return true; }
	template<class U> bool operator!=(const FirstTouchAllocator<U>&) const noexcept { return false; }
};

} } // namespace yuzen::util

#endif // YUZEN_FIRSTTOUCH_H
