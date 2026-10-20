// omp_overhead.cpp -- what an OpenMP parallel loop costs on this machine
//
//     g++ -O2 -fopenmp tools/scaling/omp_overhead.cpp -o omp_overhead
//     for t in 1 2 4 8 16 28 56 112; do
//       OMP_NUM_THREADS=$t OMP_PROC_BIND=close OMP_PLACES=cores ./omp_overhead
//     done
//
// For each loop length n it times, averaged over many repetitions:
//   empty   a parallel for whose body does nothing (fork, schedule, join)
//   axpy    y[i] += a x[i]: a loop of the length Yuzen's grid sweeps have
//   dot     a sum over the loop, with a reduction clause
// and prints one line per n: threads, n, and the three in microseconds. The
// serial time of the same axpy is printed too, so that the length below
// which a parallel loop is slower than a serial one can be read off. The
// grid sweeps of a 1-D film of 300 cells and a 60 x 60 film are n = 300 and
// 3600; a step runs some tens of such loops per Newton pass.
//
// Copyright 2017-2026 Suguru Shiratori, Tokyo City University
// SPDX-License-Identifier: Apache-2.0
#include <omp.h>
#include <chrono>
#include <cstdio>
#include <vector>

static double now() {
	return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

int main()
{
	const int ns[] = {300, 3600, 16384, 65536, 262144, 1048576};
	const int T = omp_get_max_threads();
	std::printf("# threads\tn\tempty_us\taxpy_us\tdot_us\tserial_axpy_us\n");
	for(int n : ns)
	{
		std::vector<double> x((size_t)n, 1.0), y((size_t)n, 0.5);
		#pragma omp parallel for schedule(static)
		for(int i = 0; i < n; i++) { x[(size_t)i] = 1.0 + 1e-9 * i; y[(size_t)i] = 0.5; }	// first touch
		const int reps = n <= 4000 ? 20000 : (n <= 70000 ? 2000 : 200);
		volatile double sink = 0;
		double t0, e, ax, dt, se;
		// warm up
		for(int r = 0; r < 100; r++) {
			#pragma omp parallel for schedule(static)
			for(int i = 0; i < n; i++) y[(size_t)i] += 1e-12 * x[(size_t)i];
		}
		t0 = now();
		for(int r = 0; r < reps; r++) {
			#pragma omp parallel for schedule(static)
			for(int i = 0; i < n; i++) { (void)i; }
		}
		e = (now() - t0) / reps;
		t0 = now();
		for(int r = 0; r < reps; r++) {
			#pragma omp parallel for schedule(static)
			for(int i = 0; i < n; i++) y[(size_t)i] += 1e-12 * x[(size_t)i];
		}
		ax = (now() - t0) / reps;
		t0 = now();
		for(int r = 0; r < reps; r++) {
			double s = 0;
			#pragma omp parallel for schedule(static) reduction(+:s)
			for(int i = 0; i < n; i++) s += x[(size_t)i] * y[(size_t)i];
			sink = sink + s;
		}
		dt = (now() - t0) / reps;
		t0 = now();
		for(int r = 0; r < reps; r++)
			for(int i = 0; i < n; i++) y[(size_t)i] += 1e-12 * x[(size_t)i];
		se = (now() - t0) / reps;
		std::printf("%d\t%d\t%.3f\t%.3f\t%.3f\t%.3f\n", T, n, e * 1e6, ax * 1e6, dt * 1e6, se * 1e6);
		(void)sink;
	}
	return 0;
}
