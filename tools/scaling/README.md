# SMP scaling measurement kit

Yuzen uses **shared-memory parallelism (OpenMP) only**; it does not support MPI.
This directory holds everything needed to measure, reproducibly, "how much
faster does it get on 1 machine as the thread count is raised".

```
make_cases.py     generates the measurement cases and manifest.tsv
scaling.sh        sweeps thread counts, backends and cases, appending to a TSV
summarize.py      turns the TSV into a table of speed-up and parallel efficiency (text / latex / csv)
results/          where measured results are kept (scaling-<host>.tsv)
```

---

## Procedure

### 1. Build

Build with the AMGCL backend, so that both backends can be measured: which of
the two is faster depends on the dimension, the size and the thread count
(the manual, "Which backend to ask for").

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DYUZEN_WITH_AMGCL=ON
cmake --build build -j
./build/bin/yuzen --version            # check the version
```

If the AMGCL headers are absent, CMake fetches them with git. If that
fails, the build degrades to Eigen only, and even if `amgcl` is given to
`--backends` the run actually uses `eigen`. In that case the TSV's `status`
column carries `used:eigen`, so the results cannot be tabulated unnoticed.

### 2. Prepare the cases

Generated copies are committed to the repository, so normally nothing needs doing.

```sh
python3 tools/scaling/make_cases.py --check   # do the local files match the generated ones?
python3 tools/scaling/make_cases.py           # regenerate them
```

### 3. Start with a small run to check

If the checkout has not kept the executable bit — a clone on Windows, an
unzipped archive, a copy through a filesystem without permissions — the script
gives `Permission denied`. Prefix it with `bash` and it runs either way.

```sh
tools/scaling/scaling.sh --tier quick --repeats 1 --threads "1 4"
python3 tools/scaling/summarize.py tools/scaling/results/scaling-*.tsv
```

Check that the `status` column is `ok` throughout and that the `steps` column
matches `steps` in manifest.tsv. If this breaks down, the later numbers measure
not "parallel performance" but "failure to converge".

### 4. The real run

A sweep aimed at 112 cores (2 sockets x 56 cores).

```sh
tools/scaling/scaling.sh \
    --threads  "1 2 4 8 14 28 56 84 112" \
    --backends "eigen amgcl" \
    --repeats  3
```

**Do not specify a thread count above the number of physical cores.** Even
where Hyper-Threading is enabled, keep `--threads` within the physical core
count. Oversubscription is not a small loss but an order-of-magnitude one
(there is a case where going to 4 threads on a 2-core machine made it 10 times
slower).

If it stops part-way through, **simply run the same command again**. Any
`(case, backend, threads, repeat)` already in the TSV is skipped and
measurement resumes from where it left off.

For a guide to the time required, take "seconds per step on one thread x number
of steps" from the case table (the comment at the top of `make_cases.py`) and
sum it over the thread counts. `--dry-run` shows what would be run without
running it.

### 5. Summarize

```sh
python3 tools/scaling/summarize.py tools/scaling/results/scaling-*.tsv
python3 tools/scaling/summarize.py tools/scaling/results/*.tsv --format latex
```

When passing on results, just send **the TSV in `results/` as it is**. The
machine's particulars (CPU name, socket count, core count, NUMA node count,
memory size, OS, Yuzen version, executable path, `OMP_PROC_BIND`) are written
at the top of the file prefixed with `#`, so no separate explanation is needed.

---

## How the measurement is designed

The premise of the measurement is that **the amount of work does not change
when the thread count changes**. There are three things that would break that, and
all of them are suppressed on the case side.

1. **Automatic time-step control** — `[Numerics.TimeStepControl] Enable = false`.
   With automatic control in place, the history of step sizes can diverge
   through rounding differences.

2. **Termination by time** — `Numerics.MaxSteps` fixes the number of steps, and
   `EndTime` is placed out of reach (1e30). Every run therefore takes the same
   number of steps.

3. **Field output** — `NeutralFile = false`, `VTKFile = false`, with all output
   intervals at 1e30. The output interval **also acts as an upper bound on the
   time step**, so making it small would change the step size as well. It also
   keeps I/O time out of the measurement.

The thread count is also the one asked for: `scaling.sh` sets
`YUZEN_THREAD_CAP=0`, which turns off the solver's cap on the threads of a
small grid (at most max(8, cells/1500) by default; see the manual,
sec:parallel:threads). With the cap on, every count above the cap would
measure the same run.

What is measured is the wall-clock time of **the time-step loop only**, taken
from `Time loop [s]` in JobInfo.txt. `Elapsed Time [s]` also covers grid
construction, the initial field and the initial output, so as the grid grows
the set-up stops being negligible and it cannot serve as a parallel-scaling
metric.

### How the backend is specified

None of the shipped cases writes a `Backend` key (the default is `auto`).
When `scaling.sh` copies the input into the working directory, it inserts

```toml
[Numerics.HeightEquation.LinearSolver]
Backend = "amgcl"
...
[Numerics.ConcEquation.LinearSolver]
Backend = "amgcl"
...
```

into **both** blocks. With only one of them set, the other stays on `auto` and
it becomes unclear what is being measured.

### How the grid size is chosen

Filling 112 threads takes a correspondingly large problem. The guide is the
block decomposition of the `Vector` reduction,

```
number of blocks = min(256, number of elements / 512)
```

so that

* for blocks to reach all 112 threads, **element count >= 57,344**
* for the load to approach uniformity, **element count >= 131,072** (the block count is capped at 256)

is required. The 2-D 512² = 262,144 and the 1-D Nx = 131,072 are at that
level. Put the other way round, **a small case cannot in principle be made
faster on 112 cores**. The smaller of the 1-D cases is included precisely to
show that "it does not get faster".

Note that because the reduction's block count is capped at 256, **the reduction
is parallelized no further beyond 256 threads**. This is the design price paid
for bit reproducibility (the same values coming out whatever the thread count),
and it is the ceiling on the performance side.

### List of cases

`cases/scaling/manifest.tsv` is the list, and `scaling.sh` reads it too.

| Column | Meaning |
|---|---|
| `case` | Case name. `cases/scaling/<case>.toml` |
| `dim` | 1 or 2 |
| `nx`, `ny` | Number of grid points |
| `nodes` | `nx*ny` (number of unknowns) |
| `dt` | Fixed time step [s] |
| `steps` | Number of steps taken |
| `tier` | `quick` (for checking operation) / `full` (for the real run) |

`dt` has been chosen **from measurements**, so that "Newton settles in 2-3
iterations at any size". Get this wrong and Newton stalls, so that what is
measured is not parallel performance but poor convergence. To change the number
of steps, edit the table at the top of `make_cases.py` and regenerate
(`scaling.sh` does not rewrite the inputs).

---

## Format of the result file

`results/scaling-<host>.tsv`. UTF-8, tab separated, one run per line. It is
written by appending, so the same file can be added to any number of times.

### Header (lines beginning with `#`)

A sequence of `key<TAB>value`. Written just once, when the file is created.

```
# Yuzen SMP scaling measurement
# format-version	1
# generated-by	tools/scaling/scaling.sh
# date	2026-09-17T00:47:05Z
# host	hpc01
# cpu	Intel(R) Xeon(R) Platinum 8480+
# sockets	2
# cores-per-socket	56
# threads-per-core	2
# numa-nodes	2
# memory-kb	527654912
# os	Linux 5.14.0
# yuzen-version	Yuzen-1.0.0
# executable	/home/you/yuzen/build/bin/yuzen
# omp-proc-bind	close
```

### Body

| Column | Type | Meaning |
|---|---|---|
| `case` | string | Case name (`case` in manifest.tsv) |
| `dim` | 1 \| 2 | Dimension |
| `nx`, `ny` | integer | Number of grid points |
| `nodes` | integer | Number of unknowns |
| `backend` | `auto`\|`eigen`\|`amgcl` | The linear solver **that was requested** |
| `threads` | integer | `OMP_NUM_THREADS` |
| `repeat` | integer | Which repetition of the same condition (1-based) |
| `steps` | integer \| `NA` | Number of steps actually taken |
| `wall_s` | real \| `NA` | Wall-clock time of **the time-step loop only** [s] |
| `nk_total` | integer \| `NA` | Total number of Newton iterations (0 for ConcZ) |
| `bcg_total` | integer \| `NA` | Total number of linear iterations |
| `status` | string | See the table below |

Values of `status`:

| Value | Meaning |
|---|---|
| `ok` | Finished normally, no warnings, step count as specified |
| `warned` | Finished normally but Warnings.txt is not empty |
| `fatal` | FatalErrors.txt is not empty |
| `exitN` | Died with exit code N (timeout is 124) |
| `nodata` | JobInfo.txt could not be read |
| `stepsN` | Stopped at N steps (differs from what was specified) |
| `.../used:X` | X was used instead of the backend that was asked for |

`summarize.py` discards anything other than `ok` by default. `--include-bad`
brings them back in.

### Conventions of the summary

* Speed-up `S(p) = T(1)/T(p)`, parallel efficiency `E(p) = S(p)/p`.
* Where there are repetitions, the **minimum** is taken. The minimum rather
  than the mean, because the scatter in the measurements is one-sided noise
  from "being interrupted by other work".
* The scatter appears in the `spread` column (max/min − 1). Do not trust a row
  where this is large.
* Because `amgcl` changes its preconditioning with the thread count, even the
  number of linear iterations shifts slightly. `summarize.py` detects this and
  prints a note. It means the numerator and the denominator did different
  amounts of work, so if the difference exceeds a few percent the speed-up
  needs reading with care. With `eigen` the same values come out bit for bit
  whatever the thread count, so the problem does not arise.
