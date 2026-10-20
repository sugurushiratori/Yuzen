# Contributing to Yuzen

Thank you for looking. This file is deliberately short: the procedures it
would otherwise repeat are in the manual, and keeping them in one place is how
they stay correct.

## Before anything else: the one rule

**Numerical output must not change by accident.** Results are bit-identical
across OpenMP thread counts, and the test suite asserts it. A change that
moves a number is not forbidden — it is *required to be documented and
quantified*. "No case changed; `compare.py` reads `0.000e+00` throughout" is
the expected statement for a refactoring or a new term. "Results change
slightly" is not an acceptable statement.

See the manual, *Reproducibility* and *Update policy*.

## Reporting a problem

Open an issue and include:

* what you ran — the input file, or the smallest edit to a shipped case that
  reproduces it;
* the version, from `yuzen --version` or the start-up banner
  (`Yuzen-1.0.0 (build at: ...)`);
* the linear-solver backend and `OMP_NUM_THREADS`, both of which are recorded
  in `JobInfo.txt`;
* the whole of `JobInfo.txt`, and `FatalErrors.txt` if it is not empty.

A wrong answer and a crash need the same information. A *slow* run needs one
thing more: the `###### Performance ######` section of `JobInfo.txt`, which
gives the time-step loop alone with the Newton and linear iteration counts.

## Proposing a change

1. Build and run the test suite: `ctest --test-dir build`. All three
   thread-count entries must pass.
2. Run the shipped cases in `cases/` before and after your change, and compare
   them with `tools/compare.py` — **fixing the backend and the thread count
   across the comparison**, so that you are measuring your change and not the
   rounding of a different linear solver.
3. Say in the change description whether the numbers moved, and by how much,
   on which case.
4. Update the manual (`docs/manual/`) in the same change, and the companion
   documents `docs/theory.md`, `docs/numerics.md` and `docs/validation.md`
   where they touch the same point. The manual is the authoritative
   description; a change that lands in the code but not in the manual leaves
   the two out of step.
5. Add an entry to the revision history, in the form the manual's
   *Revision history* appendix specifies, and the same line to
   `CHANGELOG.md`.

Making a release is a separate checklist, kept in
[`docs/releasing.md`](docs/releasing.md).

The manual's *Extending the model* chapter is the full procedure for adding a
physics term, a boundary condition, a shape function, a material-property
dependence or a linear-solver backend — including the four things in the code
that keep results bit-identical, which constrain how a term may be written.
**Read that chapter before adding physics.** It is in `docs/manual/`; the PDF
and the searchable HTML edition are attached to each release.

## Style

* C++17. Follow the surrounding file rather than any external guide.
* Everything is in English: the comments, `docs/` and the manual. Public
  methods carry a Javadoc block; `docs/Doxyfile` builds the API reference and
  is expected to stay free of warnings.
* A comment explains *why*, not *what*. The comments that earn their place in
  this code are the ones recording a trap somebody already fell into.

## License

By contributing you agree that your contribution is licensed under the
Apache License 2.0, as the rest of the project is.
