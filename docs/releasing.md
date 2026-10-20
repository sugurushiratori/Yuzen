# Releasing Yuzen

The checklist for the maintainer, in the order the steps are done. A release
is a tag `vMAJOR.MINOR.PATCH` on a commit for which every step below has
passed; the release workflow (`.github/workflows/release.yml`) then builds
the manual and publishes the release text. What the version number promises
is set out in the manual, *Versions and updates*.

## 1. The version

Set the new version in every place that carries it:

- [ ] `yuzen/src/Version.h`: `YUZEN_VERSION_Major`, `_Minor`, `_Rev`.
- [ ] `CITATION.cff`: `version`. Add `date-released` on the day of the tag.
- [ ] `docs/Doxyfile`: `PROJECT_NUMBER`.
- [ ] `yuzen/src/yuzen.manifest`: `assemblyIdentity version` (four fields).
- [ ] For a minor release: add the new `Yuzen-MAJOR.MINOR` series to the
      accepted list in `Input::ReadRoot()`. Leave the old entries in place
      unless the meaning of a key has changed (manual, *Versions and
      updates*).

`yuzen --version` must then print `Yuzen-MAJOR.MINOR.PATCH`.

## 2. The release text

- [ ] Add a section `## [MAJOR.MINOR.PATCH]` to `CHANGELOG.md`. The release
      workflow publishes that section as the text of the release, and it
      stops if the section is missing.
- [ ] Add the same entry to the manual's *Revision history* appendix
      (`docs/manual/chapters/app-revisions.tex`), in the form that appendix
      specifies. A change that moved a number names the case and the size of
      the difference. A release in which nothing moved says so.

## 3. Checks on the release commit

Run these on a clean checkout of the commit that will be tagged:

- [ ] `ctest` on a Release build with the Eigen backend
      (`-DYUZEN_WITH_AMGCL=OFF`) and on one with AMGCL. Both pass all
      entries, including the 1- and 4-thread re-runs.
- [ ] The CI workflow passes on every runner of its matrix (Linux, macOS,
      Windows) and the sanitizer job passes.
- [ ] `cd docs/manual && make && make check`: no undefined reference, no
      overfull line, and both editions build.
- [ ] `python3 tools/doc_keys.py --exe build/bin/yuzen --workdir /tmp/dk`
      reports no mismatch between the manual's key tables and the parser.
      (`ctest` runs it too.)
- [ ] `reuse lint` reports the project compliant.
- [ ] Doxygen (`cd docs && doxygen Doxyfile`) runs free of warnings.
- [ ] If the numerical output changed in this release:
  - [ ] Regenerate the example figures with
        `python3 docs/figures/make_examples.py --run`.
  - [ ] Bring the numbers quoted in `docs/validation.md` and the manual's
        *Verification* chapter up to date by re-running the commands each
        section gives.
- [ ] If the run time of the shipped cases is to be quoted anew, measure it
      with `tools/case_runtimes.py` (and the parallel scaling with
      `tools/scaling/`). Then redraw `docs/manual/figures/parallel-scaling.*`
      with `docs/manual/figures/make_scaling.py`.
- [ ] Read the diff since the previous release for anything that should not be
      public: credentials, personal data other than the author's contact,
      machine-local paths, unpublished data.

## 4. Tag and publish

To try the release workflow first, push a tag with a suffix, for example
`v1.0.0-test`. It builds and checks the manual and publishes a pre-release
with the four files; delete that release and its tag
(`git push origin :v1.0.0-test`) afterwards.

```sh
git tag -a vMAJOR.MINOR.PATCH -m "Yuzen MAJOR.MINOR.PATCH"
git push origin vMAJOR.MINOR.PATCH
```

- [ ] The release workflow passes. The release then carries
      `yuzen-manual.pdf`, `yuzen-manual.html` and the text from
      `CHANGELOG.md`, and GitHub adds the source archives.
- [ ] Download the source archive, build it and run `ctest -R yuzen_tests`.
      This checks that the archive is complete and that Eigen can be fetched
      when it is not installed.
- [ ] With the Zenodo integration enabled, Zenodo archives the release from
      `.zenodo.json` and assigns a DOI. On the first release, add the
      concept DOI (the one that resolves to the latest version) to
      `CITATION.cff` as `doi:` and to the README's *Citing Yuzen* section.
