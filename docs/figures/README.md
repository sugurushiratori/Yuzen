# The example-result figures

These are drawn from output the solver actually wrote, by
[`make_examples.py`](make_examples.py), from the cases in [`../../cases/`](../../cases/).
Nothing here is idealized or redrawn by hand: if a figure looks wrong, the run
behind it is wrong.

```sh
python3 docs/figures/make_examples.py --run      # run the cases, then draw
python3 docs/figures/make_examples.py            # draw from results already there
```

`--run` takes about fifteen minutes on two cores. Without it the script expects
the result directories under `runs/figures/` (or wherever `--workdir` points),
which is not kept in the repository.

| Stem | Case | Shown in |
|---|---|---|
| `leveling-1d` | `leveling_1d_reg` | README, manual §*What a result looks like* |
| `leveling-2d` | `allphys_2d_reg` | README, manual §*What a result looks like* |
| `wetting` | `wet_static` | README, manual ch. *Wetting and the contact line* |
| `concz` | `concz_2d` | README, manual ch. *Concentration and boundary conditions* |

Each stem is written in four forms, and each form has one consumer: `.png` for
this repository's README and for the HTML manual's fallback, `.pdf` for the
LaTeX manual, `.svg` for the HTML manual, which inlines it, and `.gif` for the
two cases that are worth animating. The manual reaches them as
`../figures/<stem>.pdf`; `docs/manual/tohtml.py` searches this directory as
well as `docs/manual/figures/`.

**These files are generated.** Do not edit them — edit `make_examples.py`, or
the case, and run it again.
