# Yuzen manual — writing contract (read this before writing any chapter)

The manual is built with **LuaLaTeX** from `docs/manual/yuzen-manual.tex`.
Every chapter is one file under `docs/manual/chapters/`, pulled in with
`\include`. You write **body text only** — no preamble, no `\documentclass`,
no `\begin{document}`, no `\usepackage`.

## Language and voice

* **English.** The manual, the source material in `docs/*.md` and the code
  comments are all English.
* American spelling (`behavior`, `normalized`, `leveling`, `meter`), to match
  the existing `README.md`.
* Address the reader as "you". Write in the present tense.
* Be concrete. The source documents state numbers and name files; keep every
  number and every file name. Never round a measured value, never invent one.
* Do not soften a caveat. Where the source says something is **not** verified,
  or a guideline is a rule of thumb rather than a theorem, say so just as
  plainly.

## File skeleton

Each chapter file starts with a comment line, then exactly one `\chapter`:

```latex
% ch-example.tex -- one-line description
\chapter{Title in sentence case}
\label{ch:example}

Opening paragraph: what this chapter covers and who needs it.
...
```

Labels follow a fixed scheme, so that other chapters can `\cref` them:

| Kind | Pattern | Example |
|---|---|---|
| chapter | `ch:<name>` | `\label{ch:wetting}` |
| section | `sec:<chapter>:<slug>` | `\label{sec:wetting:dp}` |
| equation | `eq:<slug>` | `\label{eq:height}` |
| table | `tab:<slug>` | `\label{tab:numerics-keys}` |
| figure | `fig:<slug>` | `\label{fig:stagger}` |

Reference them with `\cref{...}` / `\Cref{...}` (cleveref supplies the word
"Section", "Table" and so on — do **not** write "Section~\ref{...}").

## Inline markup — use the semantic macros, never bare `\texttt`

| Macro | For | Example |
|---|---|---|
| `\key{...}` | an input-file key | `\key{Physics.Evaporation.Coef}` |
| `\val{...}` | a literal value written in an input file | `\val{"WM"}`, `\val{true}`, `\val{1e-6}` |
| `\file{...}` | a file or directory name | `\file{JobInfo.txt}`, `\file{cases/eq/}` |
| `\code{...}` | source identifiers, types, functions, CLI commands | `\code{SetMatrixA()}`, `\code{CRS\_Matrix}` |
| `\env{...}` | environment variables | `\env{YUZEN\_DEBUG\_FLOOR}` |
| `\reqd` | the word "required" where a default would go | |

Escape `_`, `&`, `%`, `#`, `$` inside these macros (`\code{Nx\_max}`). Inside
`lstlisting`-derived environments nothing needs escaping.

## Environments available

Code, with syntax highlighting already configured:

```latex
\begin{toml}
[Grid]
Dimension = "xh"
Nx        = 300
\end{toml}
```

`toml` (an input file), `shell` (commands), `cpp` (C++), `plain` (output, logs,
anything unhighlighted). There is also an `inp` environment for the legacy
format, but the manual does not specify that format, so do not write new
examples in it.

Admonitions — use them sparingly, and only where the source material carries a
genuine warning:

```latex
\begin{warn}[Do not judge convergence on the absolute correction alone]
  ...
\end{warn}

\begin{note}[Why the ghost layer is two cells deep]
  ...
\end{note}

\begin{tip}[Sizing a contact-line case]
  ...
\end{tip}
```

The `★` warning blocks in `docs/numerics.md` map onto `warn`. The optional
argument is the box title; keep it to one short line, in sentence case.

## Tables

For a block of input keys, always use the `keytable` environment. It is a
`longtable` with four fixed columns — key, type, default, meaning — and it
breaks across pages:

```latex
\begin{keytable}{Keys of \key{Numerics.TimeStepControl}}
\label{tab:tsc-keys}
\key{Enable}   & bool   & \val{false} & Make the time step adaptive. \\
\key{Criterion}& double & \val{0.01}  & Tolerated local truncation error. Must be $> 0$. \\
\key{MinDt}    & double & \val{1e-15} & Lower bound on the step, in seconds. The run aborts if a smaller step is requested. \\
\end{keytable}
```

Put units and the allowed range in the fourth column, as a sentence, not as a
separate column. Use `\reqd` in the default column for a key with no default.

For anything else use `yzntable`, whose second argument is the column spec:

```latex
\begin{yzntable}{Measured residual velocity}{@{}d d d d@{}}
\label{tab:resid}
\textsf{\bfseries $N_x$} & \textsf{\bfseries $\Delta x$ [\si{\micro\metre}]} &
\textsf{\bfseries residual} & \textsf{\bfseries order} \\ \midrule
1000 & 0.6  & 42.38 & \\
2000 & 0.3  & 10.79 & 1.97 \\
4000 & 0.15 & 2.707 & 1.99 \\
\end{yzntable}
```

Column types beyond the standard ones: `L{w}` `C{w}` `R{w}` are ragged
paragraph columns of width `w`, and `d` is a right-aligned column set in
lining figures (use it for every column of numbers).

## Mathematics

`amsmath` and `unicode-math` are loaded. Use `\[ ... \]` or `equation` /
`align` for display maths and `$...$` inline. Number an equation only when the
text refers back to it.

Shorthands are already defined — **use them** so the notation is uniform:

`\As \Bo \Dpn \Ma \Spn \Ev \Pen \Df \Mr \Ca \Pgx \Pgy` (dimensionless groups,
upright), `\Qvec \uvec \Pvec \er` (vectors), `\hstar \Hstar \Hmax \Hf`
`\thetae \thetacap \Pihat \epsmach \dtexp`.

So the height equation is written

```latex
\begin{equation}\label{eq:height}
  \frac{\partial H}{\partial T} = -\nabla\cdot\Qvec - \Ev\,\hat{E}(H,C) + \hat{S}.
\end{equation}
```

Units in prose go through `siunitx`: `\SI{0.2}{\micro\metre}`,
`\si{\pascal\second}`, `\SI{9.80665}{\metre\per\second\squared}`. Plain numbers
such as `1e-13` stay in maths mode as $10^{-13}$, or in `\val{}` when they are
literally what you type into an input file.

## What must stay accurate

* Every default value, range and allowed string comes from
  the parser source itself — `yuzen/src/IO/Input.cpp`, `Input.h`,
  `ConfigBlock.cpp` and `CheckParameters.cpp`. Where a prose document in
  `docs/` disagrees with the parser, **the parser wins** — and say so in a
  footnote rather than silently picking one.
* Keep all cross-references to real files in the repository (`cases/...`,
  `tools/...`, `yuzen/src/...`, `yuzen/tests/...`).
* Do not promote a rule of thumb into a law, and do not attach a precision to
  a number that the source does not give it.

## Compiling

From `docs/manual/`:

```sh
latexmk -lualatex -interaction=nonstopmode yuzen-manual.tex
```

Your chapter must compile with no undefined control sequence, no missing
reference, and no overfull box wider than 10 pt. Check your own chapter before
reporting back.

## The Japanese edition (`chapters-ja/`)

`yuzen-manual-ja.tex` builds the Japanese edition from `chapters-ja/`, one
file for each file of `chapters/` under the same name. **The English chapter
is the reference.** Change it first, then carry the change into the Japanese
file; `make check` runs `tools/manual_ja_check.py`, which fails when the two
disagree on anything a reader acts on: labels and references, the contents of
`\key \val \file \code \env \url`, every formula, every code listing, every
`\SI`/`\num` quantity and every number written in the text, the figures, and
the number of tables and boxes.

What is translated, and what is not:

* Translate the prose, the titles of chapters and sections, the captions,
  the text cells of tables, the optional titles of `warn`/`note`/`tip`, the
  footnotes and the labels of `description` items.
* Copy unchanged: every LaTeX command and environment, every `\label`, every
  `\cref`/`\Cref` argument, the argument of `\key \val \file \code \env \url`,
  all maths (inline and displayed), every code listing (input files, shell
  commands, program output — the program prints English), every `\SI{}{}`,
  the column specification of a table, and the bibliography entries.
* Every number stays exactly as written, digits and format alike: `0.02`,
  `18.94`, `30{,}000`. Never round, never convert.
* Write である調 (plain form): 「〜である」「〜する」. The reader is addressed
  implicitly, without あなた.
* Japanese punctuation: 「、」「。」, full-width （） around Japanese text,
  half-width ( ) around formulas and code. No space between Japanese and an
  inline macro or formula; LuaTeX-ja spaces them.
* Cross references: cleveref prints 第3章, 3.2節, 表3.1, 図3.1, 式(3.1),
  付録A. Write `\cref{sec:x}に示す` or `\cref{sec:x}を参照`, never
  「Section」 or 「節\ref」, and no `~` before `\cref`.
* The program's own words stay English: key names, values, file names,
  messages it prints, column headings of its output files.

Write Japanese, not English in Japanese words. The English chapter fixes
*what* is said; the Japanese sentence is free to say it in its own order and
its own grammar. Concretely:

* 無生物主語を避ける。「この仮定がNS方程式を一つの式に帰着させる」→
  「この仮定により、NS方程式は一つの式に帰着する」。「\cref{x}は〜を示す」→
  「〜は\cref{x}に示す」「\cref{x}で〜を述べる」。
* 英語の語順を崩す。結論や主題を先に置き、長い修飾は文を分ける。関係節を
  名詞の前に何重にも積まない（「〜する〜である〜の値」）。
* 一文一義。英語のセミコロン、ダッシュ（---）、コロンで繋いだ文は、日本語では
  句点で切るか「〜ため」「〜が」「〜ので」で繋ぐ。本文中の「――」は使わない。
* 代名詞を訳さない。「それは」「その」「これらの」は、なくても通じるなら削る。
  主語が文脈から明らかなら省略する。
* 名詞の羅列（「〜の〜の〜の」）を避け、動詞で言う：「時間刻みの縮小の回数」→
  「時間刻みを縮めた回数」。
* 受動態を多用しない。英語の受動態は「〜する」「〜になる」で言えることが多い。
* 英語の強調や修辞（exactly, simply, in fact, note that）を逐語訳しない。
  「まさに」「単に」「実際」「〜に注意」は必要なときだけ使う。
* 技術文書らしい定型を使う：「〜とする」「〜を用いる」「〜による」「〜に従う」
  「〜の場合」「〜のとき」「〜に限る」「〜を参照」。
* 読者への指示は論文調で書く。「〜を読む。」「〜を見る。」「〜を読むとよい。」
  「〜を見ればよい。」「〜を読むこと。」→「〜を参照されたい。」。内容の所在を
  示すだけなら「〜は\cref{x}に示す」「\cref{x}で述べる」。括弧内の
  「（\cref{x}を参照）」はそのままでよい。「〜すればよい」「〜だろう」
  「〜しよう」「〜てほしい」「ください」などの話し言葉は使わない。
* カタカナ語は定着したものに限る（ソルバー、ケース、スレッド）。訳語表にある語を使う。
* 数字・キー・式・参照の位置は文の中で動かしてよい。消したり変えたりしてはいけない
  （`make check`が見る）。

Terms (use these throughout):

| English | 日本語 |
|---|---|
| thin liquid film | 液膜 |
| lubrication approximation / long-wave | 潤滑近似 / 長波長近似 |
| film height, thickness | 膜厚 |
| free surface / substrate | 自由表面 / 基板 |
| capillary pressure / capillary leveling | 毛管圧 / 毛管レベリング |
| leveling | レベリング |
| disjoining pressure | 分離圧 |
| precursor film | 前駆膜 |
| contact line / contact angle | 接触線 / 接触角 |
| wetting / dewetting / spreading | 濡れ / 脱濡れ / 濡れ広がり |
| Marangoni (solutal) | （溶質）マランゴニ |
| evaporation / solute / solvent / concentration | 蒸発 / 溶質 / 溶媒 / 濃度 |
| viscosity / surface tension | 粘度 / 表面張力 |
| gravity (in-plane, vertical) / centrifugal force | 重力（面内、鉛直）/ 遠心力 |
| topography | 基板形状 |
| slip | すべり |
| flux / mobility / driving potential | 流束 / 移動度 / 駆動ポテンシャル |
| free energy / dissipation | 自由エネルギー / 散逸 |
| well-mixed model / concentration-in-z model | 鉛直一様モデル（WM）/ 鉛直分布モデル（CZ） |
| Newton–Kantorovich iteration | ニュートン・カントロビッチ反復 |
| Crank–Nicolson | クランク・ニコルソン法 |
| block Gauss–Seidel | ブロック・ガウス＝ザイデル法 |
| time step / adaptive time stepping | 時間刻み / 適応時間刻み |
| local truncation error | 局所打ち切り誤差 |
| residual / convergence / stagnation | 残差 / 収束 / 停滞 |
| linear solver / preconditioner / backend | 線形ソルバー / 前処理 / バックエンド |
| grid, mesh / cell / face / ghost cell | 格子 / セル / セル界面 / ゴーストセル |
| boundary condition / periodic / symmetry | 境界条件 / 周期 / 対称 |
| convection scheme / limiter | 対流スキーム / リミッター |
| conservative / conservation | 保存形 / 保存 |
| thread / reproducibility / bit-identical | スレッド / 再現性 / ビット単位で一致 |
| verification / validation | 検証 / 妥当性確認 |
| manufactured solution | 製造解 |
| order of accuracy | 精度次数 |
| steady state | 定常状態 |
| input file / key / block / default | 入力ファイル / キー / ブロック / 既定値 |
| field file / restart | 場ファイル / リスタート |
| case / test suite | ケース / テストスイート |
| rule of thumb | 目安 |

Build the edition with `make ja` and check it with `make check`.
