# WCN Lab Assignment: MATLAB Experiments 1–15

This is a LaTeX submission of the Wireless and Cellular Networks Laboratory manual (23ET7PCWCN). It covers all 15 experiments in the manual's order.

## Project layout

```
config.m                 <- the ONLY place the USN is defined
run_all.m                <- runs every experiment and saves the outputs
experiments/             <- exp01.m ... exp15.m plus helper functions
output_images/           <- generated PNGs: <USN>_EXPnn_FIGk.png, <USN>_EXPnn_OUTPUT.png
generated/               <- generated files (usn.tex, expNN_code.m, expNN_outputs.tex)
main.tex                 <- LaTeX source of the report
WCN_Lab_Assignment.pdf   <- compiled report
```

## 1. Set your USN (one place only)

Open `config.m` and change this line:

```matlab
USN = "YOUR_USN_HERE";
```

Nothing else holds the USN. `run_all.m` reads it from `config.m` and uses it for:

* every output-image filename (`<USN>_EXP01_FIG1.png`, `<USN>_EXP01_OUTPUT.png`, …),
* the PDF title page and page headers, through `generated/usn.tex`.

The images for the old USN are deleted automatically when you re-run.

## 2. Run the experiments and generate the images

**MATLAB:** open the project folder in MATLAB and run:

```matlab
run_all          % all 15 experiments
run_all([3 7])   % only experiments 3 and 7
```

**GNU Octave (the version used to build the included PDF):**

```bash
octave --no-gui --eval "run_all"
```

Use `octave` (Qt graphics), not `octave-cli`. The `octave-cli` binary only has the FLTK/gnuplot toolkits and draws poorer figures. A full run takes about 5–10 minutes. The slowest experiments are 5, 6, 10, 13 and 14, which are Monte Carlo simulations.

For each experiment, `run_all`:

1. runs `experiments/expNN.m` unchanged,
2. records the Command Window output, saves it as `output_images/<USN>_EXPNN_OUTPUT.txt` and renders it as `<USN>_EXPNN_OUTPUT.png`,
3. saves every figure at 200 dpi as `<USN>_EXPNN_FIGk.png` (shown in the PDF as Figure N.k),
4. writes `generated/expNN_outputs.tex`, which `main.tex` includes.

It also writes `generated/expNN_code.m`: the program followed by the helper functions it calls, merged into the single code listing shown in the PDF.

You can also run a single script on its own, e.g. `cd experiments; exp04`. That only shows the figures. It does not save them.

## 3. Compile the PDF

```bash
pdflatex -jobname=WCN_Lab_Assignment main.tex
```

```bash
pdflatex -jobname=WCN_Lab_Assignment main.tex
```

Run it twice so the table of contents is filled in. You can also use `latexmk -pdf -jobname=WCN_Lab_Assignment main.tex`. The packages it needs are listings, titlesec, fancyhdr, caption, float, chngcntr, beramono, lmodern and hyperref, all in standard TeX Live/MiKTeX.

**Order matters:** change the USN → `run_all` → `pdflatex`.

## Notes on the transcription

* The manual is a scanned document. The code was transcribed from high-resolution renders of the scans. Every program was then executed, and the results were compared with the "check values" printed in the manual.
* `prctile_local.m` and `ternary_local.m` are called by Experiments 9 and 12 but **are not printed in the manual**. Minimal implementations were added in `experiments/`.
* The only edit to the printed code is in Experiments 8 and 12. The table-heading strings `'99 %% BW (kHz)'` and `'PAPR 99 %% (dB)'` were passed as `%s` arguments and printed a literal `%%`, so they were changed to a single `%`. No computation was changed.
* Experiment 12 uses `sinc`. In MATLAB that function belongs to the Signal Processing Toolbox (Octave has it built in).
* The PDF in this folder was produced with **GNU Octave**, since MATLAB was not installed. The title page states the execution environment. Monte Carlo values depend on the random number generator, so a MATLAB run will give slightly different simulated numbers.
