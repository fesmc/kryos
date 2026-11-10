# kryos documentation

The site [fesmc.github.io/kryos](https://fesmc.github.io/kryos/) is rendered from the Markdown pages in this folder with
Quarto (`_quarto.yml`; the workflow `.github/workflows/docs.yml` publishes it on every push to `main`). Preview locally
with `quarto preview docs`.

Not rendered, kept as written:

- [design.md](design.md): the design proposal (vocabulary, architecture, contracts, exchange, regional vs global,
  time, the decisions).
- [dev/](dev/): design notes, plans, gate reports, backend notes (e.g. [dev/plastic-backend.md](dev/plastic-backend.md),
  [dev/vilma2-paleomist.md](dev/vilma2-paleomist.md)).
- [study/](study/): the codebases read for the design, with file:line references as of 2026-10-24.
