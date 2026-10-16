# bioconda submission recipes

These are the **submission-ready** recipes. They differ from the in-repo
`recipes/<tool>/meta.yaml` in exactly one way: the source.

| | in-repo recipe | this copy |
|---|---|---|
| `source:` | `path: ../../` (local checkout) | `url:` the tagged archive + `sha256:` |

The in-repo recipes point at the local checkout so the pixitainer image builds
with no network. bioconda cannot fetch a local path, so it needs a versioned
archive and a checksum. Everything else — the `cargo install --locked --path
tools/<tool>` build, the `--version` test, the synthetic-fixture run, the
`about:` block — is identical, and every tool commits its `Cargo.lock` so
`--locked` pins the dependency graph.

**One version per tool (D-2).** Each recipe takes its own tag, `<tool>-v<version>`,
so bumping one tool does not force a release of the other five.

## Steps

**1. Push the tags.** Only after this do the URLs resolve — before that every
one is a 404, including in `fill-sha256.sh`.

```fish
for t in fastcheck-v0.1.0 fastdedup-v1.3.0 fastprep-v0.2.0 \
         fastsim-v0.1.0 fastmaps-v0.1.0 fasthelper-v0.1.0
    git tag $t
    git push origin $t
end
```

**2. Fill the checksums.** The placeholder in every `meta.yaml` is `FILL_ME`.

```fish
./recipes/bioconda/fill-sha256.sh
```

It fetches each archive, writes the real `sha256`, and prints what it wrote. It
exits non-zero if any tag is missing, so a partial push is caught rather than
half-filled.

**3. Verify each recipe builds.** From a checkout of `bioconda-recipes`, with
the recipes copied in:

```fish
cp -r recipes/bioconda/<tool> /path/to/bioconda-recipes/recipes/<tool>
cd /path/to/bioconda-recipes
bioconda-utils build --packages <tool> --keep-old-work
```

`bioconda-utils lint` will complain about anything else it wants — most likely
`file: recipe/meta.yaml` recipe-name mismatches or a missing `doc_url` — but
the six here already carry `home`, `dev_url`, `doc_url`, `license`,
`license_file` and `recipe-maintainers`.

**4. Open one PR per package.** bioconda takes one recipe per PR; a batch of
six usually draws a request to split anyway.

## Notes for the reviewer

- **The recipes download nothing at build time.** `--locked` against a
  committed `Cargo.lock`, and the archive is a tag, not a moving ref. The
  deliberately-documented exception is the pixitainer image, which is not part
  of these recipes.
- **`fastprep` and `fastdedup` carry the Mandrillus Project attribution**
  (D-9) in their `CITATION.cff` and `LICENSE`. The other four are Raphaël
  Ribes' own work.
- **No data ships.** Every `test:` section runs `--version` and a synthetic
  fixture built inline in the recipe (D-11).
- `fastprep` is at **0.2.0**, not 0.1.0: changing a shipped default
  (`--head-dropout` 4 → 0) is a behaviour change, so it is a minor bump. The
  other five are unchanged.

## When a tool is bumped

Edit `recipes/<tool>/meta.yaml` and this copy together — the version appears in
the `{% set version %}` line and in the `--version` test. Then re-run steps
1–4 for that tool alone.
