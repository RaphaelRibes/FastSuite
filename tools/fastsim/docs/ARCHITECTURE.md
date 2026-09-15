# Architecture

```
profile JSON ──► profile.rs ──► typed CDFs
                                    │
cli.rs ──► Config ──► main.rs loop  │
                    │                ▼
                    │        randutil.rs (seeded ChaCha8)
                    │                │
                    │   ┌────────────┼──────────────┐
                    │   ▼            ▼              ▼
                    │ biology.rs  artefact.rs   quality.rs
                    │ (insert)    (h/t/adapters) (Q + errors)
                    │   │            │              │
                    │   └─────► build_read() ──────┘
                    │                │
                    ▼                ▼
              fastq.rs (R1/R2)   truth.rs (gt.tsv)
```

## Data flow per pair

1. **Insert** — `biology.rs`: length `I ~ N(insert_mean, insert_sd)` (floor
   `insert_min`); sequence = random DNA at GC% or a slice of a reference
   FASTA contig. The insert is the true double-stranded fragment; R1 reads
   its 5′ end, R2 the reverse complement of its 3′ end. With probability
   `dup_rate` the pair is a **PCR duplicate**: it re-sequences a fragment
   stored in the `History` ring buffer (only fresh inserts enter the pool),
   so the template is real and earlier in the file. With probability
   `drop_rate` (mutually exclusive with dup) the pair is a **must-drop**:
   `I ~ U[1, drop_min_len)` — a fragment too short to survive trimming.
   Reads stay exactly `read_len` bases in every class (adapter readthrough
   fills the window); only the GT keep window shrinks.
2. **Artefact decision** — `artefact.rs::decide`:
   - sample `(t1, t2)` from the **joint** R1×R2 tail CDF (correlated tails —
     the key finding), then `h1`, `h2` from the head CDFs.
   - if `h + t >= read_len` on either mate, the pair is entirely artefact →
     `full_flag=1` (the correct answer is to drop it). Full-length reads are
     NOT sampled from a separate pair-level spike: the profile's head AND
     tail histograms already contain them (an all-G read has both a 150 bp
     head and a 150 bp tail), so the marginal CDFs reproduce them
     per-mate. This matters because the measured full-length rate is
     R2-only (mandril: R1 = 0.00000, R2 = 0.00217) — a pair-level spike
     would fabricate all-G R1 reads that don't exist in the data.
   - Artefacts are drawn fresh even for duplicates: dropout is a per-cluster
     failure, so a copy gets its own polyG runs.
3. **Read construction** — `build_read`: `[h Gs][biology][adapter?][t Gs]`,
   exactly `read_len` bases. Adapter appears only when `I < read_len − h − t`
   (readthrough), cycling the adapter sequence.
4. **Quality** — `quality.rs::qualify`:
   - per-read baseline `q_base` sampled from `read_meanq` CDF; offset
     `δ = q_base − mean(pos_q)`;
   - `q_i = pos_q[i] + δ`, clamped [0, 41];
   - bases inside a G-run lose `tail_q_collapse` Phred;
   - substitution error per base at `10^(−q/10)` (uniform over other 3
     bases). The emitted base reflects the error → self-consistent FASTQ.
   - full-length reads use flat `full_length_q` (all bases are artefact).
5. **Write** — `fastq.rs` (NovaSeq-style headers, per-pair unique
   coordinates, Phred+33), `truth.rs` (one GT row per pair).

## Ground-truth semantics

The GT row `(insert, h1, t1, h2, t2, full_flag)` fully determines the correct
trimmer output per mate:

```
[0, h)                5′ G-run            → trim
[h, min(I, L−t))      biology             → keep
[min(I, L−t), L−t)    adapter (I < L−t)   → trim
[L−t, L)              3′ G-run            → trim
```

`full_flag=1` → drop the pair. A scorer compares any trimmer's output window
to `[h, min(I, L−t))`; length-only scoring is insufficient (a tool can remove
the right number of bases from the wrong end).

## Determinism

ChaCha8 seeded from `--seed`; every draw flows from that stream, so the same
seed + profile + args produce byte-identical outputs. The validation sbatch
re-profiles simulated reads with `polyg-profiler` to check the closed loop.
