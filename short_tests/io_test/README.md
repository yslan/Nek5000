# io_test — Nek5000 I/O / restart CI test

Exercises the checkpoint **read** path (`mfi`, `mfi_getv`, `mfi_gets` in
`core/ic.f`) and the field **write** path (`prepost`) across the many
restart/IO option combinations, plus the `full_restart` round-trip.

## How to run

From `short_tests/`:
```
NEK_SOURCE_ROOT=<repo>/Nek5000 FC=mpif77 CC=mpicc MPI=1 PARALLEL_PROCS=2 \
EXAMPLES_ROOT=<repo>/Nek5000/short_tests TOOLS_BIN=<repo>/Nek5000/bin \
LOG_ROOT=<an existing dir> \
python3 -m pytest NekTests.py -k IO_Test -v
```
Each sub-test writes `PASSED`/`FAILED` to the run log; `userchk` prints
**`All I/O tests PASSED`** iff every sub-test passed, and `NekTests.py` asserts
that phrase. (Do not change those literal strings — the harness greps them.)

## Test variants (`NekTests.py::IO_Test`)

Seven variants (all run at the default whole-field/single-round batching unless
noted — `userParam06=nrst_rd`/`userParam07=nrst_cm` are reset to 0 in every non-
multi-round test so they exercise the default path, not a leftover bounded config):

| Method | np | What it forces |
|---|---|---|
| `test_PnPn2_Parallel`      | 2 | Crystal-router restart (default); runs the full sub-test suite over hrefine 1/2/3/2,2, then the `full_restart` round-trip (PnPn-2 and PnPn). |
| `test_PnPn2_Parallel_RMA`  | 2 | MPI-RMA restart (`ifcrrs=.false.` via `userParam02=1`), hrefine 1/2/3/2,2. |
| `test_PnPn2_Serial`        | 1 | np==1 direct-read bypass, hrefine 1/2/3/2,2. |
| `test_PnPn2_Parallel_HighOrderFP64` | 2 | Writes an FP64 order-10 checkpoint (`userParam03=1`), reads it into an lx1=6 build (`userParam03=2`) — the PR #900 high-order acceptance case (CR). |
| `test_PnPn2_Parallel_RMA_HighOrderFP64` | 2 | Same PR #900 high-order FP64 read, but via the MPI-RMA bounded window (`userParam02=1`). |
| `test_PnPn2_Parallel_MultiRound` | 2 | Forces batched CR (`nrst_rd=3`, `nrst_cm=2`, `loglevel=3`) on hrefine 1/2/2,2; asserts max rounds/batch == 2. |
| `test_PnPn2_Parallel_RMA_MultiRound` | 2 | Same multi-round stress on the MPI-RMA path. |

Toggles (wired in `io_test.usr:usrdat`/`userchk`): `userParam02` 0=crystal
router / 1=MPI-RMA; `userParam03` 1=hio write / 2=hio read; `userParam06`→`nrst_rd`
(read-batch), `userParam07`→`nrst_cm` (per-round receive cap); `[MESH] hrefine`.

## What `io_test.usr` does

`userchk` dispatches one of three modes (each ends in `call exitt`):
`userParam03>=1` → `hio_highorder_test`; `hrefine>0` → h-refine read tests;
else → the normal suite (`read_old_file_test`, `write_read_mpiio_file_test`,
`write_read_multipe_file_test` [np even], `write_read_fld_file_test`,
`read_gfldr_test`). Fields compared: mesh, velocity, pressure, temperature, 2
passive scalars — analytic reference cached by `fields_init`, compared by
`fields_compare` (L2 norm).

Shared helpers (post-cleanup): `report` (pass/fail tail), `chk_restart_hdr`
(restart-header verification via a 14-char flag spec), `run_write` (the write
prologue), and `REFDATA` (the `/REFDATA/` reference-field cache header).

## Fixtures (committed checkpoint inputs)

| File | Format |
|---|---|
| `m12io_test0.f00001` | single-prec, 1 file, order 6 |
| `m22io_test{0,1}.f00001` | single-prec, 2 files, order 10 (→ interp to 6) |
| `m22io_test{0,1}.f00002` | single-prec, 2 files, byte-swapped |
| `m22io_test{0,1}.f00003` | double-prec, 2 files, byte-swapped |
| `fb2io_test.fld01` | binary `.fld` (p67=4) |
| `fa2io_test.fld01` | ASCII `.fld` (p67=-4) |
| `io_test.re2`, `io_test_rs.re2` | 36-element 3D meshes |

Everything else in this directory (`nek5000`, `obj/`, `*.log.*`, run-produced
`io_test0.f*`, `hr2io_test*`, `hio*`, `io_test.fld*`, `SIZE`, `.ma2`) is
generated and not committed. `NekTests.py::_clean_generated_flds` removes the
run-produced checkpoints at each method start.

## Coverage matrix

| Axis | Value | Covered |
|---|---|---|
| Reader | MPIIO single-file | ✅ |
|        | non-MPIIO multi-file | ✅ (np even) |
| Redistribution | crystal router (`ifcrrs=T`) | ✅ |
|                | MPI-RMA (`ifcrrs=F`) | ✅ |
|                | np==1 bypass | ✅ |
| Order | copy (`nxr==lx1`) | ✅ |
|       | interp (`nxr≠lx1`) | ✅ |
| wdsize | 4 / 8 | ✅ / ✅ |
| hrefine | off / 2 / 3 / 2,2 | ✅ |
| Fields | mesh/vel/pressure/temp/2 PS | ✅ |
| Format | std `.f` / multi-file / bin `.fld` / ASCII `.fld` / gfldr | ✅ |
| Pressure | `if_full_pres` (PnPn / PnPn-2) | ✅ |
| Full restart | `full_restart` round-trip | ✅ |

### Known gaps (not yet covered)
- **MHD B-field restart** — no `bx/by/bz` read or compare (`ifmhd` only skips
  that field in `useric`). Main gap.
- **Multi-file read on odd `np`** — `write_read_multipe_file_test` is gated on
  `mod(NP,2)==0`, so the serial (np==1) path never exercises multi-file.
- **Dedicated `iskip` case** — present-but-skipped fields are exercised only
  implicitly via the restart-flag checks, not a dedicated assertion.
- **Passive scalars in `full_restart`** — the rs case is velocity+temperature
  only (`io_test_rs.par` has no `[SCALAR*]`).

(The `io_test.usr` readability refactor and coverage history are recorded in the
project's design notes, kept outside this source tree / in the PR discussion.)

## Testing the bounded-receive CR path (batches / rounds)

The crystal-router restart splits the read into **batches** (`nrst_rd`) and each
batch's redistribution into **rounds** bounded by a receive cap (`nrst_cm`, runtime;
0 = auto). Two knobs (wired in `io_test.usr:usrdat`): `uparam06 -> nrst_rd`,
`uparam07 -> nrst_cm`. Set `general:loglevel = 3` to print, per `mfi_gets/getv`
call, two verbose lines: batch/round shape and buffer sizes.

**Run the dedicated multi-round CI test** (hrefine 1 then 2,2, `nrst_rd=3`,
`nrst_cm=2` → `nrounds=2`; asserts fields pass AND max rounds == 2):
```
cd short_tests
NEK_SOURCE_ROOT=<repo>/Nek5000 FC=mpif77 CC=mpicc MPI=1 PARALLEL_PROCS=2 \
EXAMPLES_ROOT=<repo>/Nek5000/short_tests TOOLS_BIN=<repo>/Nek5000/bin \
LOG_ROOT=<an existing dir> \
python3 -m pytest NekTests.py -k "IO_Test and MultiRound" -v
```
(drop `and MultiRound` to run all seven IO_Test variants.)

**See the results** in the run log:
```
grep -E 'rounds/batch|cap=' io_test/io_test.log.2.pn_pn_2.parallel
```

**Demo output** (`nrst_rd=3`, `nrst_cm=2`, `loglevel=3`; each field read prints 4
lines — per-call stats, then the `avail`/`batch`/`round` sizing review.
`ifcrrs= T` is the crystal-router path, `F` is MPI-RMA):
```
 mfi_getv: nbatch= 6  rounds/batch min/max/avg= 2 2 2.0  recvmax= 3  ifcrrs= F
 mfi_getv  avail: wk= 6976368  w2= 9966240   need/elem: nxyzr= 648  li= 650
 mfi_getv  batch: nbe= 3  (nrst_rd= 3  w2fit= 15380  wkfit= 10749 )
 mfi_getv  round: cap= 2  bound by nrst_cm= 2
 mfi_gets: nbatch= 6  rounds/batch min/max/avg= 2 2 2.0  recvmax= 3  ifcrrs= F
 mfi_gets  avail: wk= 6976368  w2= 9966240   need/elem: nxyzr= 216  li= 218
 mfi_gets  batch: nbe= 3  (nrst_rd= 3  w2fit= 46140  wkfit= 32149 )
 mfi_gets  round: cap= 2  bound by nrst_cm= 2
 ...
 PASSED
```
Reading: `nbatch=6` (the field's `nelr` split into 6 read batches of `nbe=3`),
`rounds/batch = 2` (each batch's receive split into 2 rounds of `cap=2`),
`recvmax=3` (largest per-rank incoming per batch). The `avail` line is the
ceilings (`wk=2*lwk`, `w2=lrbs`) vs per-element need (`nxyzr`, tuple `li`); the
`batch` line shows which limit binds `nbe` (`nrst_rd` / `w2fit` / `wkfit`); the
`round` line shows the receive `cap` and its binding limit. `PASSED` confirms
the fields are byte-clean across the split. At the default `loglevel=2` these
lines are silent.

**Try it by hand** on any case: in the `.par` set `[GENERAL] loglevel = 3`,
`userparam06 = <nrst_rd>`, `userparam07 = <nrst_cm>`, then run and grep as above.
`nrst_rd` small ⇒ more batches; `nrst_cm` small ⇒ more rounds/batch.
