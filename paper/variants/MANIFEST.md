# Reconstructed source variants

Each version below is the exact state of the file after one recorded
Write or Edit, replayed from the session transcript. `verified` means the
final replayed state matches the file currently on disk byte for byte,
which makes the whole chain for that file faithful.

## projects\quasiperfect\paper\export_conversation.py — verified

| v | time (UTC) | operation | sha1 | label |
|---|---|---|---|---|
| 1 | 2026-09-30T05:30:24 | Write (full file) | `3570ecab2217` |  |
| 2 | 2026-09-30T06:05:31 | Edit | `167b9daf5ff8` |  |
| 3 | 2026-09-30T06:05:33 | Edit | `a771fe2a8ddd` |  |

## projects\quasiperfect\paper\reconstruct_variants.py — verified

| v | time (UTC) | operation | sha1 | label |
|---|---|---|---|---|
| 1 | 2026-09-30T06:09:40 | Write (full file) | `22a6f3cb9b0e` |  |
| 2 | 2026-09-30T06:11:48 | Edit | `89155eeda1b3` |  |
| 3 | 2026-09-30T06:13:14 | Edit | `c6d4394ecf4a` |  |
| 4 | 2026-09-30T06:13:18 | Edit | `c0078854475d` |  |
| 5 | 2026-09-30T06:13:26 | Edit | `e31adf68b7db` |  |
| 6 | 2026-09-30T06:17:32 | Edit | `560011e8a2d4` |  |

## projects\quasiperfect\paper\verify_variants.py — verified

| v | time (UTC) | operation | sha1 | label |
|---|---|---|---|---|
| 1 | 2026-09-30T06:15:40 | Write (full file) | `4e6be5228548` |  |
| 2 | 2026-09-30T06:17:29 | Edit | `ff8b80f348b4` |  |

## projects\quasiperfect\peek.py — verified

| v | time (UTC) | operation | sha1 | label |
|---|---|---|---|---|
| 1 | 2026-09-29T19:20:13 | Write (full file) | `4ef4197f638d` |  |

## projects\quasiperfect\qp.py — verified

| v | time (UTC) | operation | sha1 | label |
|---|---|---|---|---|
| 1 | 2026-09-29T18:58:07 | Write (full file) | `924feb5443cc` | **Issue 1, faulty**: `last_prime` bisects only past the maximum and flags its rising branch as unresolved instead of solving it |
| 2 | 2026-09-29T18:58:20 | Edit | `a6e73cba01ac` | `record` no longer re-factorises to re-check; callers verify |
| 3 | 2026-09-29T18:58:51 | Edit | `bb7f56533ce9` | the `lo >= 2` branch enumerates its finite prime range (Issue 2) |
| 4 | 2026-09-29T18:58:54 | Edit | `70f3b2b842c9` | coprimality generalised from c = 1 to any c (Issue 2) |
| 5 | 2026-09-29T18:58:55 | Edit | `72bcc11bc6eb` | docstring for the coprimality bound corrected |
| 6 | 2026-09-29T19:12:23 | Edit | `5afd73661794` | faulty rising-branch code removed from `last_prime` (Issue 1) |
| 7 | 2026-09-29T19:12:30 | Edit | `d6bd539c38fd` | **Issue 1, fixed**: `poly_roots` bisects every monotone stretch |
| 8 | 2026-09-29T19:20:10 | Edit | `0928f4dce930` | records the current state, for the branch-sampling diagnostic |
| 9 | 2026-09-29T19:20:38 | Edit | `394379eba7d6` | **Issue 4, fixed**: refine open exponents before the last prime |
| 10 | 2026-09-29T19:20:39 | Edit | `a6f639e78882` | that threshold exposed as a parameter |
| 11 | 2026-09-29T19:23:18 | Edit | `c9e0726e7c0d` | **Issue 5, rejected**: same rule at every r, ten times slower |
| 12 | 2026-09-29T19:23:19 | Edit | `fac673b91cfb` | ratio parameter added to the rejected rule |
| 13 | 2026-09-29T19:25:15 | Edit | `8e6b002dba8a` | rejected rule adjusted |
| 14 | 2026-09-29T19:25:17 | Edit | `ff1c24c85c77` | **Issue 5, reverted** to the r = 1 rule only |
| 15 | 2026-09-30T05:16:29 | Edit | `8aa22972cef1` | **Issue 8, faulty**: early refinement keyed on `hi` alone, which does not decrease, so refinement never terminates. Not runnable by itself: it uses a parameter added only in v16 |
| 16 | 2026-09-30T05:16:31 | Edit | `de3929cc7906` | the same defect, now runnable: threshold exposed as a parameter |
| 17 | 2026-09-30T05:31:21 | Edit | `8d699b2480e6` | **Issue 8, fixed**: keyed on both ends of the range, so the condition eventually fails and the branch ends |

## projects\quasiperfect\run_parallel.py — verified

| v | time (UTC) | operation | sha1 | label |
|---|---|---|---|---|
| 1 | 2026-09-29T23:59:16 | Write (full file) | `d5bd9b93be99` |  |

## projects\quasiperfect\split_task.py — verified

| v | time (UTC) | operation | sha1 | label |
|---|---|---|---|---|
| 1 | 2026-09-30T03:46:25 | Write (full file) | `86877fb5a314` |  |

## projects\quasiperfect\test_leaves.py — verified

| v | time (UTC) | operation | sha1 | label |
|---|---|---|---|---|
| 1 | 2026-09-29T19:09:33 | Write (full file) | `9b3e2f6142b4` |  |

## projects\quasiperfect\test_qp.py — **UNVERIFIED**

The replayed chain does not reproduce the file now on disk. The usual cause is an edit made through the shell rather than the editing tool: the transcript records the command but not the resulting text, so the chain has a gap. The versions below are still exact up to that gap, and the final row is the file as it stands on disk.

| v | time (UTC) | operation | sha1 | label |
|---|---|---|---|---|
| 1 | 2026-09-29T18:58:18 | Write (full file) | `bc2621d3b623` |  |
| 2 | 2026-09-29T18:59:04 | Edit | `dc7e9445673c` | phase 2 appended; later truncated through the shell, which is why this chain does not verify |
| current | – | file on disk | `09b973661f7c` | the version actually used |

