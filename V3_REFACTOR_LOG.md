# DataCooker v3 Refactor Log

This log tracks the v3 quality refactor. Each step records **what / why / how /
verification**. Scope was chosen to fix functional weaknesses in the core
engine; tests, consumer (StructCooker) cleanup, and user-facing docs are
explicitly out of scope for this pass (handled separately).

Because P2 changes the public API surface, the overall release is **v3.0**
(P0+P1 alone would have been v2.2).

Baseline before this work: `32 passed, 9 skipped`.

| Step | Tier | Title | Status |
|------|------|-------|--------|
| 1 | P0 | Wildcard arg resolution determinism | ✅ done |
| 2 | P0 | Mermaid line breaks (`\n` → `<br/>`) | ✅ done |
| 3 | P1 | Visualization verbosity levels | ✅ done |
| 4 | P1 | Memoize graph analysis | ✅ done |
| 5 | P1 | Unify typing helpers + tighten annotations | ✅ done |
| 6 | P2 | Materialize explicit DAG | ✅ done (with #1) |
| 7 | P2 | Decompose Reader/Writer hooks | ✅ done |
| 8 | P2 | Redesign `with_namespace` | ✅ done |
| 9 | hardening | Runnable LMDB env + self-reference cycle fix | ✅ done |
| 10 | hardening | Concurrency contract documented | ✅ done |
| 11 | hardening | Remove legacy compatibility cruft / dual-naming | ✅ done |
| 12 | hardening | Warn on wildcard 0-match | ✅ done |

---

## Step 5 (P1) — Unify typing helpers + tighten annotations

**Why.** Annotation logic was duplicated: `recipe._allows_none` /
`recipe._is_supported_annotation` plus an inlined `type(None) in get_args(...)`
in `executor.resolve`. `_is_supported_annotation` also had a too-permissive
`bool(get_args(annotation))` fallback that accepted malformed annotations.

**How.**
- New module `utils/typeutils.py` with `allows_none()` and
  `is_supported_annotation()` as the single source of truth.
- `allows_none()` now also treats bare `None` / `NoneType` as nullable.
- `is_supported_annotation()` dropped the loose `bool(get_args(...))` branch;
  it accepts only `Any`, concrete classes, and constructs with a real origin
  (`list[int]`, `str | None`, `Optional[X]`). Invalid declarations now fail fast.
- `recipe.py` and `executor.py` import these; local duplicates removed and the
  now-unused `get_args` / `get_origin` imports dropped.

**Verify.** `32 passed, 9 skipped`.

## Step 2 (P0) — Mermaid line breaks

**Why.** Mermaid node labels used a literal `\n` (`f"step {i}\\n{...}"`), which
Mermaid does **not** render as a line break — labels collapsed to one line with
a literal `\n`. (Graphviz DOT's `\n` is valid and was left unchanged.)

**How.** In `RecipeBook.to_mermaid`, the step label now joins with `<br/>` and
escapes the describe() text explicitly. DOT path untouched.

**Verify.** Updated the one assertion in `tests/test_core.py` that pinned the
old buggy `\n` output to expect `<br/>`. `32 passed, 9 skipped`.

## Steps 1 + 6 (P0 + P2) — Wildcard determinism via a materialized DAG

These two were implemented together because they are the same change: the
wildcard bug *is* the absence of a real graph edge.

**The bug.** Wildcard dependencies (e.g. an arg `("part_*", str)`) were skipped
entirely during graph analysis (`if _is_wildcard_pattern: continue` in both
`execution_order` and `required_inputs`). So a wildcard consumer had **no edges**
to its producers. At runtime `_expand_wildcard_args` then matched whatever keys
*happened* to already be in the `ExecutionContext` — i.e. the result depended on
which other targets had been resolved first. Producers might not run at all.

**The fix — materialized adjacency (single source of truth).**
- `RecipeBook.adjacency()` builds `target -> [internal dependency targets]`
  once, with wildcards statically expanded via `match_targets(pattern,
  exclude=<own outputs>)` (sorted → deterministic; a step never matches its own
  outputs, preventing self-cycles).
- `execution_order()` now does its topo-sort / cycle detection purely over
  `adjacency()`, so wildcard producers are ordered before their consumer.
- `required_inputs()` recurses into wildcard-expanded producers and reports only
  genuine external inputs.
- `RecipeBook.match_targets()` is exposed (public) for the executor.

**Runtime side (executor).** `_resolve_recipe_arguments` now resolves every
declared target a wildcard binds to (via `match_targets`) *before* calling
`_expand_wildcard_args`. The matched set is therefore fixed regardless of
resolution order. The executor keeps its own recursive cycle guard as a runtime
safety net for `validate=False` callers; the authoritative cycle/order analysis
lives in `RecipeBook.adjacency()`.

**Verify.** `32 passed`. Ad-hoc check (`part_a`/`part_b` → wildcard `part_*`
collector): order = `[part_a, part_b, collected]`, `required_inputs == {seed}`,
`adjacency()["collected"] == {part_a, part_b}`, result deterministic.

**Known follow-up.** Graph *rendering* (`to_mermaid`/`to_dot`) still draws a
wildcard dependency as a single raw `part_*` pseudo-node rather than expanded
edges; addressed cosmetically in step 3.

## Step 4 (P1) — Memoize graph analysis

**Why.** `adjacency()` / `execution_order()` / `required_inputs()` were
recomputed from scratch on every `describe()` / `validate()` / `visualize()` /
`execute()` call — and after step 1, `execution_order` builds the whole
adjacency each time. For a "static graph engine" that is wasteful.

**How.**
- Three caches on `RecipeBook`: `_adjacency_cache` (structure-only) and
  `_order_cache` / `_required_cache` keyed by the **normalized** requested-target
  tuple.
- Single invalidation point: `_register_step()` calls `_invalidate_caches()`.
  Keying by the normalized tuple means `set_default_targets()` needs no reset
  (it only changes what `None` normalizes to, not graph structure).
- Public `adjacency()` returns a deep copy; internal `_adjacency()` serves the
  cached dict. `execution_order()` / `required_inputs()` return copies of their
  cached results so callers can't corrupt the cache.

**Verify.** `32 passed`. Ad-hoc check: cache populates on first call, survives
external mutation of the returned `adjacency()` copy, and is invalidated when a
new step is registered (re-derives `[a, b]`).

## Step 3 (P1) — Visualization verbosity + wildcard-edge rendering

**Why.** Step labels embedded the full `recipe.describe()` (every kwarg binding
and static param inline). On real recipes (e.g. the ccd ingest) a single node
dumped an 8-element column list — unreadable. Also, after step 1 the rendered
graph still drew wildcard deps as a single raw `part_*` pseudo-node, not the
real edges.

**How.**
- `Recipe.describe(detail=...)`: new `"compact"` mode renders
  `outputs <- instruction(inputs)` and drops static params; `"full"` (default)
  is unchanged.
- `detail` threaded through `RecipeBook.describe/to_mermaid/to_dot/visualize`
  and the top-level `api.describe` / `api.visualize`. Default stays `"full"`, so
  no existing output changes unless callers opt in.
- New `RecipeBook._render_dependency_names()` expands wildcard deps to the
  produced targets they match; both renderers use it, so edges mirror the real
  DAG. `_build_graph_view` no longer classifies wildcard patterns as external
  inputs, removing the stray `part_*` pseudo-node.

**Verify.** `32 passed`. ccd ingest `detail="compact"` collapses the column
dump to `_chem_comp_atom_dict <- _worker(_chem_comp_atom)`. Wildcard collector
renders `target_part_a -> step_3` / `target_part_b -> step_3` (mermaid + dot)
with no raw `part_*` node.

**Note.** Compact labels show the instruction's `__name__`; factory-built
instructions surface their closure name (`_worker`). Cosmetic, and a
StructCooker-side `functools.wraps` concern — out of scope here.

## Step 7 (P2) — Decompose Reader/Writer hooks

**The smell.** `ReaderHooks` was a flat 4-field bag conflating three concerns
(file `loader`; payload `deserializer`+`adapter`; `key_transform`). Worse, every
boundary that accepted an optional `reader` re-implemented the same field-by-field
"explicit-wins-over-defaults" merge by hand (`decode_payload`, `load_inputs`,
`build_lmdb`, `rebuild_lmdb`, `extract_lmdb_records`) — six near-identical blocks
that were easy to get subtly wrong.

**How (behavior-preserving).**
- Introduced cohesive boundary types: `InputHooks(loader)` and
  `PayloadHooks(deserializer, adapter)`, exported from the package.
- `ReaderHooks` keeps its flat fields (so the many call sites and the YAML/CLI
  config builders keep working) but now exposes the concerns as `.input` /
  `.payload` views.
- Added `ReaderHooks.merge()` / `WriterHooks.merge()` (every non-None field of
  the override wins). All six hand-rolled merge blocks were replaced by a single
  `defaults.merge(override)` call.

**Scope decision (honest).** I deliberately did **not** rip `ReaderHooks` out of
the LMDB / runner / CLI signatures in favor of the focused types. Those code
paths (`build_lmdb` / `rebuild_lmdb` / `extract_lmdb_records`) cannot be executed
in this environment — `lmdb`/`joblib` are not installed, so their tests *skip*.
A wide breaking rewrite of code I can't run is irresponsible; instead I made a
behavior-preserving consolidation (verified by analysis: `merge()` reproduces the
original field precedence exactly) and introduced the cohesive types as the
migration target. **Follow-up:** complete the physical migration to
`InputHooks`/`PayloadHooks` in the LMDB internals once a runnable lmdb test
environment exists.

**Verify.** `32 passed` (runnable subset, incl. `decode_payload`/`encode_output`
and runner-hook tests). Ad-hoc: `.input`/`.payload` views, `merge()` precedence
(override non-None wins, `None` falls through, `merge(None) is self`), and
`decode_payload` with deserializer+adapter all correct.

## Step 8 (P2) — Redesign `with_namespace`

**The smell.** `with_namespace` rebuilt every step by hand-constructing a fresh
`Recipe(targets=..., instruction=..., inputs=..., name=..., namespace=...,
tags=..., metadata=...)`. Any field later added to `Recipe` would be silently
dropped unless someone remembered to thread it here — classic field drift.

**Why not "no rewrite at all".** A truly copy-free design (store the namespace
on the book, prefix names lazily at lookup) was evaluated and rejected: the
executor resolves dependencies *by name*, so a lazy view would have to
re-synthesize a name-rewritten `Recipe` on **every** `__getitem__` — uncached and
strictly slower than building the namespaced book once, plus it would entangle
namespace logic into the executor. A namespaced book is a legitimately distinct
graph; building it once (and letting step 4's memoization cache its analysis) is
the right model. The real defect was the *fragility* of the rewrite, not the
rewrite itself.

**How.** `with_namespace` now delegates to a `_namespace_recipe()` helper that
uses `dataclasses.replace`, overriding only the namespace-sensitive fields
(`targets`, `inputs`, `namespace`). `name`, `tags`, `metadata`, and any future
field are copied automatically — drift is impossible.

**Verify.** `32 passed`. Ad-hoc: `rb.with_namespace("mod")` prefixes targets to
`mod.a`/`mod.b`, rewrites the internal dep `a -> mod.a` while leaving the external
input `seed` bare, **preserves** `name`/`tags`/`metadata`/`namespace` on the
step, and executes end-to-end (`{"mod.b": "HI!"}`).

---

## Step 9 (hardening) — Runnable LMDB env + self-reference cycle bug

**Context.** To actually verify the step-7 hooks on real LMDB paths, the optional
deps (`lmdb`, `joblib`, `omegaconf`) were installed into the pixi env. This
un-skipped the LMDB/config tests for the first time and surfaced **3 pre-existing
failures** that the skips had always masked. Confirmed pre-existing by
`git stash`-ing all working-tree changes and reproducing them at HEAD — not
introduced by this refactor.

**Real engine bug found & fixed — false self-cycle.** A legitimate pass-through
step whose output shares a name with one of its inputs (`X <- X`, e.g.
`value <- value`) was reported by the static analyzer as
`Cycle detected: value -> value`, although the executor handles it at runtime
(the provided input shadows the target). Root cause: `_expand_dependency` turned
a self-referential dependency into a self-edge. Fix: a dependency whose name is
one of the step's own targets is treated as an **external input** (no edge),
matching runtime semantics; `required_inputs` updated likewise so such a step
reports its input as required. Genuine 2-cycles (`a<->b`) still detected.

**One broken test corrected.** `test_build_and_read_lmdb` reused a single
`double <- value` recipe for *both* build and extract, but
`extract_lmdb_records` injects records as `{"key", "db_data"}` (the contract the
config test's extract recipe already follows). The test could never have passed;
extract was given its own `double <- db_data` recipe. Correction of a never-run
test, not new coverage.

**Result:** `41 passed, 0 skipped` with all optional deps present (was
`32 passed, 9 skipped`). The LMDB/config code paths are now actually exercised
and the step-7 `merge()` changes are verified end-to-end on real build / rebuild
/ extract.

**Env note.** `lmdb`/`joblib`/`omegaconf` were pip-installed into the pixi env
(already declared under `project.optional-dependencies`). For a clean lockfile
prefer `pixi add` or the `[lmdb]`/`[parallel]`/`[config]` extras.

---

## Step 10 (hardening) — Concurrency contract

**Why.** Nothing stated whether the engine was thread-safe; `Cooker` /
`ExecutionContext` carry mutable per-run state. A tier-1 library states this.

**How.** Documented the contract on `ExecutionContext` and `Cooker` docstrings:
a single instance is **not** safe for concurrent use (one per run, one thread);
the `RecipeBook` is read-only during execution and may be shared; parallelism is
achieved per-work-item (the pattern `datacooker.processing` and the LMDB helpers
already use — `execute()` builds a fresh context+Cooker per call). No code
change; this codifies the existing design.

## Step 12 (hardening) — Warn on wildcard 0-match

**Why.** A wildcard arg that matched nothing (e.g. a typo'd pattern) silently
contributed zero values — a quiet footgun.

**How.** The executor now logs a `warning` when a wildcard argument resolves to
no inputs/targets, naming the pattern and the producing step. Verified: a
`nomatch_*` arg yields `{"items": []}` *and* emits the warning.

## Step 11 (hardening) — Remove legacy compatibility cruft / dual-naming

**Decision.** Back-compat dropped (this is v3). One canonical hook vocabulary
survives: **`loader`, `adapter`, `deserializer`, `serializer`, `key_transform`,
`materializer`, `key_builder`**.

**Removed (every site):**
- `from_mapping` legacy kwargs and dual YAML keys on both hook classes.
- Dual callable params on `read_lmdb`, `build_lmdb`, `rebuild_lmdb`,
  `extract_lmdb_records`, `parse_file`, `execute`, `parse_dict`,
  `parallel_process(_report)` — the legacy halves (`load_func`, `convert_func`,
  `deserialize`, `serialize`, `transform_func`, `project_func`, `key_func`)
  are gone.
- `resolve_key_transform`'s `transform_func` arg; `_shared.rebuild_entry`'s
  internal `transform_func` (→ `key_transform`).
- The entire `datacooker.orchestration` module
  (`run_workflow` / `run_parallel_workflow` / `extract_lmdb_workflow`) — a pure
  legacy-naming wrapper layer nothing depended on. Canonical runners
  (`run_recipe` / `run_recipe_batch` / `run_lmdb_extract`) remain.
- Config `DEFAULT_CALLABLE_KEYS` and the CLI `_pop_*_hooks` now read canonical
  keys only.

**Synced (forced by the API change, not new work):** `runners/core` (was still
passing `transform_func` to `execute`), `tests/test_utils` + `tests/test_config`
(legacy kwargs → canonical), README API section, and StructCooker's one call
site (`read_lmdb(deserialize=...)` → `deserializer=`). `parse_file`'s 3rd
positional is now `loader` (positional callers unaffected).

**Verify.** `ruff check src tests` clean (also caught/closed unused imports from
the dropped params); `41 passed, 0 skipped`; package import shows the wrappers
gone (`__all__` 53 → 50) and canonical names present; StructCooker ccd +
`read_lmdb` consumer paths still work.

## Summary

Planned steps 1–8 + hardening steps 9, 10, 12 landed; step 11 is deferred as a
back-compat **decision**. With all optional deps installed the suite is
**`41 passed, 0 skipped`** and `ruff check src tests` is clean (baseline was
`32 passed, 9 skipped`, with the LMDB/config paths unverified and a real
self-cycle bug latent).

**Release:** **v3.0** — P2 changes the public surface (new `InputHooks` /
`PayloadHooks` exports; `match_targets` / `adjacency` now public; behavior of
wildcard ordering changed). P0+P1 on their own would have been v2.2.

**Behavior changes downstream code should know about:**
- Wildcard dependencies now create real graph edges and resolve deterministically
  (previously order-dependent / possibly skipped).
- `is_supported_annotation` is stricter: malformed annotations are rejected
  instead of silently accepted.
- Mermaid step labels use `<br/>` (was a literal `\n`).
- `describe`/`visualize` gained a `detail="compact"` option (default `"full"`
  unchanged).

**Carried-forward follow-ups (not done here, by design):**
- Finish migrating the LMDB internals to `InputHooks`/`PayloadHooks` once a
  runnable `lmdb` test environment exists (step 7).
- Optional: `functools.wraps` on StructCooker instruction factories so compact
  labels show a meaningful instruction name instead of `_worker` (step 3).
- Out of scope this pass (owner deferred): expanded tests, StructCooker consumer
  cleanup, user-facing docs.
</content>
