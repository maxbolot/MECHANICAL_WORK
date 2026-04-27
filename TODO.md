# Project TODO

This document tracks pending architecture and workflow improvements. Items below are intentionally scoped as investigation-first tasks; implementation decisions will be made later.

## Scope note from README

- MATLAB remains the authoritative analysis stack for this project.
- Notebooks stay a dissemination layer; there is no planned migration away from MATLAB workflows.

## 1) Histogram clamping options in precipitation binning

- Add configurable clamping behavior in the histogram subroutine for precipitation values that fall outside bin edges.
- Support the following modes:
  - No clamping
  - Clamp to lower edge only
  - Clamp to upper edge only
  - Clamp to both lower and upper edges
- Define a clear interface for selecting mode (e.g., enum/integer flag or named option).
- Investigate impact on diagnostics reproducibility and interpretation.
- Decision pending: default mode and backward-compatibility behavior.

## 2) Fork precipitation-rate threshold computation (grid-cell-by-grid-cell)

- Create a fork/variant of the precipitation-rate threshold computation where thresholds are computed independently at each grid cell.
- Compare against current approach:
  - Computational cost
  - Memory footprint
  - Numerical behavior and scientific interpretation
- Identify required I/O/schema changes for storing per-grid-cell thresholds.
- Decision pending: whether to keep both approaches or replace current workflow.

## 3) Decouple histogram logic from work/lift grid computation

- Decision: use a modular single-kernel design as the primary path.
- Keep work/lift grid computation and histogram accumulation callable independently via runtime mode switches:
  - `gridded_only`
  - `hist_only`
  - `gridded_plus_hist`
- Support multiple histogram regions in one run (regional boxes, latitude-band partitions, optional mask-based regions).
- Keep one optional standalone histogram driver only as a thin wrapper over the same shared module APIs.

### 3a) Implementation phases (Fortran)

- Phase 1: Refactor current histogram accumulation/writing into module procedures with no output changes.
- Phase 2: Add region-list configuration and loop over region sinks during accumulation.
- Phase 3: Add runtime mode switches for gridded output vs histogram output vs both.
- Phase 4: Add optional standalone histogram executable reusing the same modules.

### 3b) Build/layout decision

- Compile histogram module procedures first as an internal static library under `local/` scope (e.g., `local/lib` or archive in `local/obj`).
- Defer extraction to project-root/shared library until API stability and reuse requirements are validated.

## 4) Roadmap integration: lib extraction

- Shape `matlab/lib` as a stable shared API layer while keeping wrappers/runners repository-local.
- Prepare a versioned reusable MATLAB utility package for:
  - NetCDF readers/validators
  - Time-axis normalization/conversion
  - Time-weighting/statistics helpers
- Identify and remove duplicate helper logic across analyses before extraction.
- Define extraction boundaries and ownership:
  - What stays local to this repository (`matlab/entries`, `matlab/analysis`)
  - What is published in the shared library
- Plan migration path for dependent workflows/publication repos with minimal disruption.
- Decision pending: release/versioning strategy and extraction timeline.
