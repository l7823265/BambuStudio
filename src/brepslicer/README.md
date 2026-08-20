# BrepSlicer integration (BambuStudio)

Parametric B-rep STEP slicing wired into BambuStudio at slice time.

## Branch

```text
feature/brepslicer-step-integration
```

https://github.com/l7823265/BambuStudio/tree/feature/brepslicer-step-integration

## Build (for others)

1. Clone this fork and check out the branch above.
2. Build BambuStudio the normal way (deps + RelWithDebInfo / Release).
3. OpenCASCADE is already part of BambuStudio `deps/OCCT` — no extra install.
4. BrepSlicer **core sources are vendored** under `src/brepslicer/core`.  
   CMake picks them up automatically; optional override:

```text
-DBREPSLICER_CORE_DIR=/path/to/external/BrepSlicer-core
```

5. After build, use menu **导入STEP** / Import STEP, then slice.  
   Small STEP files use B-rep contours for G-code; oversized architectural CAD falls back to mesh.

## Layout

| Path | Role |
|------|------|
| `src/brepslicer/core/` | Independent BrepSlicer L2 sources |
| `src/libslic3r/BrepSlice.*` | Studio adapter (transform, units, ExPolygons) |
| `src/brepslicer/CMakeLists.txt` | Static `brepslicer` lib linked into `libslic3r` |
