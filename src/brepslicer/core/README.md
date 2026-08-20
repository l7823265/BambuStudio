# BrepSlicer

Parametric slicing on B-rep (no mesh). L0 + L1 + L2: STEP in, typed closed contours out
(line / arc / ellipse / cubic B-spline).

## Build (MSVC / VS 2019)

OCCT 7.5 is the same install used by `D:\OCC`:

`D:\lib\opencascade-7.5.0\build_vs2019\win64` (Debug `vc14/libd` + `bind`).

```bat
cd D:\project\BrepSlicer\core
build.bat
```

Or:

```bat
cmake --preset msvc-2019-x64
cmake --build --preset debug --parallel
```

Binaries: `_bin\x64_Debug\brep_slice.exe`, `slice_verify.exe`, `slice_selftest.exe`.

## CLI

```bat
brep_slice model.step --normal 0,0,1 --layer-height 0.05 --out slices.json
brep_slice model.step --normal 0,0,1 --layer-height 0.05 --svg svg_out --dxf slices.dxf --out slices.json
slice_verify model.step --normal 0,0,1 --layer-height 0.05
slice_selftest
```

SVG uses path `A` for circular arcs. DXF writes `LINE` / `ARC` / `CIRCLE` / `ELLIPSE` (not polylines), one CAD layer per slice (`SLICE_0000`, …). Outer contours are blue, holes red.

## Layout

See `DESIGN.md` for tolerances, loop closing, and degenerate-case decisions.
