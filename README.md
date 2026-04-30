# Bambu Arc Overhang

Post-process a Bambu Studio `.3mf` plate so the bridge layers below top surfaces become arc-overhang spirals — the trick where filament wraps around itself in concentric arcs to support an unsupported area without support material.

This repo is a Bambu-focused fork of [Wasupmacuz/arc-overhang-prusaslicer-integration](https://github.com/Wasupmacuz/arc-overhang-prusaslicer-integration), which itself builds on Steven McCulloch's original [arc-overhang algorithm](https://github.com/stmcculloch/arc-overhang) and Nic's PrusaSlicer integration. The upstream targets PrusaSlicer / OrcaSlicer; this fork adds:

- A Bambu Studio `.3mf` wrapper ([`bambu_arc_overhang.py`](bambu_arc_overhang.py)) that extracts plate gcode, runs the post-processor, refreshes the md5 sidecar Bambu printers verify, and repacks the archive.
- A "supports a top surface" filter so only bridges that actually hold up a top solid surface become arcs (no internal-web bridges that span over voids).
- A coverage filter (`--min-top-coverage`, default 50%) that skips bridges whose immediate next layer is mostly *not* solid above them.
- New defaults tuned for fast prints (5 mm/s arcs, no above-arc cooling slowdown, no Hilbert cross-hatch).
- A coverage diagnostic that surfaces total unsupported surface area per plate so you can see when arcs failed to fill a bridge.
- Multi-plate `.3mf` support with parallel processing.
- A scipy-cKDTree fast path for the post-processor's main perf hotspot — single plate runtime down ~4.4× from baseline.

## Setup

Requires Python 3.12+ (the upstream uses nested-quote f-strings).

```sh
python3.12 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
pip install scipy        # optional but strongly recommended; ~10× faster distance queries
```

## Usage

```sh
.venv/bin/python bambu_arc_overhang.py <plate.3mf>
```

Output is written next to the input as `<name>_arc.3mf`. Multi-plate `.3mf` files (Bambu Studio "Send all to printer" / "Export all plates") are auto-detected and each plate is processed in parallel.

```sh
# basic
.venv/bin/python bambu_arc_overhang.py /path/to/plate_10.gcode.3mf

# specify output
.venv/bin/python bambu_arc_overhang.py /path/to/plate.3mf -o /path/to/out.3mf

# tune for a specific part
.venv/bin/python bambu_arc_overhang.py /path/to/plate.3mf \
    --arc-speed 8 --min-top-coverage 0.7 --workers 4
```

### What it converts

For each layer's `; FEATURE: Bridge` regions, the wrapper qualifies a bridge for arc replacement only if:

1. There's an unbroken stack of `Internal solid infill` (or `Top surface`) directly above the bridge XY, capped by a `Top surface` within `--bridge-look-ahead` mm above. (Filters out bridges that aren't supporting anything.)
2. At least `--min-top-coverage` (default 50%) of the bridge's area is directly under solid material in the *immediately next* layer. (Filters out internal-web bridges that span over sparse infill.)

Bridges that don't qualify are left as-is; bridges that qualify get their bridge gcode deleted and arc gcode injected in its place.

### CLI knobs

```
--arc-speed FLOAT             Arc print speed mm/s. Default 5.
--arc-min-speed FLOAT         Floor speed for small slowdown arcs mm/s. Default 2.
--min-top-coverage FLOAT      Min fraction of bridge area under solid in next
                              layer (0..1). Default 0.5. Set 0 to disable.
--arc-center-offset FLOAT     Override ArcCenterOffset (mm). Lower values let
                              the BFS reach into thin sections; higher values
                              improve arc bonding. Unset = upstream default
                              (1.5 × nozzle_diameter).
--min-fill-ratio FLOAT        If a bridge's BFS fills less than this fraction
                              of its area, reject the bridge entirely — keep
                              the bridge gcode, no arcs. Use this if the
                              coverage diagnostic shows unfilled area on a
                              load-bearing layer. Try 0.85 or 0.95. Default
                              0 = always replace.
--above-arcs-zdist FLOAT      Re-enable above-arc cooling pass within this
                              vertical distance (mm). Default 0 = OFF.
--enable-hilbert-cooling      Restore the Hilbert-curve cross-hatch fill above
                              arcs. Requires --above-arcs-zdist > 0.
--fan-boost-whole-layer       If above-arcs cooling is on, force fan boost
                              across whole layers near arcs.
--workers INT                 Max plates to process in parallel. Default = CPU
                              count.
-o OUTPUT                     Output .3mf path. Default: <input>_arc.3mf.
```

### Coverage diagnostic

Each plate's run prints how much of the bridge area got arc support and which layers have the largest gaps:

```
[ 11.1s] Metadata/plate_10.gcode: modified
    Coverage: 17541 mm^2 of bridges converted, 507 mm^2 unsupported (2.9%)
      Layers with unfilled area > 0.5 mm^2 (worst first):
        layer 129: 389.7 mm^2 unfilled (9% of 4526 mm^2 bridge)
        layer 54:  100.7 mm^2 unfilled (2% of 5228 mm^2 bridge)
```

"Unsupported" means the bridge gcode was deleted but the arc BFS couldn't fill that area — there's no support there for the layer above. Causes are usually thin sections narrower than `2 × ArcCenterOffset` (try lowering `--arc-center-offset`) or unusual geometry the BFS skips (preview the layer in Bambu Studio).

If you see significant unsupported area on a part you care about, options are:

- **Use `--min-fill-ratio 0.95`** to refuse to arc-replace any bridge that the BFS can't ~fully cover. The bridge gcode is kept intact (printed as bridge), so the layer above is still supported — you just lose arcs on that bridge. Safest mode; use for load-bearing parts.
- Lower `--arc-center-offset` to let arcs reach into thinner sections (try 0.2; not always better — measure with the metric).
- Raise `--min-top-coverage` (e.g. 0.7) to skip bridges that would have low coverage anyway.
- Slice the part at a different orientation so bridges become more uniform.

### Print settings

The arc-overhang technique relies on each filament strand cooling before the next is laid on top of it. Recommended settings in Bambu Studio for parts where this matters:

- Bridge speed ≤ 5 mm/s, or use `--arc-speed 3` to slow only the arcs.
- Part fan at 100%.
- Lower nozzle temp than usual (PLA at ~190 °C works well).
- The script disables Bambu Studio's slicer-side bridging slowdown above arcs by default; if you want it back, use `--above-arcs-zdist 3`.

## Architecture

```
bambu_arc_overhang.py    # Bambu .3mf wrapper, multi-plate, marker translation, md5 refresh
arc_overhangs_v1.0.0.py  # The post-processor (forked from upstream).
                         # Modifications:
                         #  - OnlyBridgesSupportingTopSurfaces filter
                         #  - MinTopSurfaceCoverageRatio filter
                         #  - Inner-wall fallback in makeStartLineString
                         #  - cKDTree-based distance for the BFS hotspot
                         #  - Lazy pre-pass for chain-check polys
                         #  - Coverage diagnostic
```

## Why this is a fork

The upstream post-processor only knows about PrusaSlicer / OrcaSlicer marker syntax (`;TYPE:Bridge`, `;LAYER_CHANGE`). Bambu Studio uses `; FEATURE: Bridge` / `; CHANGE_LAYER`. Rather than patching every detection site, this fork's wrapper translates Bambu markers to Orca on the way in and back to Bambu on the way out, so the post-processor sees what it expects.

The post-processor itself was also extended for two needs upstream didn't address:

- **Internal bridges that don't touch the outer perimeter.** The original arc generator anchors the first arc on the previous layer's outer wall. Internal bridges span over voids in the model interior and never touch that wall, so arc generation failed. We added a fallback in `makeStartLineString` that also considers inner-wall rings — fixes layer-129-style geometry.
- **"Internal web" bridges over sparse infill.** Bambu marks these as `; FEATURE: Bridge` even though they're not really supporting a top surface. The new top-surface chain check + coverage ratio filter both narrow the scope to bridges that genuinely hold up solid material.

## Limitations

- Layer 129 of the test plate still has ~9% unfilled area at default settings — thin elongated bridge sections (< 1.7 mm wide) that can't fit the BFS arc geometry. Use the coverage diagnostic to spot this kind of region before printing.
- Single-plate runtime on a 12 MB gcode is ~11 s; Bambu Studio's slicer is ~4 s on the same input. We're 2.5–3× behind their C++ pipeline; remaining headroom is mostly inside shapely's per-call overhead which is hard to chip away at without a deeper rewrite.
- Print correctness comes from preview, not just from the script running clean. Always inspect the modified `.3mf` in Bambu Studio's gcode preview before sending to the printer.

## License

Inherits the upstream's license — see [LICENSE](LICENSE).
