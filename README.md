# Bambu Arc Overhang

Post-process a Bambu Studio `.3mf` plate so the bridge layers below top surfaces become arc-overhang spirals — the trick where filament wraps around itself in concentric arcs to support an unsupported area without support material.

This repo is a Bambu-focused fork of [Wasupmacuz/arc-overhang-prusaslicer-integration](https://github.com/Wasupmacuz/arc-overhang-prusaslicer-integration), which itself builds on Steven McCulloch's original [arc-overhang algorithm](https://github.com/stmcculloch/arc-overhang) and Nic's PrusaSlicer integration. The upstream targets PrusaSlicer / OrcaSlicer; this fork adds:

- A Bambu Studio `.3mf` wrapper ([`bambu_arc_overhang.py`](bambu_arc_overhang.py)) that extracts plate gcode, runs the post-processor, refreshes the md5 sidecar Bambu printers verify, and repacks the archive.
- Multi-plate `.3mf` support — every `Metadata/plate_*.gcode` in the archive is processed in parallel.
- A "supports a top surface" filter so only bridges that actually hold up a top solid surface become arcs (no internal-web bridges that span over voids).
- A coverage filter (`--min-top-coverage`, default 50%) that skips bridges whose immediate next layer is mostly *not* solid above them.
- A closing operation on bridge polygons (`--bridge-closing`, default 1.0 mm) so the arc BFS sees the bridge's full surface — independent of the slicer's chosen infill direction — instead of a comb of thin parallel polys with gaps.
- **Preserved bridge gcode in unfilled regions.** Wherever the arc BFS can't reach inside a converted bridge, the script re-emits the original bridge gcode at `bridge_speed` so the layer above is never left without support. This is on by default; removing supporting bridge is never desired.
- A coverage diagnostic that splits each bridge into arc-filled / preserved-as-bridge / truly-unsupported buckets so problems are easy to spot.
- Arcs inherit the slicer's bridge profile by default — `bridge_speed`, `travel_speed`, and `default_acceleration` from the Bambu CONFIG_BLOCK control arc print speed, travel feedrate, and M204 emitted before each arc block. Tune in Bambu Studio; CLI overrides remain available.
- Fan is always forced to 100% on the arc layer (overrides slicer profile) — arcs need maximum cooling regardless of how the rest of the print is tuned.
- Other defaults tuned for fast prints (no above-arc cooling slowdown, no Hilbert cross-hatch).
- A scipy-cKDTree fast path for the BFS distance hotspot — single plate runtime down ~4.4× from baseline. A 10-plate `.3mf` (~660 MB unzipped) processes in ~60 seconds wall-time.

## Setup

Requires Python 3.12+ (the upstream uses nested-quote f-strings).

```sh
python3.12 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

`requirements.txt` includes scipy, which the post-processor uses for fast point-to-polyline distance via `cKDTree`. Without it the script falls back to a slower segment-loop, but install scipy unless you have a reason not to.

## Usage

```sh
.venv/bin/python bambu_arc_overhang.py <plate.3mf>
```

For an input named `<name>.gcode.3mf` the script writes:

- `<name>_arc.gcode.3mf` — multi-plate, for previewing all plates in one Bambu Studio session.
- `<name>_arc_plate_N.gcode.3mf` (one per plate) — single-plate, for **uploading** to a printer.

Use the per-plate file when sending to a printer. Bambu Studio's `Send to printer` short-circuits and uploads the entire loaded `.3mf` whenever it sees a gcode-only archive (`m_exported_file == true` in [`Plater::send_gcode`](https://github.com/bambulab/BambuStudio/blob/master/src/slic3r/GUI/Plater.cpp)), so a multi-plate archive uploads ALL plates regardless of which one you select. The single-plate files are renumbered to `plate_1` and have their `model_settings.config` / `slice_info.config` / `_rels/model_settings.config.rels` trimmed to one plate, so each upload is just that plate's gcode.

Pass `--no-multi-plate` or `--no-per-plate` to skip either side. Multi-plate `.3mf` files (Bambu Studio "Send all to printer" / "Export all plates") are auto-detected and each plate is processed in parallel.

### Auto-run as a Bambu Studio post-processing script

Bambu Studio can run external scripts on each plate's G-code at slice time. The script is invoked once per plate, the slicer's temp `.gcode` path is appended as the last argument, and the script must **modify the file in place** (see `libslic3r/GCode/PostProcessor.cpp::run_post_process_scripts`). When the input is a raw `.gcode` file (not a `.3mf`), `bambu_arc_overhang.py` enters in-place mode automatically.

To enable:

1. In Bambu Studio, open the **Process** preset → **Others** tab → **Post-processing scripts**.
2. Paste this into the textarea (one script per line):

   ```
   /absolute/path/to/repo/.venv/bin/python /absolute/path/to/repo/bambu_arc_overhang.py
   ```

   Add CLI flags after the script path the same way you would on the command line, e.g. `--arc-speed 5 --min-top-coverage 0.7`. **Don't quote the gcode path** — Bambu appends it itself.

3. Save the preset (the modified-preset indicator appears next to the preset name).

Caveats:

- Slicing time goes up by however long `bambu_arc_overhang.py` would take on that plate's gcode (~3–20 s per plate, depending on bridge complexity).
- The setting lives on the active Process preset. Add it to every Process preset you slice with.
- After post-processing, Bambu re-parses the gcode for the preview. There's an upstream bug where this re-parse path breaks the H2C/H2D nozzle auto-mapping table (`The printer failed to build the nozzle auto-mapping table`); X1/P1 are unaffected. See the comment block in `BackgroundSlicingProcess.cpp::finalize_gcode`.
- Auto-running the post-processor produces a multi-plate `.gcode.3mf` on export, which still uploads the entire archive on Send-to-Printer. Use the wrapper's `.3mf` mode (see above) afterward to emit per-plate uploadable files, or skip post-processing in Bambu and process the exported `.gcode.3mf` end-to-end with this script.
- The first slice after a script is configured shows a security-warning modal listing the script text. Click **Yes** to run, **No** to skip post-processing for that session, or close the dialog to cancel the slice. The choice is remembered until the project is closed.

```sh
# basic
.venv/bin/python bambu_arc_overhang.py /path/to/plate.3mf

# specify output
.venv/bin/python bambu_arc_overhang.py /path/to/plate.3mf -o /path/to/out.3mf

# tune for a specific part
.venv/bin/python bambu_arc_overhang.py /path/to/plate.3mf \
    --arc-speed 8 --min-top-coverage 0.7 --workers 4
```

## What it converts

For each layer's `; FEATURE: Bridge` regions, the wrapper qualifies a bridge for arc replacement only if:

1. There's an unbroken stack of `Internal solid infill` (or `Top surface`) directly above the bridge XY, capped by a `Top surface` within ~4 mm above. (Filters out bridges that aren't supporting anything.)
2. At least `--min-top-coverage` (default 50%) of the bridge's area is directly under solid material in the *immediately next* layer. (Filters out internal-web bridges that span over sparse infill.)

Bridges that don't qualify are left as-is. For bridges that qualify:

1. The bridge polygon is built by buffering the slicer's gcode path and then applying a closing operation (`--bridge-closing`, default 1.0 mm) so adjacent gcode lines merge into one solid polygon covering the bridge's actual surface.
2. The arc BFS fills as much of that polygon as the algorithm can reach.
3. **Wherever the BFS leaves a gap, the original bridge gcode is re-emitted in that gap region** (`PreserveBridgeInUnfilled`, default on). The layer above always has either an arc or the original bridge underneath it.
4. Bridge gcode in arc-filled regions is deleted; bridge gcode outside the polygon footprint is untouched.

## CLI knobs

```
--arc-speed FLOAT             Arc print speed mm/s. Default = the slicer's
                              bridge_speed (from CONFIG_BLOCK). Tune the bridge
                              profile in Bambu Studio and the arcs follow.
                              Pass an explicit value to override.
--arc-min-speed FLOAT         Floor speed for short arcs that the BFS slows
                              down. Default = --arc-speed (no slowdown — arcs
                              print at bridge_speed regardless of length).
                              Set explicitly to enable per-arc slowdown for
                              small radii.
--min-top-coverage FLOAT      Min fraction of bridge area under solid in next
                              layer (0..1). Default 0.5. Set 0 to disable.
--bridge-closing FLOAT        Closing-operation radius (mm) applied to bridge
                              polygons. Merges adjacent buffered gcode lines
                              so the BFS sees the bridge's full surface
                              regardless of the slicer's infill direction.
                              Default 1.0 mm. Set 0 to disable.
--arc-center-offset FLOAT     Override ArcCenterOffset (mm). Lower values let
                              the BFS reach into thinner sections; higher
                              values improve arc bonding on small radii.
                              Unset = upstream default (1.5 × nozzle_diameter).
--min-fill-ratio FLOAT        If a bridge's BFS fills less than this fraction
                              of its area, reject the bridge entirely — the
                              whole bridge stays as bridge gcode, no arcs.
                              Use this if you'd rather a load-bearing bridge
                              skip arcs entirely than have them mixed with
                              preserved bridge segments. Default 0 (off).
--above-arcs-zdist FLOAT      Re-enable above-arc cooling pass within this
                              vertical distance (mm). Default 0 = OFF. Upstream
                              default was 3 mm, which slowed perimeters and
                              infill on ~15 layers above each arc.
--enable-hilbert-cooling      Restore the Hilbert-curve cross-hatch fill that
                              overwrites internal solid infill above arcs.
                              Requires --above-arcs-zdist > 0. Off by default.
--fan-boost-whole-layer       If above-arcs cooling is on, force fan boost
                              across whole layers near arcs.
--workers INT                 Max plates to process in parallel. Default =
                              one worker per plate. Pass an integer to cap it.
--no-multi-plate              Skip writing the multi-plate <name>_arc.gcode.3mf.
--no-per-plate                Skip writing the per-plate
                              <name>_arc_plate_N.gcode.3mf files. Bambu Studio
                              uploads the entire .3mf when you Send-to-Printer
                              from a multi-plate gcode 3MF; the per-plate files
                              exist so each upload is one plate.
-o OUTPUT                     Multi-plate output path. Default:
                              <input-name-with-extension>_arc.gcode.3mf.
                              Per-plate files derive from this by inserting
                              _plate_N before the extension.
```

## Coverage diagnostic

Each plate's run prints the bridge area split into three buckets:

```
[ 12.6s] Metadata/plate_2.gcode: modified
    Coverage: 18190 mm^2 of bridges -> 16025 arc-filled (88.1%),
              2080 preserved as bridge (11.4%),
              85 truly unsupported (0.5%)
      Layers with truly-unsupported area > 0.5 mm^2 (worst first):
        layer 128: 85.1 mm^2 unsupported (2% of 4561 mm^2 bridge)
```

| Bucket | Meaning |
|---|---|
| **arc-filled** | Bridge area replaced by arcs. |
| **preserved as bridge** | Original bridge gcode re-emitted because arcs couldn't reach. Still supports the layer above; just at bridge speed instead of arc speed. |
| **truly unsupported** | Bridge area covered by *neither* arcs *nor* preserved bridge gcode. This is the only number that means a real gap. With default settings it's near zero — only the closing-radius extension (the 1 mm strip beyond where the original bridge gcode actually was) can fall in this bucket. |

If "truly unsupported" is non-zero on a load-bearing layer, options are:

- Lower `--bridge-closing` (try 0.5) so the bridge polygon stays closer to the original gcode footprint and the BFS doesn't need to reach into the closing extension.
- Use `--min-fill-ratio 0.95` to refuse to arc-replace any bridge the BFS can't ~fully cover. The whole bridge stays as bridge — you lose arcs on it, but support is intact.
- Lower `--arc-center-offset` (try 0.2) to let arcs reach into thinner sections. Sometimes worse — measure with the diagnostic.
- Slice the part at a different orientation so bridges become more uniform.

## Print settings

The arc-overhang technique relies on each filament strand cooling before the next is laid on top of it. Recommended settings in Bambu Studio for parts where this matters:

- Part fan is forced to 100% on the layer with arcs (the script always emits `M106 S255` regardless of slicer profile).
- Lower nozzle temp than usual (PLA at ~190 °C works well).
- The script disables Bambu Studio's slicer-side bridging slowdown above arcs by default. If you want it back, use `--above-arcs-zdist 3`.
- **Bridge speed in your Bambu profile drives both the *preserved* bridge segments and the *arcs themselves* by default.** Set it ≤ 10 mm/s if your part has thin bridge regions. Pass `--arc-speed N` to override only the arc speed without touching the slicer profile.

## Performance

| Workload | Wall time |
|---|---|
| Single plate, 12 MB gcode, 6 bridge layers | ~11 s |
| 10-plate `.3mf`, ~660 MB unzipped, parallel | ~26 s |
| 8-plate `.3mf`, ~410 MB unzipped, parallel | ~22 s |

The biggest hotspot was per-point shapely.distance during the BFS — replaced with a cKDTree query against a densified copy of the polygon boundary, which gave a ~4.4× speedup on the test plate. Bambu Studio's slicer can re-slice the same input in ~4 s; we're ~2.5–3× behind their C++ pipeline. Most of the remaining time is inside shapely's `intersection` / `intersects` / `buffer` calls, which are inherent to the BFS algorithm.

The repack step uses zlib level 1 (vs. the Python default 6) and writes the multi-plate and per-plate files in parallel via `ThreadPoolExecutor`. On the 10-plate test that dropped repack time from ~18 s to ~5 s; output `.3mf` is ~10–15% larger but transfer-time-equivalent since printers read the file once and decompression is cheap. Wall time is bounded by the slowest single plate (parallelism stops at one core per plate); an intra-plate split would be possible but isn't currently implemented.

## Architecture

```
bambu_arc_overhang.py    # .3mf wrapper: multi-plate processing, marker
                         # translation, md5 refresh, parallel repack, and
                         # per-plate single-plate .gcode.3mf emission so
                         # Send-to-Printer uploads one plate at a time.
arc_overhangs_v1.0.0.py  # post-processor (forked from upstream).
                         # Modifications:
                         #  - OnlyBridgesSupportingTopSurfaces filter
                         #  - MinTopSurfaceCoverageRatio filter
                         #  - BridgePolyClosingRadius (closing op on bridge polys)
                         #  - PreserveBridgeInUnfilled (re-emit bridge gcode in BFS gaps)
                         #  - Inner-wall fallback in makeStartLineString
                         #  - cKDTree-based distance for the BFS hotspot
                         #  - Lazy pre-pass for chain-check polys
                         #  - 3-bucket coverage diagnostic
```

## Why this is a fork

The upstream post-processor only knows about PrusaSlicer / OrcaSlicer marker syntax (`;TYPE:Bridge`, `;LAYER_CHANGE`). Bambu Studio uses `; FEATURE: Bridge` / `; CHANGE_LAYER`. Rather than patching every detection site, this fork's wrapper translates Bambu markers to Orca on the way in and back to Bambu on the way out, so the post-processor sees what it expects.

The post-processor itself was also extended for needs upstream didn't address:

- **Internal bridges that don't touch the outer perimeter.** The original arc generator anchors the first arc on the previous layer's outer wall. Internal bridges span over voids in the model interior and never touch that wall, so arc generation failed. We added a fallback in `makeStartLineString` that also considers inner-wall rings.
- **"Internal web" bridges over sparse infill.** Bambu marks these as `; FEATURE: Bridge` even though they're not really supporting a top surface. The top-surface chain check + coverage ratio filter narrow the scope to bridges that genuinely hold up solid material.
- **Bridge polygon shape sensitive to slicer infill direction.** Buffering the gcode LineString alone leaves gaps perpendicular to the bridge fill direction; the closing operation merges them.
- **Removed bridge gcode left layers above unsupported.** When the BFS only filled part of a converted bridge, the rest had its bridge gcode deleted and nothing replacing it. Preserved-bridge fallback fixes this — the layer above always has either an arc or original bridge underneath.

## Limitations

- Single-plate runtime is ~3× slower than Bambu's slicer. The remaining headroom is mostly inside shapely's per-call overhead which is hard to chip away at without a deeper rewrite.
- The closing-radius extension (the ~1 mm strip beyond where the original bridge gcode actually was) can end up in the "truly unsupported" bucket on rare layers. Lower `--bridge-closing` if it bothers you, or accept that a 1 mm edge typically rests on the inner wall and not on bridge fill anyway.
- Print correctness comes from preview, not just from the script running clean. Always inspect the modified `.3mf` in Bambu Studio's gcode preview before sending to the printer.

## License

Inherits the upstream's license — see [LICENSE](LICENSE).
