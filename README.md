# Arc Overhangs - PrusaSlicer & OrcaSlicer Integration

<p align="center">
<img src="https://github.com/nicolai-wachenschwan/arc-overhang-prusaslicer-integration/blob/main/examples/ExampleCatchImage.png" width=600>
</p>

This script modifies G-code to generate arc-based overhangs, enabling successful 3D printing of large 90° overhangs. Originally developed by [Steven McCulloch](https://github.com/stmcculloch), with further improvements by [Nic](https://github.com/nicolai-wachenschwan).

## 🚀 Features
- Converts overhangs into arc structures to improve printability.
- Works as a post-processing script for **PrusaSlicer** and **OrcaSlicer**.
- Automatically extracts slicer parameters from G-code.
![possible usecases](examples/UseCasesGallery.png)

## 🔧 Installation & Usage

### **Option A: Command Line Execution**
Run the script using the full path to both Python and the script:
```sh
"C:\path\to\python.exe" "C:\path\to\arc_overhangs_v1.0.0.py" "C:\path\to\input.gcode"
```

### **Option B: Slicer Integration**  
- **PrusaSlicer**: Go to **Print Settings → Output Options → Post-processing Scripts**.
- **OrcaSlicer**: Enable **Advanced Mode**, then navigate to **Others → Post-processing Scripts**.

Enter the full path to your Python executable followed by the script path.
```sh
"C:\path\to\python.exe" "C:\path\to\arc_overhangs_v1.0.0.py"
```

## 📌 Requirements
- **Python** 3.10+
- **Libraries**: `shapely`, `numpy`, `numpy-hilbert-curve`, `matplotlib`
- **Tested on**: PrusaSlicer 2.5–2.9, OrcaSlicer 2.2.0

## ⚠️ Known Issues
- Certain parameter values may cause failures (see documentation for details).

## 💡 Tips
- For advanced configuration, edit the **Parameter** section in the script.
- **Print as cold as possible.** Nic used 190 degrees for PLA. You can probably go even lower. If you require a higher temp for the rest of the print, you could insert some temp-change gcode before and after the arcs are printed. Might waste a lot of time though.
- **Maximize cooling.** Set your fans to full blast. This technique will probably not work too well with ABS and materials that can't use cooling fans, but it hasn't been tested.
- **Print slowly.** I use around 2 mm/s. Even that is too fast sometimes for the really tiny arcs since they have almost no time to cool before the next layer begins.
