# 3DPrint Antenna Inductor Maker

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python: 3.10+](https://img.shields.io/badge/Python-3.10%2B-brightgreen.svg)](https://www.python.org/)
[![Amateur Radio](https://img.shields.io/badge/Amateur%20Radio-KB1U-red.svg)](https://www.qrz.com/)

A parametric 3D CAD designer and instant **direct-to-STL exporter** for high-efficiency RF loading coil forms. Built specifically for amateur radio experimenters, antenna builders, and makers.

Whether you are building a shortened delta loop, vertical, mobile whip, or portable wire antenna, this tool calculates the exact winding dimensions and generates a clean, **100% watertight, manifold 3D-printable model** in seconds—without needing SolidWorks, OpenSCAD, or any CAD middleman.

---

## Features

- **Direct 3D Printable STL Export**: Uses vectorized Continuous Signed Distance Fields (SDF) and VTK Flying Edges to produce sub-millimeter, manifold STL models with **0 open edges** in under 1 second.
- **Interactive 3D Preview**:
  - Full hardware-accelerated 3D viewport rendered via PyVista & VTK.
  - Studio 3-point lighting rig with hardware shadow mapping to clearly reveal helical grooves and through-hole depth.
  - Smooth, gimbal-lock-free turntable orbit (Left-Click Drag), pan (Middle/Shift Drag), and zoom (Right-Click Drag or Scroll).
  - Dual view modes: **Form + Wire Winding** (realistic copper wire preview) or **Form Only (STL)**.
  - Quick-view camera presets: *Isometric*, *Front*, *Side*, *Top*, and *Reset*.
- **Wire Anchoring Holes**: Automatically generates radial through-the-wall anchor holes at the start and end of the coil. Push wire through from the inside bore to lock windings firmly under field tension.
- **Suspension Through-Holes**: Generates full through-holes at both end margins to suspend the coil inline with Dacron rope, paracord, or zip-ties.
- **Switch / Tap Mounting Hole**: Dedicated radial mounting hole in the lower margin for band switches, banana plugs, or tap posts.
- **Standard Conductor Presets**: Built-in AWG wire gauges (10–24 AWG) with typical bare and insulated (PVC/PTFE/THHN) diameters, plus custom wire OD support.
- **Multi-CAD Compatibility**: In addition to direct STL export, also exports:
  - **OpenSCAD (`.scad`)** parametric scripts.
  - **SolidWorks VBA Macro (`.bas`)** for automated native feature generation.
  - One-click clipboard copy for SolidWorks.

---

## RF Inductance Calculation

Inductance calculations use Wheeler's classical single-layer air-core solenoid equation:

$$L (\mu\text{H}) = \frac{d^2 \cdot n^2}{18d + 40l}$$

Where:
- $d$ = mean coil diameter in inches (form diameter + wire OD)
- $l$ = winding length in inches ($n \times \text{pitch}$)
- $n$ = number of turns

The quadratic formula solves for the required turns ($n$) given your target inductance and groove pitch factor.

---

## Installation & Requirements

### 1. Clone the Repository
```bash
git clone https://github.com/DC-Light/3DPrint-Antenna-Inductor-Maker.git
cd 3DPrint-Antenna-Inductor-Maker
```

### 2. Install Dependencies
Make sure you have Python 3.10+ installed. Install the required libraries using pip:

```bash
pip install -r requirements.txt
```

*(Tkinter is included with standard Python installations on Windows and macOS).*

### 3. Run the Application
```bash
python 3DPrint_Inductor_Maker.py
```

---

## 3D Printing Recommendations

| Parameter | Recommended Setting |
| :--- | :--- |
| **Material** | **PETG**, **ASA**, or **ABS** (Avoid PLA for outdoor antennas due to UV and heat creep). |
| **Orientation** | **Vertical on the print bed** (standing up). No supports required! |
| **Perimeters / Walls** | **4 to 6 walls** (ensures solid mechanical strength against wire tension). |
| **Infill** | **30% – 50% Gyroid** for uniform rigidity. |
| **Layer Height** | **0.16 mm – 0.20 mm** for smooth groove seating. |
| **Seam Placement** | Set to **Random** or aligned away from the wire grooves. |

---

## Controls & Navigation

- **Left-Click + Drag**: Rotate / orbit around the form.
- **Middle-Click + Drag** (or **Shift + Left-Drag**): Pan camera along view plane.
- **Right-Click + Drag** or **Mouse Wheel**: Zoom in / out.
- **Buttons**:
  - **Isometric / Front / Side / Top / Reset**: Instant camera view realignment.
  - **⚡ Calculate Coil & Update 3D Preview**: Recalculates turns, lengths, and regenerates 3D model.
  - **💾 Export 3D Printable STL**: Exports production-ready binary `.stl` file.
  - **ℹ️ About / License**: View author credits and MIT license.
  - **❓ Help / Guide**: Full technical guide, parameter explanations, and antenna theory.

---

## Author & License

- **Author**: Andrew Freeston / **KB1U** / DC-LIGHT LLC
- **License**: [MIT License](LICENSE)

*Created to encourage amateur radio operators, experimenters, and makers everywhere to build, print, and tune their own antennas!*
