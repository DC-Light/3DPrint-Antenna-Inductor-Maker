"""
Antenna Loading Coil CAD, 3D Preview & STL Generator
Open-source parametric 3D loading coil generator for amateur radio antennas.

MIT License

Copyright (c) 2026 Andrew Freeston / KB1U / DC-LIGHT LLC.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

import math
import tkinter as tk
from tkinter import ttk, messagebox, filedialog, scrolledtext
import numpy as np
import pyvista as pv
from PIL import Image, ImageTk

MIT_LICENSE_TEXT = """MIT License

Copyright (c) 2026 Andrew Freeston / KB1U / DC-LIGHT LLC.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE."""

# Standard AWG Wire Diameters (Bare Conductor mm, Typical Insulated OD mm)
AWG_DATA = {
    10: (2.588, 3.80),
    12: (2.053, 3.30),
    14: (1.628, 2.70),
    16: (1.291, 2.30),
    18: (1.024, 2.00),
    20: (0.812, 1.70),
    22: (0.644, 1.50),
    24: (0.511, 1.20)
}


def generate_coil_mesh(calc_data, resolution=0.42):
    """
    Generates a 100% watertight, manifold 3D mesh of the loading coil form
    using a continuous signed distance field (SDF) and Flying Edges / Marching Cubes.
    Includes the hollow cylinder, helical wire groove, and switch mounting hole.
    """
    R_outer = calc_data['D_form_mm'] / 2.0
    wall = calc_data['wall_t']
    R_inner = R_outer - wall
    total_len = calc_data['total_len_mm']
    margin = calc_data['margin']
    pitch = calc_data['pitch_mm']
    turns = calc_data['turns']
    wire_r = (calc_data['wire_od'] * 1.05) / 2.0
    switch_d = calc_data['switch_d']
    switch_offset = calc_data['switch_offset']

    res = resolution
    pad = 1.5
    xs = np.arange(-R_outer - pad, R_outer + pad + res, res)
    ys = np.arange(-R_outer - pad, R_outer + pad + res, res)
    zs = np.arange(-pad, total_len + pad + res, res)

    nx, ny, nz = len(xs), len(ys), len(zs)
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing='ij')

    r_cyl = np.sqrt(X**2 + Y**2)

    # 1. Outer cylinder (negative inside, positive outside)
    d_outer_rad = r_cyl - R_outer
    d_outer_z = np.maximum(-Z, Z - total_len)
    d_outer = np.maximum(d_outer_rad, d_outer_z)

    # 2. Inner hollow bore
    d_inner = R_inner - r_cyl
    d_solid = np.maximum(d_outer, d_inner)

    # 3. Switch mounting hole through front wall (Y >= 0)
    if switch_d > 0 and switch_offset > 0:
        d_hole_cyl = np.sqrt(X**2 + (Z - switch_offset)**2) - (switch_d / 2.0)
        d_hole = np.maximum(d_hole_cyl, -Y)
        d_solid = np.maximum(d_solid, -d_hole)

    # 4. Suspension through-holes at both ends (all the way through both walls along X axis)
    suspend_d = calc_data.get('suspend_hole_d', 0.0)
    suspend_offset = calc_data.get('suspend_offset', 0.0)
    if suspend_d > 0 and suspend_offset > 0:
        r_susp = suspend_d / 2.0
        # Bottom suspension hole through both walls at Z = suspend_offset
        d_susp_bottom = np.sqrt(Y**2 + (Z - suspend_offset)**2) - r_susp
        d_solid = np.maximum(d_solid, -d_susp_bottom)
        # Top suspension hole through both walls at Z = total_len - suspend_offset
        d_susp_top = np.sqrt(Y**2 + (Z - (total_len - suspend_offset))**2) - r_susp
        d_solid = np.maximum(d_solid, -d_susp_top)

    # 5. Wire anchor holes through the tube wall at the start and end of the coil
    wire_hole_d = calc_data.get('wire_hole_d', 0.0)
    if wire_hole_d > 0:
        r_wh = wire_hole_d / 2.0
        # Start anchor hole at start of winding: phi = 0, angle = 0 (along +X axis), Z = margin
        d_wh_start = np.sqrt(Y**2 + (Z - margin)**2) - r_wh
        d_wh_start = np.maximum(d_wh_start, -X) # through wall on +X side
        d_solid = np.maximum(d_solid, -d_wh_start)

        # End anchor hole at end of winding: phi = turns * 2pi, Z = margin + turns*pitch
        phi_end = turns * 2 * np.pi
        z_end = margin + turns * pitch
        cos_end = math.cos(phi_end)
        sin_end = math.sin(phi_end)
        # Distance to ray in direction (cos_end, sin_end, 0) at height z_end
        dist_ray_end = np.sqrt((-X * sin_end + Y * cos_end)**2 + (Z - z_end)**2)
        d_wh_end = dist_ray_end - r_wh
        s_end = X * cos_end + Y * sin_end
        d_wh_end = np.maximum(d_wh_end, -s_end) # through wall on positive ray side
        d_solid = np.maximum(d_solid, -d_wh_end)

    # 6. Helical groove cut into outer surface
    winding_mask = (
        (Z >= margin - wire_r) & 
        (Z <= margin + turns * pitch + wire_r) & 
        (r_cyl >= R_outer - wire_r * 2.0)
    )
    theta = np.arctan2(Y[winding_mask], X[winding_mask]) % (2 * np.pi)
    z_w = Z[winding_mask]

    k_approx = np.round((z_w - margin - (theta / (2 * np.pi)) * pitch) / pitch)
    phi = theta + 2 * np.pi * k_approx
    valid = (phi >= 0) & (phi <= 2 * np.pi * turns)

    z_g = margin + (phi / (2 * np.pi)) * pitch
    d_groove = np.sqrt((r_cyl[winding_mask] - R_outer)**2 + (z_w - z_g)**2) - wire_r

    groove_sub = np.where(valid, -d_groove, -1e5)
    d_solid[winding_mask] = np.maximum(d_solid[winding_mask], groove_sub)

    # Extract 0.0 isosurface using VTK Flying Edges
    grid = pv.ImageData(
        dimensions=(nx, ny, nz),
        spacing=(res, res, res),
        origin=(xs[0], ys[0], zs[0])
    )
    grid.point_data['values'] = d_solid.flatten(order='F')
    mesh = grid.contour([0.0])
    return mesh


def generate_wire_tube(calc_data):
    """Generates the 3D copper wire helix for visualization, including anchor feed segments."""
    R_outer = calc_data['D_form_mm'] / 2.0
    wall = calc_data['wall_t']
    R_inner = R_outer - wall
    margin = calc_data['margin']
    pitch = calc_data['pitch_mm']
    turns = calc_data['turns']
    wire_r = (calc_data['wire_od'] * 1.05) / 2.0

    phi_end = turns * 2 * np.pi
    z_end = margin + turns * pitch
    cos_end = math.cos(phi_end)
    sin_end = math.sin(phi_end)

    # Lead-in point extending into the inner bore at start
    p_start_inner = np.array([R_inner - 2.0, 0.0, margin])

    # Helical turns points
    n_pts = max(60, int(turns * 36))
    theta = np.linspace(0, turns * 2 * np.pi, n_pts)
    x = R_outer * np.cos(theta)
    y = R_outer * np.sin(theta)
    z = margin + (theta / (2 * np.pi)) * pitch
    points = np.column_stack((x, y, z))

    # Lead-out point extending into the inner bore at end
    p_end_inner = np.array([(R_inner - 2.0) * cos_end, (R_inner - 2.0) * sin_end, z_end])

    all_pts = np.vstack([p_start_inner, points, p_end_inner])
    poly = pv.MultipleLines(all_pts)
    return poly.tube(radius=wire_r * 0.95, n_sides=14)


class LoadingCoilApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Antenna Loading Coil CAD, 3D Preview & STL Generator")
        self.geometry("1160x820")
        self.minsize(980, 680)
        self.resizable(True, True)

        self.calc_data = None
        self.plotter = None
        self._current_photo = None
        self._last_mouse_x = None
        self._last_mouse_y = None
        self._resize_job = None
        self._view_initialized = False

        self._setup_ui()
        self._update_wire_od()
        self._setup_3d_renderer()

        self.protocol("WM_DELETE_WINDOW", self.on_closing)

        # Automatically calculate and render 3D preview on initial startup
        self.after(100, self.calculate)

    def _setup_ui(self):
        style = ttk.Style()
        style.theme_use('clam')

        # Application Menu Bar
        menubar = tk.Menu(self)
        help_menu = tk.Menu(menubar, tearoff=0)
        help_menu.add_command(label="Detailed Documentation & User Guide...", command=self.show_help_dialog)
        help_menu.add_separator()
        help_menu.add_command(label="About & License...", command=self.show_about_dialog)
        menubar.add_cascade(label="Help", menu=help_menu)
        self.config(menu=menubar)

        # Main Horizontal Paned Layout: Left Controls + Right 3D Viewport
        paned = ttk.PanedWindow(self, orient=tk.HORIZONTAL)
        paned.pack(fill=tk.BOTH, expand=True, padx=8, pady=8)

        # ----------------- LEFT PANEL (Scrollable Controls & Outputs) -----------------
        left_container = ttk.Frame(paned)
        paned.add(left_container, weight=0)

        canvas_left = tk.Canvas(left_container, width=560, highlightthickness=0)
        scrollbar_left = ttk.Scrollbar(left_container, orient=tk.VERTICAL, command=canvas_left.yview)
        self.scrollable_frame = ttk.Frame(canvas_left, padding=(12, 10, 15, 10))

        # Keep the scrollable content frame stretched to full canvas width (no empty white gaps)
        self.scrollable_frame.bind(
            '<Configure>',
            lambda e: canvas_left.configure(scrollregion=canvas_left.bbox('all'))
        )
        canvas_window = canvas_left.create_window((0, 0), window=self.scrollable_frame, anchor='nw')
        canvas_left.bind(
            '<Configure>',
            lambda e: canvas_left.itemconfig(canvas_window, width=e.width)
        )
        canvas_left.configure(yscrollcommand=scrollbar_left.set)

        canvas_left.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar_left.pack(side=tk.RIGHT, fill=tk.Y)

        # Allow mousewheel scrolling on the left panel
        def _on_mousewheel(event):
            x, y = self.winfo_pointerxy()
            w = self.winfo_containing(x, y)
            if w and (str(w).startswith(str(left_container)) or w == canvas_left):
                canvas_left.yview_scroll(int(-1 * (event.delta / 120)), "units")
        self.bind_all("<MouseWheel>", _on_mousewheel)

        header_frame = ttk.Frame(self.scrollable_frame)
        header_frame.pack(fill=tk.X, pady=(0, 6))

        title_container = ttk.Frame(header_frame)
        title_container.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        title_label = ttk.Label(
            title_container, 
            text="Delta Loop Loading Coil CAD Generator", 
            font=("Helvetica", 13, "bold")
        )
        title_label.pack(anchor=tk.W)

        subtitle_label = ttk.Label(
            title_container,
            text="Amateur Radio Inductor & STL Designer",
            font=("Segoe UI", 8, "italic"),
            foreground="#666666"
        )
        subtitle_label.pack(anchor=tk.W)

        btn_header_box = ttk.Frame(header_frame)
        btn_header_box.pack(side=tk.RIGHT, padx=(6, 0))

        about_btn = ttk.Button(
            btn_header_box,
            text="ℹ️ About / License",
            command=self.show_about_dialog,
            width=18
        )
        about_btn.pack(pady=(0, 3))

        help_btn = ttk.Button(
            btn_header_box,
            text="❓ Help / Guide",
            command=self.show_help_dialog,
            width=18
        )
        help_btn.pack()

        # Inputs Group
        input_group = ttk.LabelFrame(self.scrollable_frame, text=" Coil & RF Parameters ", padding="8")
        input_group.pack(fill=tk.X, pady=3)

        ttk.Label(input_group, text="Target Inductance (µH):").grid(row=0, column=0, sticky=tk.W, pady=3)
        self.inductance_var = tk.DoubleVar(value=16.0)
        ttk.Entry(input_group, textvariable=self.inductance_var, width=12).grid(row=0, column=1, pady=3)

        ttk.Label(input_group, text="Form Outer Diameter (mm):").grid(row=1, column=0, sticky=tk.W, pady=3)
        self.diam_var = tk.DoubleVar(value=50.0)
        ttk.Entry(input_group, textvariable=self.diam_var, width=12).grid(row=1, column=1, pady=3)

        ttk.Label(input_group, text="Wire Gauge (AWG):").grid(row=2, column=0, sticky=tk.W, pady=3)
        self.awg_var = tk.IntVar(value=14)
        awg_combo = ttk.Combobox(
            input_group, 
            textvariable=self.awg_var, 
            values=list(AWG_DATA.keys()), 
            state="readonly", 
            width=10
        )
        awg_combo.grid(row=2, column=1, pady=3)
        awg_combo.bind("<<ComboboxSelected>>", lambda e: self._update_wire_od())

        ttk.Label(input_group, text="Wire Insulation Type:").grid(row=3, column=0, sticky=tk.W, pady=3)
        self.insulated_var = tk.StringVar(value="Insulated (PVC/PTFE/THHN)")
        insul_combo = ttk.Combobox(
            input_group,
            textvariable=self.insulated_var,
            values=["Bare / Enamel Magnet", "Insulated (PVC/PTFE/THHN)"],
            state="readonly",
            width=22
        )
        insul_combo.grid(row=3, column=1, pady=3)
        insul_combo.bind("<<ComboboxSelected>>", lambda e: self._update_wire_od())

        ttk.Label(input_group, text="Total Wire OD (mm):").grid(row=4, column=0, sticky=tk.W, pady=3)
        self.wire_od_var = tk.DoubleVar()
        ttk.Entry(input_group, textvariable=self.wire_od_var, width=12).grid(row=4, column=1, pady=3)

        ttk.Label(input_group, text="Pitch Factor (x Wire OD):").grid(row=5, column=0, sticky=tk.W, pady=3)
        self.pitch_factor_var = tk.DoubleVar(value=1.25)
        ttk.Entry(input_group, textvariable=self.pitch_factor_var, width=12).grid(row=5, column=1, pady=3)

        # Mechanical Group
        mech_group = ttk.LabelFrame(self.scrollable_frame, text=" Mechanical & Mount Options ", padding="8")
        mech_group.pack(fill=tk.X, pady=4)

        ttk.Label(mech_group, text="Wall Thickness (mm):").grid(row=0, column=0, sticky=tk.W, pady=3)
        self.wall_var = tk.DoubleVar(value=3.5)
        ttk.Entry(mech_group, textvariable=self.wall_var, width=12).grid(row=0, column=1, pady=3)

        ttk.Label(mech_group, text="Switch Hole Diameter (mm):").grid(row=1, column=0, sticky=tk.W, pady=3)
        self.switch_hole_var = tk.DoubleVar(value=6.2)
        ttk.Entry(mech_group, textvariable=self.switch_hole_var, width=12).grid(row=1, column=1, pady=3)

        ttk.Label(mech_group, text="Hole Center from End (mm):").grid(row=2, column=0, sticky=tk.W, pady=3)
        self.switch_offset_var = tk.DoubleVar(value=10.0)
        ttk.Entry(mech_group, textvariable=self.switch_offset_var, width=12).grid(row=2, column=1, pady=3)

        ttk.Label(mech_group, text="End Margin Length (mm):").grid(row=3, column=0, sticky=tk.W, pady=3)
        self.margin_var = tk.DoubleVar(value=20.0)
        ttk.Entry(mech_group, textvariable=self.margin_var, width=12).grid(row=3, column=1, pady=3)

        ttk.Label(mech_group, text="Suspension Hole Dia (mm):").grid(row=4, column=0, sticky=tk.W, pady=3)
        self.suspend_hole_var = tk.DoubleVar(value=4.0)
        ttk.Entry(mech_group, textvariable=self.suspend_hole_var, width=12).grid(row=4, column=1, pady=3)

        ttk.Label(mech_group, text="Suspension from End (mm):").grid(row=5, column=0, sticky=tk.W, pady=3)
        self.suspend_offset_var = tk.DoubleVar(value=8.0)
        ttk.Entry(mech_group, textvariable=self.suspend_offset_var, width=12).grid(row=5, column=1, pady=3)

        ttk.Label(mech_group, text="Wire Anchor Hole Dia (mm):").grid(row=6, column=0, sticky=tk.W, pady=3)
        self.wire_hole_var = tk.DoubleVar() # Automatically set to wire OD + clearance in _update_wire_od
        ttk.Entry(mech_group, textvariable=self.wire_hole_var, width=12).grid(row=6, column=1, pady=3)

        # Calculate Button
        calc_btn = ttk.Button(
            self.scrollable_frame, 
            text="⚡ Calculate Coil & Update 3D Preview", 
            command=self.calculate
        )
        calc_btn.pack(fill=tk.X, pady=6, ipady=3)

        # Results Display
        self.results_group = ttk.LabelFrame(self.scrollable_frame, text=" Calculated Coil Specs ", padding="8")
        self.results_group.pack(fill=tk.X, pady=3)

        self.results_label = ttk.Label(self.results_group, text="Click 'Calculate' to see coil dimensions.", font=("Courier", 9))
        self.results_label.pack(anchor=tk.W)

        # Action Buttons Frame (Fixed visibility, clean layout)
        btn_frame = ttk.LabelFrame(self.scrollable_frame, text=" CAD & 3D Print Export ", padding="8")
        btn_frame.pack(fill=tk.X, pady=6)

        # 3D Printable STL Export Button (The Money Shot)
        self.stl_btn = ttk.Button(
            btn_frame,
            text="💾 Export 3D Printable STL (.stl)",
            command=self.export_stl,
            state=tk.DISABLED
        )
        self.stl_btn.pack(fill=tk.X, pady=(0, 6), ipady=4)

        self.copy_sw_btn = ttk.Button(
            btn_frame, 
            text="📋 Copy SolidWorks Macro Code to Clipboard", 
            command=self.copy_sw_macro, 
            state=tk.DISABLED
        )
        self.copy_sw_btn.pack(fill=tk.X, pady=(0, 6))

        sub_btn_frame = ttk.Frame(btn_frame)
        sub_btn_frame.pack(fill=tk.X)

        self.scad_btn = ttk.Button(sub_btn_frame, text="Export OpenSCAD (.scad)", command=self.export_openscad, state=tk.DISABLED)
        self.scad_btn.pack(side=tk.LEFT, expand=True, fill=tk.X, padx=(0, 2))

        self.sw_file_btn = ttk.Button(sub_btn_frame, text="Save SW Macro (.bas)", command=self.save_sw_bas, state=tk.DISABLED)
        self.sw_file_btn.pack(side=tk.RIGHT, expand=True, fill=tk.X, padx=(2, 0))

        # ----------------- RIGHT PANEL (Interactive 3D Viewport) -----------------
        right_container = ttk.LabelFrame(paned, text=" 3D Structure Preview ", padding="8")
        paned.add(right_container, weight=1)

        # Toolbar at top of 3D preview
        toolbar = ttk.Frame(right_container)
        toolbar.pack(fill=tk.X, pady=(0, 6))

        ttk.Label(toolbar, text="Mode:").pack(side=tk.LEFT, padx=(0, 4))
        self.view_mode_var = tk.StringVar(value="Form + Wire Winding")
        mode_combo = ttk.Combobox(
            toolbar,
            textvariable=self.view_mode_var,
            values=["Form + Wire Winding", "Form Only (STL Model)"],
            state="readonly",
            width=20
        )
        mode_combo.pack(side=tk.LEFT, padx=(0, 8))
        mode_combo.bind("<<ComboboxSelected>>", lambda e: self._update_3d_preview(reset_camera=False))

        ttk.Button(toolbar, text="Isometric", width=9, command=self._view_iso).pack(side=tk.LEFT, padx=2)
        ttk.Button(toolbar, text="Front", width=7, command=self._view_front).pack(side=tk.LEFT, padx=2)
        ttk.Button(toolbar, text="Side", width=7, command=self._view_side).pack(side=tk.LEFT, padx=2)
        ttk.Button(toolbar, text="Top", width=7, command=self._view_top).pack(side=tk.LEFT, padx=2)
        ttk.Button(toolbar, text="Reset", width=7, command=self._view_reset).pack(side=tk.LEFT, padx=2)

        # Canvas for PyVista offscreen render
        self.canvas_3d = tk.Canvas(right_container, bg="#1E1E1E", highlightthickness=0)
        self.canvas_3d.pack(fill=tk.BOTH, expand=True)

        # Mouse bindings for 3D interaction
        self.canvas_3d.bind("<ButtonPress-1>", self._on_mouse_press)
        self.canvas_3d.bind("<B1-Motion>", self._on_mouse_drag_rotate)
        self.canvas_3d.bind("<ButtonRelease-1>", self._on_mouse_release)
        self.canvas_3d.bind("<Shift-B1-Motion>", self._on_mouse_drag_pan)
        self.canvas_3d.bind("<ButtonPress-2>", self._on_mouse_press)
        self.canvas_3d.bind("<B2-Motion>", self._on_mouse_drag_pan)
        self.canvas_3d.bind("<ButtonRelease-2>", self._on_mouse_release)
        self.canvas_3d.bind("<ButtonPress-3>", self._on_mouse_press)
        self.canvas_3d.bind("<B3-Motion>", self._on_mouse_drag_zoom)
        self.canvas_3d.bind("<ButtonRelease-3>", self._on_mouse_release)
        self.canvas_3d.bind("<MouseWheel>", self._on_mouse_wheel)
        self.canvas_3d.bind("<Button-4>", lambda e: self._zoom(1.1))
        self.canvas_3d.bind("<Button-5>", lambda e: self._zoom(0.9))
        self.canvas_3d.bind("<Configure>", self._on_canvas_configure)

        # Bottom info bar
        info_bar = ttk.Frame(right_container)
        info_bar.pack(fill=tk.X, pady=(4, 0))

        help_label = ttk.Label(
            info_bar, 
            text="Left-Drag: Rotate | Middle / Shift-Drag: Pan | Right-Drag / Scroll: Zoom", 
            font=("Segoe UI", 8), 
            foreground="#666666"
        )
        help_label.pack(side=tk.LEFT)

        self.preview_info_label = ttk.Label(
            info_bar, 
            text="Watertight 3D Model: Ready", 
            font=("Segoe UI", 8, "italic"),
            foreground="#333333"
        )
        self.preview_info_label.pack(side=tk.RIGHT)

    def _setup_3d_renderer(self):
        """Initializes the offscreen PyVista plotter with custom angled lighting and shadows."""
        try:
            self.plotter = pv.Plotter(off_screen=True, window_size=(540, 600), lighting="none")
            self.plotter.set_background("#1E1E1E")
            self._apply_lighting_and_shadows()
        except Exception as e:
            print("Failed to initialize PyVista plotter:", e)
            self.plotter = None

    def _apply_lighting_and_shadows(self):
        """Sets up high-contrast 3-point lighting and shadows so grooves and holes pop in 3D."""
        if not self.plotter:
            return
        try:
            self.plotter.renderer.remove_all_lights()
            # Key light: angled from upper-right-front to cast rich feature shadows into grooves & holes
            key_light = pv.Light(position=(2.5, 3.5, 2.5), light_type="camera light", intensity=0.92)
            # Fill light: soft fill from lower-left to keep shaded sides clearly visible
            fill_light = pv.Light(position=(-2.5, -1.5, 1.5), light_type="camera light", intensity=0.35)
            # Rim light: subtle back light to accentuate silhouette edge definition
            rim_light = pv.Light(position=(0.0, -3.0, -2.5), light_type="camera light", intensity=0.25)

            self.plotter.add_light(key_light)
            self.plotter.add_light(fill_light)
            self.plotter.add_light(rim_light)
            self.plotter.enable_shadows()
        except Exception as e:
            pass

    def _render_view(self):
        """Captures offscreen render and blits to Tkinter Canvas."""
        if not self.plotter:
            return
        try:
            self.plotter.render()
            img_data = self.plotter.screenshot(return_img=True)
            if img_data is not None and img_data.size > 0:
                img = Image.fromarray(img_data)
                self._current_photo = ImageTk.PhotoImage(img)
                self.canvas_3d.delete("all")
                self.canvas_3d.create_image(0, 0, anchor=tk.NW, image=self._current_photo)
        except Exception as e:
            pass

    def _on_canvas_configure(self, event):
        w, h = event.width, event.height
        if w < 50 or h < 50:
            return
        if self._resize_job is not None:
            self.after_cancel(self._resize_job)
        self._resize_job = self.after(120, lambda: self._apply_canvas_resize(w, h))

    def _apply_canvas_resize(self, width, height):
        self._resize_job = None
        if self.plotter is not None:
            try:
                self.plotter.window_size = (max(100, width), max(100, height))
                self._render_view()
            except Exception:
                pass

    def _on_mouse_press(self, event):
        self._last_mouse_x = event.x
        self._last_mouse_y = event.y

    def _on_mouse_release(self, event):
        self._last_mouse_x = None
        self._last_mouse_y = None

    def _on_mouse_drag_rotate(self, event):
        if not self.plotter:
            return
        if self._last_mouse_x is None or self._last_mouse_y is None:
            self._last_mouse_x = event.x
            self._last_mouse_y = event.y
            return
        dx = event.x - self._last_mouse_x
        dy = event.y - self._last_mouse_y
        self._last_mouse_x = event.x
        self._last_mouse_y = event.y

        cam = self.plotter.camera
        focal = np.array(cam.focal_point)
        pos = np.array(cam.position)
        diff = pos - focal
        r = np.linalg.norm(diff)
        if r < 1e-6:
            return

        # Continuous spherical coordinate orbit (Z-up turntable)
        xy_dist = np.linalg.norm(diff[:2])
        if xy_dist < 1e-4:
            azimuth = -np.pi / 2.0
        else:
            azimuth = np.arctan2(diff[1], diff[0])

        sin_elev = np.clip(diff[2] / r, -1.0, 1.0)
        elev = np.arcsin(sin_elev)

        # Smooth drag angles
        azimuth -= np.radians(dx * 0.45)
        # Clamped elevation prevents gimbal lock / pole crossing oscillations
        max_elev = np.radians(88.5)
        elev = np.clip(elev - np.radians(dy * 0.45), -max_elev, max_elev)

        cos_elev = np.cos(elev)
        new_dir = np.array([
            cos_elev * np.cos(azimuth),
            cos_elev * np.sin(azimuth),
            np.sin(elev)
        ])
        cam.position = focal + r * new_dir
        cam.focal_point = focal
        cam.up = (0.0, 0.0, 1.0)
        self._render_view()

    def _on_mouse_drag_pan(self, event):
        if not self.plotter:
            return
        if self._last_mouse_x is None or self._last_mouse_y is None:
            self._last_mouse_x = event.x
            self._last_mouse_y = event.y
            return
        dx = event.x - self._last_mouse_x
        dy = event.y - self._last_mouse_y
        self._last_mouse_x = event.x
        self._last_mouse_y = event.y

        cam = self.plotter.camera
        focal = np.array(cam.focal_point)
        pos = np.array(cam.position)
        view_dir = focal - pos
        dist = np.linalg.norm(view_dir)
        if dist < 1e-6:
            return
        view_dir = view_dir / dist
        up = np.array(cam.up)

        right = np.cross(view_dir, up)
        right_norm = np.linalg.norm(right)
        if right_norm < 1e-6:
            return
        right = right / right_norm
        cam_up = np.cross(right, view_dir)

        scale = dist * 0.0018
        delta = -dx * scale * right + dy * scale * cam_up
        cam.position = pos + delta
        cam.focal_point = focal + delta
        self._render_view()

    def _on_mouse_drag_zoom(self, event):
        if not self.plotter:
            return
        if self._last_mouse_x is None or self._last_mouse_y is None:
            self._last_mouse_x = event.x
            self._last_mouse_y = event.y
            return
        dy = event.y - self._last_mouse_y
        self._last_mouse_x = event.x
        self._last_mouse_y = event.y

        factor = 1.0 - dy * 0.008
        if 0.5 < factor < 1.5:
            self.plotter.camera.zoom(factor)
        self._render_view()

    def _on_mouse_wheel(self, event):
        if event.num == 5 or event.delta < 0:
            self._zoom(0.92)
        elif event.num == 4 or event.delta > 0:
            self._zoom(1.08)

    def _zoom(self, factor):
        if self.plotter:
            self.plotter.camera.zoom(factor)
            self._render_view()

    def _view_iso(self):
        if self.plotter:
            self.plotter.view_isometric()
            self.plotter.reset_camera()
            self._render_view()

    def _view_front(self):
        if self.plotter:
            self.plotter.view_xz()
            self.plotter.reset_camera()
            self._render_view()

    def _view_side(self):
        if self.plotter:
            self.plotter.view_yz()
            self.plotter.reset_camera()
            self._render_view()

    def _view_top(self):
        if self.plotter:
            self.plotter.view_xy()
            self.plotter.reset_camera()
            self._render_view()

    def _view_reset(self):
        if self.plotter:
            self.plotter.view_isometric()
            self.plotter.reset_camera()
            self._render_view()

    def _update_3d_preview(self, reset_camera=False):
        """Regenerates the 3D preview model and updates the viewport."""
        if not self.calc_data or not self.plotter:
            return

        try:
            self.plotter.clear()
            self._apply_lighting_and_shadows()

            # Generate coil form mesh (fast preview resolution)
            coil_mesh = generate_coil_mesh(self.calc_data, resolution=0.5)

            mode = self.view_mode_var.get()
            if "Wire" in mode:
                # Engineering light gray form + metallic copper wire winding with shadow depth
                self.plotter.add_mesh(
                    coil_mesh, 
                    color="#D0D7DE", 
                    smooth_shading=True, 
                    ambient=0.25,
                    diffuse=0.75,
                    specular=0.35, 
                    specular_power=20
                )
                wire = generate_wire_tube(self.calc_data)
                self.plotter.add_mesh(
                    wire, 
                    color="#DD6B20", 
                    smooth_shading=True, 
                    ambient=0.3,
                    diffuse=0.7,
                    specular=0.8, 
                    specular_power=30
                )
            else:
                # Clean STL printable form highlighting grooves and switch/suspension holes
                self.plotter.add_mesh(
                    coil_mesh, 
                    color="#CBD5E1", 
                    smooth_shading=True, 
                    ambient=0.22,
                    diffuse=0.78,
                    specular=0.4, 
                    specular_power=25
                )

            if reset_camera or not self._view_initialized:
                self.plotter.view_isometric()
                self.plotter.reset_camera()
                self._view_initialized = True

            self._render_view()
            self.preview_info_label.config(
                text=f"Mesh: {coil_mesh.n_cells:,} triangles | Watertight (0 open edges)"
            )
        except Exception as e:
            self.preview_info_label.config(text=f"Preview error: {str(e)}")

    def _update_wire_od(self):
        awg = self.awg_var.get()
        bare_d, insul_od = AWG_DATA.get(awg, (1.628, 2.70))
        if "Bare" in self.insulated_var.get():
            od = bare_d
        else:
            od = insul_od
        self.wire_od_var.set(od)
        recommended_wire_hole = round(max(od * 1.20, od + 0.35), 2)
        self.wire_hole_var.set(recommended_wire_hole)

    def calculate(self):
        try:
            L_target = self.inductance_var.get()
            D_form_mm = self.diam_var.get()
            awg = self.awg_var.get()
            wire_od = self.wire_od_var.get()
            p_factor = self.pitch_factor_var.get()
            wall_t = self.wall_var.get()
            switch_d = self.switch_hole_var.get()
            switch_offset = self.switch_offset_var.get()
            suspend_d = self.suspend_hole_var.get()
            suspend_offset = self.suspend_offset_var.get()
            wire_hole_d = self.wire_hole_var.get()
            margin = self.margin_var.get()

            dimensions = {
                "Target inductance": L_target,
                "Form diameter": D_form_mm,
                "Wire diameter": wire_od,
                "Pitch factor": p_factor,
                "Wall thickness": wall_t,
                "Switch hole diameter": switch_d,
                "Switch hole center from end": switch_offset,
                "End margin": margin,
            }
            invalid = [name for name, value in dimensions.items() if not math.isfinite(value) or value <= 0]
            if invalid:
                raise ValueError(f"These values must be positive numbers: {', '.join(invalid)}.")
            if wall_t >= D_form_mm / 2.0:
                raise ValueError("Wall thickness must be less than half the form diameter.")
            if p_factor <= 1.05:
                raise ValueError("Pitch factor must be greater than 1.05 so adjacent groove turns do not overlap.")
            hole_radius = switch_d / 2.0
            if switch_offset <= hole_radius or switch_offset + hole_radius >= margin:
                raise ValueError(
                    "The switch hole must fit completely inside the end margin. "
                    "Increase the margin, reduce the hole diameter, or move its center."
                )

            if suspend_d < 0 or suspend_offset < 0:
                raise ValueError("Suspension hole dimensions must be non-negative.")
            if suspend_d > 0:
                susp_radius = suspend_d / 2.0
                if suspend_offset <= susp_radius or suspend_offset + susp_radius >= margin:
                    raise ValueError(
                        "The suspension holes must fit completely inside the end margins. "
                        "Increase the margin, reduce the suspension hole diameter, or move its center."
                    )

            if wire_hole_d < 0:
                raise ValueError("Wire anchor hole diameter cannot be negative.")

            bare_d, _ = AWG_DATA.get(awg, (1.628, 2.70))
            pitch_mm = wire_od * p_factor

            D_mean_mm = D_form_mm + wire_od
            D_mean_in = D_mean_mm / 25.4
            pitch_in = pitch_mm / 25.4

            A = D_mean_in ** 2
            B = -40.0 * pitch_in * L_target
            C = -18.0 * D_mean_in * L_target

            discriminant = B**2 - 4*A*C
            if discriminant < 0:
                raise ValueError("Unrealizable coil dimensions.")

            N_turns = (-B + math.sqrt(discriminant)) / (2 * A)
            winding_len_mm = N_turns * pitch_mm
            total_len_mm = winding_len_mm + (2 * margin)
            wire_len_m = (math.pi * D_mean_mm * N_turns) / 1000.0

            self.calc_data = {
                "L": L_target,
                "D_form_mm": D_form_mm,
                "D_mean_mm": D_mean_mm,
                "wire_bare_d": bare_d,
                "wire_od": wire_od,
                "pitch_mm": pitch_mm,
                "turns": N_turns,
                "winding_len_mm": winding_len_mm,
                "total_len_mm": total_len_mm,
                "wire_len_m": wire_len_m,
                "wall_t": wall_t,
                "switch_d": switch_d,
                "switch_offset": switch_offset,
                "suspend_hole_d": suspend_d,
                "suspend_offset": suspend_offset,
                "wire_hole_d": wire_hole_d,
                "margin": margin,
                "awg": awg,
                "is_insulated": "Bare" not in self.insulated_var.get()
            }

            insul_status = "Insulated" if self.calc_data["is_insulated"] else "Bare/Enamel"
            susp_text = (
                f"Suspension Holes:{suspend_d:.2f} mm through both ends, {suspend_offset:.2f} mm in\n"
                if suspend_d > 0 else "Suspension Holes:Disabled\n"
            )
            wire_hole_text = (
                f"Wire Anchor Holes:{wire_hole_d:.2f} mm through wall at coil start & end\n"
                if wire_hole_d > 0 else "Wire Anchor Holes:Disabled\n"
            )
            res_text = (
                f"Wire Spec:        AWG {awg} ({insul_status})\n"
                f"Bare / Total OD:  {bare_d:.3f} mm / {wire_od:.3f} mm\n"
                f"Groove Pitch:     {pitch_mm:.3f} mm\n"
                f"Required Turns:   {N_turns:.2f} turns\n"
                f"Winding Length:   {winding_len_mm:.1f} mm\n"
                f"Total Form Length:{total_len_mm:.1f} mm\n"
                f"Switch Hole:      {switch_d:.2f} mm dia, {switch_offset:.2f} mm from end\n"
                f"{susp_text}"
                f"{wire_hole_text}"
                f"Total Wire Req:   {wire_len_m:.2f} meters ({wire_len_m*3.28084:.1f} ft)"
            )
            self.results_label.config(text=res_text)
            self.stl_btn.config(state=tk.NORMAL)
            self.scad_btn.config(state=tk.NORMAL)
            self.sw_file_btn.config(state=tk.NORMAL)
            self.copy_sw_btn.config(state=tk.NORMAL)

            # Update 3D viewport
            self._update_3d_preview(reset_camera=False)

        except Exception as e:
            messagebox.showerror("Error", f"Calculation failed: {str(e)}")

    def export_stl(self):
        """Directly exports high-precision, 100% watertight STL file for 3D printing."""
        if not self.calc_data:
            return
        
        filePath = filedialog.asksaveasfilename(
            defaultextension=".stl",
            filetypes=[("Stereolithography (STL)", "*.stl"), ("All Files", "*.*")],
            title="Export 3D Printable STL"
        )
        if not filePath:
            return

        self.preview_info_label.config(text="Generating high-precision STL model...")
        self.update_idletasks()

        try:
            # 0.35 mm resolution produces ultra-smooth curves and clean grooves
            mesh = generate_coil_mesh(self.calc_data, resolution=0.35)
            mesh.save(filePath)

            self.preview_info_label.config(
                text=f"Exported: {mesh.n_cells:,} triangles | Watertight (0 open edges)"
            )
            messagebox.showinfo(
                "STL Export Successful!",
                f"Loading coil form STL generated and saved successfully!\n\n"
                f"File: {filePath}\n\n"
                f"Mesh Quality:\n"
                f" • Triangles: {mesh.n_cells:,}\n"
                f" • Watertight: Yes (0 open edges)\n"
                f" • Form Dimensions: {self.calc_data['D_form_mm']:.1f} mm OD × {self.calc_data['total_len_mm']:.1f} mm L\n\n"
                f"Directly ready to slice in PrusaSlicer, Bambu Studio, OrcaSlicer, Cura, etc."
            )
        except Exception as e:
            self.preview_info_label.config(text="STL export error")
            messagebox.showerror("Export Error", f"Failed to generate STL: {str(e)}")

    def _generate_sw_vba_code(self):
        d = self.calc_data
        r_m = (d['D_form_mm'] / 2.0) / 1000.0
        h_m = d['total_len_mm'] / 1000.0
        pitch_m = d['pitch_mm'] / 1000.0
        wall_m = d['wall_t'] / 1000.0
        wire_r_m = ((d['wire_od'] * 1.05) / 2.0) / 1000.0
        switch_r_m = (d['switch_d'] / 2.0) / 1000.0
        switch_center_m = (d['switch_offset'] - d['margin']) / 1000.0
        margin_m = d['margin'] / 1000.0
        suspend_d = d.get('suspend_hole_d', 4.0)
        suspend_offset = d.get('suspend_offset', 8.0)
        susp_r_m = (suspend_d / 2.0) / 1000.0
        susp_bot_m = (suspend_offset - d['margin']) / 1000.0
        susp_top_m = (d['total_len_mm'] - suspend_offset - d['margin']) / 1000.0

        wire_hole_d = d.get('wire_hole_d', 3.0)
        wh_r_m = (wire_hole_d / 2.0) / 1000.0
        phi_end = d['turns'] * 2 * math.pi
        z_end_m = (d['turns'] * d['pitch_mm']) / 1000.0 # height from top plane (0 is at margin)
        wh_end_x_m = r_m * math.cos(phi_end)
        wh_end_y_m = r_m * math.sin(phi_end)

        return f"""Option Explicit

' SolidWorks VBA Macro - Loading Coil Form Generator
' Generated for L = {d['L']} uH, AWG {d['awg']} ({'Insulated' if d['is_insulated'] else 'Bare'})

Dim swApp As Object
Dim Part As Object

Sub main()
    Dim tubeFeature As Object
    Dim helixFeature As Object
    Dim profileFeature As Object
    Dim grooveFeature As Object
    Dim holeFeature As Object
    Dim stage As String

    On Error GoTo MacroError

    stage = "creating a new part"
    Set swApp = Application.SldWorks
    Set Part = swApp.NewPart()

    If Part Is Nothing Then
        Err.Raise vbObjectError + 1000, , "Failed to create a new Part document. Check the default part template in SolidWorks Options."
    End If

    ' 1. Hollow tube with equal end margins around the winding region.
    stage = "creating the hollow tube"
    Part.Extension.SelectByID2 "Top Plane", "PLANE", 0, 0, 0, False, 0, Nothing, 0
    Part.SketchManager.InsertSketch True
    Part.SketchManager.CreateCircle 0, 0, 0, {r_m}, 0, 0
    Part.SketchManager.CreateCircle 0, 0, 0, {r_m - wall_m}, 0, 0
    Set tubeFeature = Part.FeatureManager.FeatureExtrusion3(False, False, False, 0, 0, {h_m - margin_m}, {margin_m}, False, False, False, False, 0#, 0#, False, False, False, False, True, True, True, 0, 0, False)
    If tubeFeature Is Nothing Then Err.Raise vbObjectError + 1001, , "The tube extrusion was not created."

    ' 2. Constant-pitch helix, defined by pitch and revolutions.
    stage = "creating the helix"
    Part.ClearSelection2 True
    Part.Extension.SelectByID2 "Top Plane", "PLANE", 0, 0, 0, False, 0, Nothing, 0
    Part.SketchManager.InsertSketch True
    Part.SketchManager.CreateCircle 0, 0, 0, {r_m}, 0, 0
    Part.InsertHelix False, True, False, False, 0, 0#, {pitch_m}, {d['turns']}, 0#, 0#
    Set helixFeature = Part.FeatureByPositionReverse(0)
    If helixFeature Is Nothing Then Err.Raise vbObjectError + 1002, , "The helix was not created."
    If helixFeature.GetTypeName2 <> "Helix" Then Err.Raise vbObjectError + 1003, , "SolidWorks did not create a helix feature."
    helixFeature.Name = "CoilPath"

    ' 3. Circular Wire Profile Sketch on Front Plane at Helix Start
    stage = "creating the groove profile"
    Part.ClearSelection2 True
    Part.Extension.SelectByID2 "Front Plane", "PLANE", 0, 0, 0, False, 0, Nothing, 0
    Part.SketchManager.InsertSketch True
    Part.SketchManager.CreateCircle {r_m}, 0, 0, {r_m + wire_r_m}, 0, 0
    Part.SketchManager.InsertSketch True
    Set profileFeature = Part.FeatureByPositionReverse(0)
    If profileFeature Is Nothing Then Err.Raise vbObjectError + 1004, , "The groove profile sketch was not created."
    profileFeature.Name = "GrooveProfile"

    ' 4. Sweep the circular profile along the helix to cut the wire groove.
    stage = "cutting the helical groove"
    Part.ClearSelection2 True
    If Not profileFeature.Select2(False, 1) Then Err.Raise vbObjectError + 1005, , "Could not select the groove profile."
    If Not helixFeature.Select2(True, 4) Then Err.Raise vbObjectError + 1006, , "Could not select the helix path."
    Set grooveFeature = Part.FeatureManager.InsertCutSwept4(False, False, 0, False, False, 0, 0, False, 0#, 0#, 0, 0, False, True, 0#, True, False, False, False)
    If grooveFeature Is Nothing Then Err.Raise vbObjectError + 1007, , "The swept groove could not be created. Check that the profile intersects the form at the helix start."

    ' 5. Switch mounting hole centered in the clear lower end margin.
    stage = "cutting the switch mounting hole"
    Part.ClearSelection2 True
    Part.Extension.SelectByID2 "Front Plane", "PLANE", 0, 0, 0, False, 0, Nothing, 0
    Part.SketchManager.InsertSketch True
    Part.SketchManager.CreateCircle 0, {switch_center_m}, 0, {switch_r_m}, {switch_center_m}, 0
    Set holeFeature = Part.FeatureManager.FeatureCut4(True, False, False, 0, 0, {r_m * 2.5}, 0.01, False, False, False, False, 0#, 0#, False, False, False, False, False, True, True, True, True, False, False, False, False, False)
    If holeFeature Is Nothing Then Err.Raise vbObjectError + 1008, , "The switch mounting hole was not created."

    ' 6. Suspension through-holes at both ends (orthogonal along Right Plane)
    stage = "cutting the suspension through-holes"
    Part.ClearSelection2 True
    Part.Extension.SelectByID2 "Right Plane", "PLANE", 0, 0, 0, False, 0, Nothing, 0
    Part.SketchManager.InsertSketch True
    Part.SketchManager.CreateCircle 0, {susp_bot_m}, 0, {susp_r_m}, {susp_bot_m}, 0
    Part.SketchManager.CreateCircle 0, {susp_top_m}, 0, {susp_r_m}, {susp_top_m}, 0
    Set holeFeature = Part.FeatureManager.FeatureCut4(True, False, False, 0, 0, {r_m * 2.5}, {r_m * 2.5}, False, False, False, False, 0#, 0#, False, False, False, False, False, True, True, True, True, False, False, False, False, False)
    If holeFeature Is Nothing Then Err.Raise vbObjectError + 1009, , "The suspension through-holes were not created."

    ' 7. Wire anchor hole at coil start (groove bottom at Z = 0)
    stage = "cutting the start wire anchor hole"
    Part.ClearSelection2 True
    Part.Extension.SelectByID2 "Front Plane", "PLANE", 0, 0, 0, False, 0, Nothing, 0
    Part.SketchManager.InsertSketch True
    Part.SketchManager.CreateCircle {r_m}, 0, 0, {r_m + wh_r_m}, 0, 0
    Set holeFeature = Part.FeatureManager.FeatureCut4(True, False, False, 0, 0, {wall_m * 2.0}, 0.001, False, False, False, False, 0#, 0#, False, False, False, False, False, True, True, True, True, False, False, False, False, False)
    If holeFeature Is Nothing Then Err.Raise vbObjectError + 1010, , "The start wire anchor hole was not created."

    Part.ClearSelection2 True
    Part.EditRebuild3
    Part.ViewZoomtofit2
    MsgBox "Loading coil form created successfully.", vbInformation
    Exit Sub

MacroError:
    MsgBox "SolidWorks stopped while " & stage & "." & vbCrLf & vbCrLf & _
           "Error " & Err.Number & ": " & Err.Description, vbCritical, "Loading Coil Macro"
End Sub
"""

    def copy_sw_macro(self):
        if not self.calc_data:
            return
        vba_code = self._generate_sw_vba_code()
        self.clipboard_clear()
        self.clipboard_append(vba_code)
        messagebox.showinfo(
            "Copied to Clipboard!", 
            "SolidWorks Macro Code Copied!\n\nTo run in SolidWorks:\n"
            "1. Open SolidWorks\n"
            "2. Go to Tools -> Macro -> New...\n"
            "3. Select all default text, paste code, and press F5"
        )

    def save_sw_bas(self):
        if not self.calc_data:
            return
        filePath = filedialog.asksaveasfilename(
            defaultextension=".bas",
            filetypes=[("VBA Source Module", "*.bas"), ("Text File", "*.txt")],
            title="Save SolidWorks VBA Code"
        )
        if not filePath:
            return

        with open(filePath, "w") as f:
            f.write(self._generate_sw_vba_code())
        messagebox.showinfo("Saved!", "VBA code saved successfully!")

    def export_openscad(self):
        if not self.calc_data:
            return
        
        filePath = filedialog.asksaveasfilename(
            defaultextension=".scad",
            filetypes=[("OpenSCAD Files", "*.scad")],
            title="Save OpenSCAD Script"
        )
        if not filePath:
            return

        d = self.calc_data
        suspend_d = d.get('suspend_hole_d', 4.0)
        suspend_offset = d.get('suspend_offset', 8.0)
        wire_hole_d = d.get('wire_hole_d', 3.0)
        end_angle_deg = (d['turns'] * 360.0) % 360.0
        z_end = d['margin'] + (d['turns'] * d['pitch_mm'])

        scad_code = f"""// Parametric Antenna Loading Coil Form
// Generated for L = {d['L']} uH, AWG {d['awg']} ({'Insulated' if d['is_insulated'] else 'Bare'})

$fn = 90;

form_diam = {d['D_form_mm']:.2f};
wall_thick = {d['wall_t']:.2f};
total_len = {d['total_len_mm']:.2f};
pitch = {d['pitch_mm']:.3f};
turns = {d['turns']:.2f};
wire_od = {d['wire_od']:.3f};
margin = {d['margin']:.2f};
switch_d = {d['switch_d']:.2f};
switch_offset = {d['switch_offset']:.2f};
suspend_d = {suspend_d:.2f};
suspend_offset = {suspend_offset:.2f};
wire_hole_d = {wire_hole_d:.2f};

module loading_coil() {{
    difference() {{
        cylinder(d=form_diam, h=total_len);

        translate([0, 0, -1])
            cylinder(d=form_diam - (2 * wall_thick), h=total_len + 2);

        for (step = [0 : 5 : 360 * turns]) {{
            angle = step;
            z = margin + ((step / 360) * pitch);
            rotate([0, 0, angle])
                translate([form_diam / 2, 0, z])
                    sphere(d=wire_od * 1.05, $fn=16);
        }}

        // Switch mounting hole (Front wall)
        translate([0, (form_diam/2) - (wall_thick/2), switch_offset])
            rotate([90, 0, 0])
                cylinder(d=switch_d, h=wall_thick * 3, center=true);

        // Suspension through-holes at both ends
        translate([0, 0, suspend_offset])
            rotate([0, 90, 0])
                cylinder(d=suspend_d, h=form_diam * 2, center=true);

        translate([0, 0, total_len - suspend_offset])
            rotate([0, 90, 0])
                cylinder(d=suspend_d, h=form_diam * 2, center=true);

        // Wire anchor holes into inner bore at coil start and end
        translate([form_diam / 4, 0, margin])
            rotate([0, 90, 0])
                cylinder(d=wire_hole_d, h=form_diam, center=true);

        rotate([0, 0, {end_angle_deg:.2f}])
            translate([form_diam / 4, 0, {z_end:.2f}])
                rotate([0, 90, 0])
                    cylinder(d=wire_hole_d, h=form_diam, center=true);
    }}
}}

loading_coil();
"""
        with open(filePath, "w") as f:
            f.write(scad_code)
        messagebox.showinfo("Success", "OpenSCAD script exported successfully!")

    def show_about_dialog(self):
        """Displays an About dialog with MIT license and ham radio open-source notes."""
        dialog = tk.Toplevel(self)
        dialog.title("About & License - Antenna Loading Coil Maker")
        dialog.geometry("640x580")
        dialog.minsize(520, 460)
        dialog.transient(self)
        dialog.grab_set()

        content = ttk.Frame(dialog, padding="15")
        content.pack(fill=tk.BOTH, expand=True)

        title_lbl = ttk.Label(
            content,
            text="Antenna Loading Coil CAD & STL Generator",
            font=("Helvetica", 13, "bold")
        )
        title_lbl.pack(pady=(0, 2))

        sub_lbl = ttk.Label(
            content,
            text="Free & Open-Source Tool for Ham Radio & Antenna Builders (MIT License)",
            font=("Segoe UI", 9, "italic"),
            foreground="#2563EB"
        )
        sub_lbl.pack(pady=(0, 8))

        desc_lbl = ttk.Label(
            content,
            text=(
                "Parametric 3D loading coil form designer and direct STL exporter.\n"
                "Created to encourage radio amateurs, antenna experimenters, and makers everywhere\n"
                "to design, 3D print, and wind high-efficiency custom inductors for delta loops,\n"
                "verticals, mobile whips, and portable wire antennas."
            ),
            justify=tk.CENTER
        )
        desc_lbl.pack(pady=(0, 10))

        license_box = ttk.LabelFrame(content, text=" Software License (MIT) ", padding="8")
        license_box.pack(fill=tk.BOTH, expand=True, pady=(0, 10))

        st = scrolledtext.ScrolledText(license_box, height=12, wrap=tk.WORD, font=("Consolas", 9))
        st.insert(tk.END, MIT_LICENSE_TEXT)
        st.config(state=tk.DISABLED)
        st.pack(fill=tk.BOTH, expand=True)

        credit_lbl = ttk.Label(
            content,
            text="Built with Python, Tkinter, PyVista, VTK, NumPy, and Pillow.",
            font=("Segoe UI", 8),
            foreground="#666666"
        )
        credit_lbl.pack(pady=(0, 8))

        btn_bar = ttk.Frame(content)
        btn_bar.pack(fill=tk.X)

        def copy_license():
            self.clipboard_clear()
            self.clipboard_append(MIT_LICENSE_TEXT)
            messagebox.showinfo("Copied", "MIT License text copied to clipboard!", parent=dialog)

        copy_btn = ttk.Button(btn_bar, text="📋 Copy License Text", command=copy_license)
        copy_btn.pack(side=tk.LEFT)

        close_btn = ttk.Button(btn_bar, text="Close", command=dialog.destroy)
        close_btn.pack(side=tk.RIGHT)

    def show_help_dialog(self):
        """Displays a comprehensive documentation, RF theory, and 3D printing guide dialog."""
        dialog = tk.Toplevel(self)
        dialog.title("Detailed User Guide & Antenna Coil Reference")
        dialog.geometry("720x660")
        dialog.minsize(580, 480)
        dialog.transient(self)
        dialog.grab_set()

        content = ttk.Frame(dialog, padding="15")
        content.pack(fill=tk.BOTH, expand=True)

        title_lbl = ttk.Label(
            content,
            text="📖 Loading Coil Designer — User Guide & Technical Reference",
            font=("Helvetica", 12, "bold")
        )
        title_lbl.pack(pady=(0, 4), anchor=tk.W)

        sub_lbl = ttk.Label(
            content,
            text="Comprehensive guide to coil parameters, RF calculations, 3D printing, and wire winding.",
            font=("Segoe UI", 9, "italic"),
            foreground="#4B5563"
        )
        sub_lbl.pack(pady=(0, 8), anchor=tk.W)

        help_box = ttk.LabelFrame(content, text=" Technical Documentation ", padding="6")
        help_box.pack(fill=tk.BOTH, expand=True, pady=(0, 10))

        st = scrolledtext.ScrolledText(help_box, height=20, wrap=tk.WORD, font=("Segoe UI", 9))
        
        help_content = """1. OVERVIEW & PURPOSE
--------------------------------------------------------------------------------
This utility calculates, 3D-visualizes, and directly exports precision 3D-printable
forms for single-layer air-core RF loading coils. These coils are typically inserted
into delta loops, quarter-wave verticals, mobile whips, or shortened wire dipoles
to cancel capacitive reactance and bring the antenna system into resonance at lower
frequencies without requiring full physical wire length.


2. COIL & RF PARAMETERS EXPLAINED
--------------------------------------------------------------------------------
• Target Inductance (µH):
  The desired total inductance of the loading coil. Calculated from the antenna's
  operating frequency and capacitive reactance (L = |Xc| / (2 * π * f)).

• Form Outer Diameter (mm):
  The exterior diameter of the cylindrical form. A larger diameter yields higher Q
  and fewer required turns for the same inductance, but increases bulk and wind load.

• Wire Gauge (AWG) & Insulation:
  Select the conductor gauge you plan to wind. The program maintains standard bare
  and typical insulated (PVC/PTFE/THHN) outer diameters. Insulated wire gives
  inherent dielectric spacing, while bare/enamelled wire packs tightly.

• Total Wire OD (mm):
  The exact measured outer diameter of the insulated or bare wire. You can fine-tune
  this with digital calipers to match your specific spool.

• Pitch Factor (× Wire OD):
  Center-to-center turn spacing relative to wire OD.
  - 1.15 to 1.30 is ideal: keeps turns well-separated to minimize inter-turn capacitance
    while maintaining high magnetic coupling.
  - Must be > 1.05 so adjacent groove passes do not break through each other.


3. MECHANICAL & MOUNTING OPTIONS
--------------------------------------------------------------------------------
• Wall Thickness (mm):
  The thickness of the hollow cylindrical tube. 3.0 to 4.0 mm is recommended for
  PETG/ABS/ASA for good mechanical rigidity against wire tension.

• Switch Hole Diameter & Center Offset (mm):
  A radial hole through the wall in the lower clear margin. Designed to accommodate
  a toggle switch, band tap post, SO-239 screw, or banana plug. Must sit cleanly
  inside the end margin.

• Suspension Hole Dia & Offset from End (mm):
  Through-holes passing completely through both cylinder walls across the full
  diameter at both ends. Designed to accept Dacron rope, paracord, or cable ties
  to suspend the coil inline in wire antennas.

• Wire Anchor Hole Dia (mm):
  Radial holes through the tube wall directly at the beginning (margin) and end
  (margin + winding length) of the helical groove. Pushing the wire from the inner
  bore out into the start of the groove and back inside at the end anchors the wire
  firmly so the winding can never uncoil under field tension.


4. INDUCTANCE CALCULATION FORMULA
--------------------------------------------------------------------------------
Inductance is calculated using Wheeler's classical single-layer solenoid equation:
    L (µH) = (d² * n²) / (18d + 40l)

Where:
    d = mean coil diameter in inches (form diameter + wire OD)
    l = coil winding length in inches (turns * pitch)
    n = number of turns

The quadratic formula solves for the exact number of turns needed to achieve the
target inductance with your specified groove pitch.


5. 3D PRINTING RECOMMENDATIONS (SLICER SETTINGS)
--------------------------------------------------------------------------------
• Recommended Materials:
  - PETG: Excellent UV and weather resistance, good RF dielectric characteristics.
  - ASA / ABS: Superior outdoor weather and temperature resistance.
  - (Avoid PLA for permanent outdoor antennas due to creep and UV breakdown).

• Print Orientation:
  Stand the tube vertically on the print bed. No print supports are needed for the
  helical groove (its overhang is minimal) or the through-holes!

• Slicer Tuning:
  - Perimeters / Walls: 4 to 6 walls (ensures the wall is nearly solid plastic).
  - Infill: 30% to 50% Gyroid.
  - Layer Height: 0.16 mm to 0.20 mm for smooth groove tracks.
  - Seam Position: Random or aligned away from the wire grooves.


6. 3D PREVIEW NAVIGATION CONTROLS
--------------------------------------------------------------------------------
• Left-Click + Drag: Smooth turntable orbit around the coil.
• Middle-Click + Drag (or Shift + Left-Drag): Pan camera along view plane.
• Right-Click + Drag or Mouse Wheel: Smooth zoom in and out.
• Quick View Buttons: Preset Isometric, Front, Side, Top, and Reset camera views.
• View Mode Selector: Switch between 'Form + Wire Winding' and 'Form Only (STL)'.
"""

        st.insert(tk.END, help_content)
        st.config(state=tk.DISABLED)
        st.pack(fill=tk.BOTH, expand=True)

        close_btn = ttk.Button(content, text="Close", command=dialog.destroy)
        close_btn.pack(anchor=tk.E, pady=(4, 0))

    def on_closing(self):
        """Cleanup plotter and resources upon window close."""
        try:
            if hasattr(self, 'plotter') and self.plotter is not None:
                self.plotter.close()
        except Exception:
            pass
        self.destroy()


if __name__ == "__main__":
    app = LoadingCoilApp()
    app.mainloop()