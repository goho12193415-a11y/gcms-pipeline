#!/usr/bin/env python3
"""
GC-MS Auto-Processing GUI v1.0
==============================
One-click desktop application for GC-MS data processing.
Wraps the v2.1 pipeline with a simple Tkinter interface.
"""
import sys, os, threading, json
from pathlib import Path
from datetime import datetime

# Add pipeline directory to path
SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))

import tkinter as tk
from tkinter import filedialog, ttk, messagebox


class GCMSApp:
    def __init__(self, root):
        self.root = root
        self.root.title("GC-MS Auto-Processing Pipeline")
        self.root.geometry("650x500")
        self.root.resizable(True, True)

        # ---- Style ----
        self.bg = "#f5f5f5"
        self.fg = "#333333"
        self.accent = "#2F5496"
        self.root.configure(bg=self.bg)

        # Header
        header = tk.Label(root, text="GC-MS Auto-Processing Pipeline",
                         font=("Segoe UI", 16, "bold"), fg=self.accent, bg=self.bg)
        header.pack(pady=(20, 5))

        subtitle = tk.Label(root, text="RAW → Peak Detection → NIST Search → Excel Report",
                           font=("Segoe UI", 9), fg="#666", bg=self.bg)
        subtitle.pack(pady=(0, 15))

        # ---- File Selection ----
        file_frame = tk.LabelFrame(root, text="Input Files", font=("Segoe UI", 10, "bold"),
                                    bg=self.bg, fg=self.fg, padx=10, pady=10)
        file_frame.pack(fill="x", padx=20, pady=(0, 10))

        self.file_listbox = tk.Listbox(file_frame, height=4, font=("Consolas", 9))
        self.file_listbox.pack(fill="x", pady=(0, 5))

        btn_frame = tk.Frame(file_frame, bg=self.bg)
        btn_frame.pack(fill="x")

        tk.Button(btn_frame, text="+ Add .RAW Files", command=self.add_raw_files,
                 bg="#4a90d9", fg="white", font=("Segoe UI", 9),
                 padx=10).pack(side="left", padx=(0, 5))
        tk.Button(btn_frame, text="+ Add .mzML Files", command=self.add_mzml_files,
                 bg="#6c757d", fg="white", font=("Segoe UI", 9),
                 padx=10).pack(side="left", padx=(0, 5))
        tk.Button(btn_frame, text="Clear", command=self.clear_files,
                 bg="#dc3545", fg="white", font=("Segoe UI", 9),
                 padx=10).pack(side="left")

        # ---- Settings ----
        settings_frame = tk.LabelFrame(root, text="Settings", font=("Segoe UI", 10, "bold"),
                                        bg=self.bg, fg=self.fg, padx=10, pady=10)
        settings_frame.pack(fill="x", padx=20, pady=(0, 10))

        row1 = tk.Frame(settings_frame, bg=self.bg)
        row1.pack(fill="x", pady=2)
        tk.Label(row1, text="Min S/N Ratio:", bg=self.bg, width=14, anchor="w").pack(side="left")
        self.sn_var = tk.StringVar(value="10")
        tk.Entry(row1, textvariable=self.sn_var, width=8).pack(side="left", padx=(0, 20))
        tk.Label(row1, text="Min RMF:", bg=self.bg, width=10, anchor="w").pack(side="left")
        self.rmf_var = tk.StringVar(value="700")
        tk.Entry(row1, textvariable=self.rmf_var, width=8).pack(side="left")

        row2 = tk.Frame(settings_frame, bg=self.bg)
        row2.pack(fill="x", pady=2)
        tk.Label(row2, text="Output Directory:", bg=self.bg, width=14, anchor="w").pack(side="left")
        self.out_var = tk.StringVar(value=str(SCRIPT_DIR / "output"))
        tk.Entry(row2, textvariable=self.out_var, width=40).pack(side="left", padx=(0, 5))
        tk.Button(row2, text="Browse", command=self.browse_output,
                 bg="#e0e0e0", font=("Segoe UI", 8), padx=5).pack(side="left")

        row3 = tk.Frame(settings_frame, bg=self.bg)
        row3.pack(fill="x", pady=2)
        self.nist_var = tk.BooleanVar(value=True)
        tk.Checkbutton(row3, text="Use NIST MS Search engine (recommended)",
                      variable=self.nist_var, bg=self.bg, anchor="w").pack(side="left")
        self.deconv_var = tk.BooleanVar(value=False)
        tk.Checkbutton(row3, text="Enable deconvolution (co-eluting peaks)",
                      variable=self.deconv_var, bg=self.bg, anchor="w").pack(side="left", padx=(20, 0))

        # ---- Run Button ----
        self.run_btn = tk.Button(root, text="▶  Run Pipeline", command=self.run_pipeline,
                                  bg="#2F5496", fg="white", font=("Segoe UI", 12, "bold"),
                                  padx=30, pady=10, cursor="hand2")
        self.run_btn.pack(pady=(10, 5))

        # ---- Progress ----
        self.status_var = tk.StringVar(value="Ready. Add files and click Run.")
        status_label = tk.Label(root, textvariable=self.status_var,
                                font=("Segoe UI", 9), fg="#666", bg=self.bg)
        status_label.pack(pady=(0, 5))

        self.progress = ttk.Progressbar(root, mode="indeterminate", length=400)
        self.progress.pack(pady=(0, 10))

        # ---- Output ----
        self.output_text = tk.Text(root, height=8, font=("Consolas", 8),
                                    bg="#1e1e1e", fg="#d4d4d4", state="disabled")
        self.output_text.pack(fill="both", expand=True, padx=20, pady=(0, 10))

        # Redirect stdout to the text widget
        self.stdout_redirect = TextRedirector(self.output_text)
        sys.stdout = self.stdout_redirect

    def add_raw_files(self):
        files = filedialog.askopenfilenames(
            title="Select .RAW files",
            filetypes=[("Thermo RAW files", "*.raw *.RAW"), ("All files", "*.*")]
        )
        for f in files:
            if f not in self.file_listbox.get(0, tk.END):
                self.file_listbox.insert(tk.END, f)

    def add_mzml_files(self):
        files = filedialog.askopenfilenames(
            title="Select .mzML files",
            filetypes=[("mzML files", "*.mzML *.mzml"), ("All files", "*.*")]
        )
        for f in files:
            if f not in self.file_listbox.get(0, tk.END):
                self.file_listbox.insert(tk.END, f)

    def clear_files(self):
        self.file_listbox.delete(0, tk.END)

    def browse_output(self):
        d = filedialog.askdirectory(title="Select output directory")
        if d:
            self.out_var.set(d)

    def run_pipeline(self):
        files = list(self.file_listbox.get(0, tk.END))
        if not files:
            messagebox.showwarning("No Files", "Please add at least one .RAW or .mzML file.")
            return

        raw_files = [f for f in files if f.lower().endswith('.raw')]
        mzml_files = [f for f in files if f.lower().endswith('.mzml') or f.lower().endswith('.mzml')]

        if not raw_files and not mzml_files:
            messagebox.showwarning("No Valid Files", "Please add .RAW or .mzML files.")
            return

        # Disable button during run
        self.run_btn.config(state="disabled", text="Running...")
        self.progress.start(10)
        self.status_var.set("Processing...")

        # Run in background thread
        thread = threading.Thread(target=self._run_thread,
                                  args=(raw_files, mzml_files), daemon=True)
        thread.start()

    def _run_thread(self, raw_files, mzml_files):
        try:
            from pipeline import run_gcms_pipeline

            config = {
                'min_sn': float(self.sn_var.get()),
                'min_rmf': int(self.rmf_var.get()),
                'use_nist': self.nist_var.get(),
                'deconv_enabled': self.deconv_var.get(),
            }
            out_dir = self.out_var.get() or str(SCRIPT_DIR / "output")

            result = run_gcms_pipeline(
                raw_files=raw_files if raw_files else None,
                mzml_files=mzml_files if mzml_files else None,
                output_dir=out_dir,
                config=config,
            )

            self.root.after(0, self._on_complete, result)
        except Exception as e:
            import traceback
            self.root.after(0, self._on_error, str(e) + "\n" + traceback.format_exc())

    def _on_complete(self, result):
        self.run_btn.config(state="normal", text="▶  Run Pipeline")
        self.progress.stop()
        if result and result.get('output_file'):
            self.status_var.set(f"Done! Output: {result['output_file']}")
            # Ask to open output folder
            if messagebox.askyesno("Complete", f"Processing complete!\n\nOutput: {result['output_file']}\n\nOpen output folder?"):
                os.startfile(str(Path(result['output_file']).parent))
        else:
            self.status_var.set("Done (no results generated)")

    def _on_error(self, error_msg):
        self.run_btn.config(state="normal", text="▶  Run Pipeline")
        self.progress.stop()
        self.status_var.set("Error occurred!")
        print(f"\n[ERROR]\n{error_msg}")
        messagebox.showerror("Error", error_msg[:500])


class TextRedirector:
    """Redirect print/output to the GUI text widget."""
    def __init__(self, text_widget):
        self.text_widget = text_widget
        self._stdout = sys.stdout

    def write(self, s):
        self._stdout.write(s)
        self.text_widget.configure(state="normal")
        self.text_widget.insert(tk.END, s)
        self.text_widget.see(tk.END)
        self.text_widget.configure(state="disabled")

    def flush(self):
        self._stdout.flush()


def main():
    root = tk.Tk()
    app = GCMSApp(root)

    # Center window
    root.update_idletasks()
    w, h = root.winfo_width(), root.winfo_height()
    sw, sh = root.winfo_screenwidth(), root.winfo_screenheight()
    x, y = (sw - w) // 2, (sh - h) // 2
    root.geometry(f"+{x}+{y}")

    root.mainloop()


if __name__ == "__main__":
    main()
