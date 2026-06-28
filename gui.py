#!/usr/bin/env python3
"""
GC-MS 自动数据处理系统 v2.1
============================
基于 NIST 官方引擎的 GC-MS 全自动鉴定软件。
一键完成：RAW → 峰检测 → NIST 搜索 → Excel 报告

版权所有 (c) 2026 苟昊
"""
import sys, os, threading, json
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))

import tkinter as tk
from tkinter import filedialog, ttk, messagebox


# ---- 软件信息 ----
APP_NAME = "GC-MS 自动数据处理系统"
APP_VERSION = "v2.1"
APP_AUTHOR = "苟昊"
APP_COPYRIGHT = f"© 2026 {APP_AUTHOR}"


class GCMSApp:
    def __init__(self, root):
        self.root = root
        self.root.title(f"{APP_NAME} {APP_VERSION}")
        self.root.geometry("680x540")
        self.root.resizable(True, True)

        # Style
        self.bg = "#f5f5f5"
        self.fg = "#333333"
        self.accent = "#2F5496"
        self.root.configure(bg=self.bg)

        # ---- 标题栏 ----
        title_frame = tk.Frame(root, bg=self.accent, height=60)
        title_frame.pack(fill="x")
        title_frame.pack_propagate(False)

        tk.Label(title_frame, text=APP_NAME,
                font=("微软雅黑", 16, "bold"), fg="white", bg=self.accent).pack(pady=(10, 0))
        tk.Label(title_frame, text=f"{APP_VERSION}  |  {APP_COPYRIGHT}",
                font=("微软雅黑", 8), fg="#b0c4de", bg=self.accent).pack()

        # ---- 文件选择 ----
        file_frame = tk.LabelFrame(root, text="输入文件", font=("微软雅黑", 10, "bold"),
                                    bg=self.bg, fg=self.fg, padx=10, pady=10)
        file_frame.pack(fill="x", padx=20, pady=(10, 8))

        self.file_listbox = tk.Listbox(file_frame, height=4, font=("Consolas", 9))
        self.file_listbox.pack(fill="x", pady=(0, 5))

        btn_frame = tk.Frame(file_frame, bg=self.bg)
        btn_frame.pack(fill="x")
        tk.Button(btn_frame, text="＋ 添加 .RAW 文件", command=self.add_raw_files,
                 bg="#4a90d9", fg="white", font=("微软雅黑", 9), padx=10).pack(side="left", padx=(0, 5))
        tk.Button(btn_frame, text="＋ 添加 .mzML 文件", command=self.add_mzml_files,
                 bg="#6c757d", fg="white", font=("微软雅黑", 9), padx=10).pack(side="left", padx=(0, 5))
        tk.Button(btn_frame, text="清空列表", command=self.clear_files,
                 bg="#dc3545", fg="white", font=("微软雅黑", 9), padx=10).pack(side="left")

        # ---- 参数设置 ----
        settings_frame = tk.LabelFrame(root, text="参数设置", font=("微软雅黑", 10, "bold"),
                                        bg=self.bg, fg=self.fg, padx=10, pady=10)
        settings_frame.pack(fill="x", padx=20, pady=(0, 8))

        row1 = tk.Frame(settings_frame, bg=self.bg)
        row1.pack(fill="x", pady=3)
        tk.Label(row1, text="最低信噪比 (S/N):", bg=self.bg, width=16, anchor="w",
                font=("微软雅黑", 9)).pack(side="left")
        self.sn_var = tk.StringVar(value="10")
        tk.Entry(row1, textvariable=self.sn_var, width=8, font=("微软雅黑", 9)).pack(side="left", padx=(0, 25))
        tk.Label(row1, text="最低匹配因子 (RMF):", bg=self.bg, width=18, anchor="w",
                font=("微软雅黑", 9)).pack(side="left")
        self.rmf_var = tk.StringVar(value="700")
        tk.Entry(row1, textvariable=self.rmf_var, width=8, font=("微软雅黑", 9)).pack(side="left")

        row2 = tk.Frame(settings_frame, bg=self.bg)
        row2.pack(fill="x", pady=3)
        tk.Label(row2, text="输出目录:", bg=self.bg, width=16, anchor="w",
                font=("微软雅黑", 9)).pack(side="left")
        self.out_var = tk.StringVar(value=str(SCRIPT_DIR / "output"))
        tk.Entry(row2, textvariable=self.out_var, width=42, font=("微软雅黑", 9)).pack(side="left", padx=(0, 5))
        tk.Button(row2, text="浏览", command=self.browse_output,
                 bg="#e0e0e0", font=("微软雅黑", 8), padx=8).pack(side="left")

        row3 = tk.Frame(settings_frame, bg=self.bg)
        row3.pack(fill="x", pady=5)
        self.nist_var = tk.BooleanVar(value=True)
        tk.Checkbutton(row3, text="使用 NIST 官方引擎搜索（推荐，Chromeleon 级别准确度）",
                      variable=self.nist_var, bg=self.bg, anchor="w",
                      font=("微软雅黑", 9)).pack(side="left")
        self.deconv_var = tk.BooleanVar(value=False)
        tk.Checkbutton(row3, text="启用解卷积（共流出峰较多时）",
                      variable=self.deconv_var, bg=self.bg, anchor="w",
                      font=("微软雅黑", 9)).pack(side="left", padx=(20, 0))

        # ---- 运行按钮 ----
        self.run_btn = tk.Button(root, text="▶  开始处理", command=self.run_pipeline,
                                  bg="#2F5496", fg="white", font=("微软雅黑", 13, "bold"),
                                  padx=40, pady=12, cursor="hand2")
        self.run_btn.pack(pady=(8, 5))

        # ---- 状态 ----
        self.status_var = tk.StringVar(value="就绪 — 请添加文件后点击"开始处理"")
        tk.Label(root, textvariable=self.status_var,
                font=("微软雅黑", 9), fg="#666", bg=self.bg).pack(pady=(0, 3))

        self.progress = ttk.Progressbar(root, mode="indeterminate", length=450)
        self.progress.pack(pady=(0, 8))

        # ---- 输出窗口 ----
        out_label = tk.Label(root, text="运行日志", font=("微软雅黑", 9, "bold"), fg="#666", bg=self.bg)
        out_label.pack(anchor="w", padx=20)

        self.output_text = tk.Text(root, height=8, font=("Consolas", 8),
                                    bg="#1e1e1e", fg="#d4d4d4", state="disabled")
        self.output_text.pack(fill="both", expand=True, padx=20, pady=(0, 5))

        # 底部版权
        tk.Label(root, text=APP_COPYRIGHT, font=("微软雅黑", 7), fg="#aaa", bg=self.bg).pack(pady=(0, 8))

        # 重定向 stdout
        sys.stdout = TextRedirector(self.output_text)

        # 初始化就绪提示
        print(f">>> {APP_NAME} {APP_VERSION}")
        print(f">>> {APP_COPYRIGHT}")
        print(">>> 就绪，请添加文件后开始处理\n")

    def add_raw_files(self):
        files = filedialog.askopenfilenames(
            title="选择 .RAW 文件",
            filetypes=[("Thermo RAW 文件", "*.raw *.RAW"), ("所有文件", "*.*")]
        )
        for f in files:
            if f not in self.file_listbox.get(0, tk.END):
                self.file_listbox.insert(tk.END, f)

    def add_mzml_files(self):
        files = filedialog.askopenfilenames(
            title="选择 .mzML 文件",
            filetypes=[("mzML 文件", "*.mzML *.mzml"), ("所有文件", "*.*")]
        )
        for f in files:
            if f not in self.file_listbox.get(0, tk.END):
                self.file_listbox.insert(tk.END, f)

    def clear_files(self):
        self.file_listbox.delete(0, tk.END)

    def browse_output(self):
        d = filedialog.askdirectory(title="选择输出目录")
        if d:
            self.out_var.set(d)

    def run_pipeline(self):
        files = list(self.file_listbox.get(0, tk.END))
        if not files:
            messagebox.showwarning("未选择文件", "请先添加至少一个 .RAW 或 .mzML 文件。")
            return

        raw_files = [f for f in files if f.lower().endswith('.raw')]
        mzml_files = [f for f in files if not f.lower().endswith('.raw')]

        if not raw_files and not mzml_files:
            messagebox.showwarning("文件格式错误", "请添加 .RAW 或 .mzML 文件。")
            return

        self.run_btn.config(state="disabled", text="处理中，请稍候...")
        self.progress.start(10)
        self.status_var.set(f"正在处理 {len(files)} 个样品...")

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
        self.run_btn.config(state="normal", text="▶  开始处理")
        self.progress.stop()
        if result and result.get('output_file'):
            out = result['output_file']
            self.status_var.set(f"处理完成！输出：{out}")
            print(f"\n>>> 处理完成！")
            print(f">>> 输出文件：{out}")
            n = result.get('processed', 0)
            print(f">>> 成功处理 {n}/{result.get('total', '?')} 个样品\n")
            if messagebox.askyesno("处理完成",
                                   f"处理完成！\n\n成功：{result.get('processed', '?')}/{result.get('total', '?')} 个样品\n\n输出文件：{Path(out).name}\n\n打开输出文件夹？"):
                os.startfile(str(Path(out).parent))
        else:
            self.status_var.set("处理完成（未生成结果）")

    def _on_error(self, error_msg):
        self.run_btn.config(state="normal", text="▶  开始处理")
        self.progress.stop()
        self.status_var.set("处理出错！请查看日志")
        print(f"\n[错误]\n{error_msg}\n")
        messagebox.showerror("处理出错", error_msg[:500])


class TextRedirector:
    """将控制台输出重定向到 GUI 文本框。"""
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
    root.update_idletasks()
    w, h = root.winfo_width(), root.winfo_height()
    sw, sh = root.winfo_screenwidth(), root.winfo_screenheight()
    x, y = (sw - w) // 2, (sh - h) // 2
    root.geometry(f"+{x}+{y}")
    root.mainloop()


if __name__ == "__main__":
    main()
