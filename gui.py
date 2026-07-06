#!/usr/bin/env python3
import sys, os, threading, traceback, datetime
from pathlib import Path

# Global error logging
LOG_FILE = Path(__file__).parent / 'error.log'

def log_error(msg):
    with open(LOG_FILE, 'a', encoding='utf-8') as f:
        f.write(f'[{datetime.datetime.now()}] {msg}\n')

try:
    SCRIPT_DIR = Path(__file__).parent
    sys.path.insert(0, str(SCRIPT_DIR))
    import tkinter as tk
    from tkinter import filedialog, ttk, messagebox
except Exception as e:
    log_error(f'Import error: {e}\n{traceback.format_exc()}')
    raise

# ---- palette (matches the project's teal/ink identity) ----
BG   = '#eef2f5'   # page
CARD = '#ffffff'   # card surface
INK  = '#14293d'   # primary text
MUT  = '#5f7180'   # secondary text
AC   = '#0e7c66'   # accent (teal)
ACH  = '#0b6353'   # accent hover
ACD  = '#9cbcb2'   # accent disabled
LINE = '#dde5ea'   # hairline border
CON_BG = '#10202b'  # console bg
CON_FG = '#cbd8e0'  # console fg
UIFONT = 'Microsoft YaHei UI'
MONO = 'Consolas'


class GCMSApp:
    def __init__(self, root):
        self.root = root
        self.root.title('GC-MS 自动鉴定分析工具')
        self.root.geometry('760x880')
        self.root.minsize(720, 800)
        self.root.configure(bg=BG)
        self._setup_style()

        wrap = tk.Frame(root, bg=BG)
        wrap.pack(fill='both', expand=True, padx=24, pady=(20, 16))

        # ---- header ----
        head = tk.Frame(wrap, bg=BG)
        head.pack(fill='x')
        bar = tk.Frame(head, bg=AC, width=5, height=46)
        bar.pack(side='left', padx=(0, 14)); bar.pack_propagate(False)
        htxt = tk.Frame(head, bg=BG)
        htxt.pack(side='left', anchor='w')
        tk.Label(htxt, text='GC-MS 自动鉴定分析工具', bg=BG, fg=INK,
                 font=(UIFONT, 17, 'bold')).pack(anchor='w')
        tk.Label(htxt, text='原始数据 → 峰检测 · NIST 搜索 · RI 标定 · 分级审核 Excel',
                 bg=BG, fg=MUT, font=(UIFONT, 9)).pack(anchor='w', pady=(2, 0))

        # ---- input files card ----
        self._section(wrap, '输入样品')
        ff = self._card(wrap)
        self.lb = tk.Listbox(ff, height=4, font=(MONO, 9), bd=0, relief='flat',
                             bg='#f7fafb', fg=INK, highlightthickness=1,
                             highlightbackground=LINE, highlightcolor=AC,
                             selectbackground='#d7ede6', selectforeground=INK,
                             activestyle='none')
        self.lb.pack(fill='x', pady=(0, 10), ipady=3)
        bf = tk.Frame(ff, bg=CARD)
        bf.pack(fill='x')
        ttk.Button(bf, text='＋  添加样品  (.RAW / .mzML / .qgd)', style='Accent.TButton',
                   command=self.add_sample).pack(side='left')
        ttk.Button(bf, text='＋  安捷伦 .D 文件夹', style='Accent.TButton',
                   command=self.add_folder).pack(side='left', padx=(8, 0))
        ttk.Button(bf, text='清空', style='Ghost.TButton',
                   command=self.clear).pack(side='left', padx=(8, 0))

        # ---- settings card ----
        self._section(wrap, '参数设置')
        sf = self._card(wrap)

        r1 = tk.Frame(sf, bg=CARD); r1.pack(fill='x', pady=(0, 8))
        self._field_label(r1, '最低信噪比')
        self.sn = tk.StringVar(value='10')
        ttk.Entry(r1, textvariable=self.sn, width=7, style='F.TEntry').pack(side='left', padx=(0, 22))
        self._field_label(r1, '最低匹配因子')
        self.rmf = tk.StringVar(value='700')
        ttk.Entry(r1, textvariable=self.rmf, width=7, style='F.TEntry').pack(side='left', padx=(0, 22))
        self._field_label(r1, '溶剂延迟 (min)')
        self.sd = tk.StringVar(value='4.0')
        ttk.Entry(r1, textvariable=self.sd, width=7, style='F.TEntry').pack(side='left')

        rs = tk.Frame(sf, bg=CARD); rs.pack(fill='x', pady=8)
        self._field_label(rs, 'RI 标品', w=11)
        self.std = tk.StringVar(value='')
        ttk.Entry(rs, textvariable=self.std, style='F.TEntry').pack(side='left', fill='x', expand=True, padx=(0, 8))
        ttk.Button(rs, text='浏览', style='Ghost.TButton', command=self.browse_std).pack(side='left')
        ttk.Button(rs, text='✕', style='Ghost.TButton', width=3,
                   command=lambda: self.std.set('')).pack(side='left', padx=(6, 0))

        r2 = tk.Frame(sf, bg=CARD); r2.pack(fill='x', pady=8)
        self._field_label(r2, '输出目录', w=11)
        self.out = tk.StringVar(value=str(SCRIPT_DIR / 'output'))
        ttk.Entry(r2, textvariable=self.out, style='F.TEntry').pack(side='left', fill='x', expand=True, padx=(0, 8))
        ttk.Button(r2, text='浏览', style='Ghost.TButton', command=self.browse).pack(side='left')

        r3 = tk.Frame(sf, bg=CARD); r3.pack(fill='x', pady=(8, 0))
        self.nist = tk.BooleanVar(value=True)
        ttk.Checkbutton(r3, text='使用 NIST 官方引擎（推荐）', variable=self.nist,
                        style='Card.TCheckbutton').pack(side='left')
        self.dec = tk.BooleanVar(value=False)
        ttk.Checkbutton(r3, text='启用解卷积', variable=self.dec,
                        style='Card.TCheckbutton').pack(side='left', padx=(24, 0))

        r4 = tk.Frame(sf, bg=CARD); r4.pack(fill='x', pady=(10, 0))
        self.is_on = tk.BooleanVar(value=False)
        ttk.Checkbutton(r4, text='内标定向提取', variable=self.is_on,
                        style='Card.TCheckbutton').pack(side='left', padx=(0, 12))
        self._field_label(r4, 'm/z')
        self.is_mz = tk.StringVar(value='45')
        ttk.Entry(r4, textvariable=self.is_mz, width=5, style='F.TEntry').pack(side='left', padx=(0, 12))
        self._field_label(r4, 'RT')
        self.is_rt1 = tk.StringVar(value='11.8')
        ttk.Entry(r4, textvariable=self.is_rt1, width=5, style='F.TEntry').pack(side='left')
        tk.Label(r4, text='–', bg=CARD, fg=MUT).pack(side='left', padx=5)
        self.is_rt2 = tk.StringVar(value='13.2')
        ttk.Entry(r4, textvariable=self.is_rt2, width=5, style='F.TEntry').pack(side='left', padx=(0, 12))
        self._field_label(r4, '确认离子')
        self.is_conf = tk.StringVar(value='55')
        ttk.Entry(r4, textvariable=self.is_conf, width=5, style='F.TEntry').pack(side='left', padx=(0, 12))
        self._field_label(r4, '内标名')
        self.is_name = tk.StringVar(value='2-octanol')
        ttk.Entry(r4, textvariable=self.is_name, width=11, style='F.TEntry').pack(side='left')

        tk.Label(sf, text='内标定向提取：无需标品，按 m/z + RT 窗口把内标峰积分出来，'
                          '另存一张各样品内标响应表（用于归一/进样质控）。',
                 bg=CARD, fg=MUT, font=(UIFONT, 8), anchor='w',
                 wraplength=640, justify='left').pack(fill='x', pady=(6, 0))

        # ---- run button ----
        self.btn = ttk.Button(wrap, text='开始处理', style='Start.TButton', command=self.run)
        self.btn.pack(pady=(16, 8))

        self.st = tk.StringVar(value='就绪 · 添加样品后点击「开始处理」')
        tk.Label(wrap, textvariable=self.st, font=(UIFONT, 9), fg=MUT, bg=BG).pack()
        self.pb = ttk.Progressbar(wrap, mode='indeterminate', length=460, style='Teal.Horizontal.TProgressbar')
        self.pb.pack(pady=(6, 10))

        # ---- console ----
        self.ot = tk.Text(wrap, height=8, font=(MONO, 8), bd=0, relief='flat',
                          bg=CON_BG, fg=CON_FG, insertbackground=CON_FG,
                          highlightthickness=1, highlightbackground='#22323e',
                          padx=12, pady=10, state='disabled')
        self.ot.pack(fill='both', expand=True)
        tk.Label(wrap, text='Copyright (c) 2026 go ho', font=(UIFONT, 7),
                 fg='#aab6bf', bg=BG).pack(pady=(8, 0))

        sys.stdout = TextRedirector(self.ot)
        self.log('GC-MS 自动鉴定分析工具')
        self.log('Copyright (c) 2026 go ho')
        self.log('就绪。')

    # ---------- style / layout helpers ----------
    def _setup_style(self):
        s = ttk.Style()
        try:
            s.theme_use('clam')
        except tk.TclError:
            pass
        s.configure('Accent.TButton', background=AC, foreground='white',
                    font=(UIFONT, 10), borderwidth=0, relief='flat', padding=(16, 9))
        s.map('Accent.TButton', background=[('active', ACH), ('pressed', ACH)],
              foreground=[('disabled', '#eef2f5')])
        s.configure('Start.TButton', background=AC, foreground='white',
                    font=(UIFONT, 13, 'bold'), borderwidth=0, relief='flat', padding=(46, 13))
        s.map('Start.TButton', background=[('active', ACH), ('pressed', ACH), ('disabled', ACD)])
        s.configure('Ghost.TButton', background=BG, foreground=INK, font=(UIFONT, 10),
                    borderwidth=1, relief='flat', padding=(14, 8))
        s.map('Ghost.TButton', background=[('active', '#e2e8ec')],
              bordercolor=[('!active', LINE), ('active', AC)])
        s.configure('F.TEntry', fieldbackground='white', bordercolor=LINE,
                    lightcolor=LINE, darkcolor=LINE, borderwidth=1, padding=6, foreground=INK)
        s.map('F.TEntry', bordercolor=[('focus', AC)], lightcolor=[('focus', AC)])
        s.configure('Card.TCheckbutton', background=CARD, foreground=INK, font=(UIFONT, 10))
        s.map('Card.TCheckbutton', background=[('active', CARD)],
              indicatorcolor=[('selected', AC), ('!selected', 'white')])
        s.configure('Teal.Horizontal.TProgressbar', background=AC, troughcolor='#dde5ea',
                    borderwidth=0, thickness=6)

    def _section(self, parent, text):
        tk.Label(parent, text=text, bg=BG, fg=AC,
                 font=(UIFONT, 10, 'bold')).pack(anchor='w', pady=(16, 6))

    def _card(self, parent):
        c = tk.Frame(parent, bg=CARD, highlightthickness=1,
                     highlightbackground=LINE, highlightcolor=LINE)
        c.pack(fill='x')
        inner = tk.Frame(c, bg=CARD)
        inner.pack(fill='x', padx=16, pady=14)
        return inner

    def _field_label(self, parent, text, w=None):
        lbl = tk.Label(parent, text=text, bg=CARD, fg=MUT, font=(UIFONT, 10),
                       anchor='w')
        if w:
            lbl.config(width=w)
        lbl.pack(side='left', padx=(0, 10))

    def log(self, msg):
        print(msg)

    # ---------- actions (unchanged logic) ----------
    def add_sample(self):
        files = filedialog.askopenfilenames(
            title='选择样品文件 (Thermo .RAW / .mzML / 岛津 .qgd)',
            filetypes=[('GC-MS data', '*.raw *.RAW *.mzML *.mzml *.qgd *.QGD'),
                       ('All', '*.*')])
        for f in files:
            if f not in self.lb.get(0, tk.END):
                self.lb.insert(tk.END, f)

    def add_folder(self):
        d = filedialog.askdirectory(
            title='选择安捷伦 .D 文件夹（或含多个 .D 的父文件夹）')
        if not d:
            return
        existing = set(self.lb.get(0, tk.END))
        # A .D acquisition folder contains AcqData/MSScan.bin.
        def is_dot_d(p):
            return p.lower().rstrip('/\\').endswith('.d') and \
                os.path.isfile(os.path.join(p, 'AcqData', 'MSScan.bin'))
        if is_dot_d(d):
            targets = [d]
        else:                                    # a parent folder: collect all .D inside
            targets = sorted(os.path.join(d, x) for x in os.listdir(d)
                             if is_dot_d(os.path.join(d, x)))
        if not targets:
            messagebox.showwarning('未找到 .D 数据',
                                   '该文件夹不是安捷伦 .D，也没有包含任何 .D 样品。')
            return
        for t in targets:
            if t not in existing:
                self.lb.insert(tk.END, t)
        if len(targets) > 1:
            self.st.set('已添加 ' + str(len(targets)) + ' 个 .D 样品')

    def clear(self):
        self.lb.delete(0, tk.END)

    def browse(self):
        d = filedialog.askdirectory(title='选择输出目录')
        if d:
            self.out.set(d)

    def browse_std(self):
        f = filedialog.askopenfilename(
            title='选择正构烷烃 RI 标品 (.qgd / .RAW / .mzML)',
            filetypes=[('GC-MS data', '*.raw *.RAW *.mzML *.mzml *.qgd *.QGD'),
                       ('All', '*.*')])
        if f:
            self.std.set(f)

    def run(self):
        files = list(self.lb.get(0, tk.END))
        if not files:
            messagebox.showwarning('未选择文件', '请先添加样品文件（.RAW / .mzML / .qgd）。')
            return
        # .RAW -> native Thermo reader; .mzML/.qgd -> load_sample (dispatched by extension)
        raw = [f for f in files if f.lower().endswith('.raw')]
        other = [f for f in files if not f.lower().endswith('.raw')]
        self.btn.config(state='disabled', text='处理中…')
        self.pb.start(10)
        self.st.set('正在处理 ' + str(len(files)) + ' 个样品…')
        threading.Thread(target=self._run, args=(raw, other), daemon=True).start()

    def _run(self, raw, mzml):
        try:
            from pipeline import run_gcms_pipeline
            cfg = {
                'min_sn': float(self.sn.get()),
                'min_rmf': int(self.rmf.get()),
                'use_nist': self.nist.get(),
                'deconv_enabled': self.dec.get(),
            }
            try:
                cfg['solvent_delay'] = float(self.sd.get())
            except (ValueError, TypeError):
                pass
            std = self.std.get().strip()
            if std:
                cfg['standard_file'] = std
            if self.is_on.get():                 # internal-standard normalization column
                try:
                    cfg['is_ion'] = float(self.is_mz.get())
                    cfg['is_rt_min'] = float(self.is_rt1.get())
                    cfg['is_rt_max'] = float(self.is_rt2.get())
                    cfg['is_name'] = self.is_name.get().strip() or 'IS'
                except (ValueError, TypeError):
                    pass
            out_dir = self.out.get() or str(SCRIPT_DIR / 'output')
            result = run_gcms_pipeline(
                raw_files=raw if raw else None,
                mzml_files=mzml if mzml else None,
                output_dir=out_dir,
                config=cfg,
            )
            # optional: targeted internal-standard extraction (no calib standard needed)
            if self.is_on.get():
                try:
                    from is_extract import extract_internal_standard
                    nm = (self.is_name.get().strip() or 'IS')
                    print('\n[内标] 定向提取 ' + nm + ' (m/z ' + self.is_mz.get() + ') …')
                    df = extract_internal_standard(
                        list(raw) + list(mzml),
                        mz=float(self.is_mz.get()),
                        rt_min=float(self.is_rt1.get()), rt_max=float(self.is_rt2.get()),
                        confirm=float(self.is_conf.get()), name=nm, log=print)
                    is_out = os.path.join(out_dir, 'internal_standard_' + nm + '.xlsx')
                    df.to_excel(is_out, index=False)
                    print('[内标] 响应表已保存: ' + is_out)
                except Exception as ie:
                    print('[内标] 提取失败: ' + str(ie))
            self.root.after(0, self._done, result)
        except Exception as e:
            import traceback
            self.root.after(0, self._err, str(e) + '\n' + traceback.format_exc())

    def _done(self, result):
        self.btn.config(state='normal', text='开始处理')
        self.pb.stop()
        if result and result.get('output_file'):
            o = result['output_file']
            outs = result.get('output_files', [o])
            self.st.set('完成 · 生成 ' + str(len(outs)) + ' 个报告 · ' + str(Path(o).parent))
            n = result.get('processed', '?')
            t = result.get('total', '?')
            print('\nDone! ' + str(n) + '/' + str(t) + ' samples')
            print('每个样品生成一个独立报告 (' + str(len(outs)) + ' 个 Excel):')
            for f in outs:
                print('  Output: ' + f)
            print('Excel sheets: 汇总 (summary) -> 待复核 (review checklist) -> Results')
            if messagebox.askyesno('处理完成',
                                   '完成。' + str(n) + '/' + str(t) + ' 个样品 → '
                                   + str(len(outs)) + ' 个独立报告（每样品一个）。\n\n'
                                   'Excel：汇总 → 待复核（只看要审的）→ Results\n\n打开输出文件夹？'):
                os.startfile(str(Path(o).parent))
        else:
            self.st.set('完成（无输出）')

    def _err(self, msg):
        self.btn.config(state='normal', text='开始处理')
        self.pb.stop()
        self.st.set('出错了！')
        print('\n[ERROR]\n' + msg)
        messagebox.showerror('错误', msg[:500])


class TextRedirector:
    def __init__(self, widget):
        self.widget = widget
        self.stdout = sys.stdout
    def write(self, s):
        self.stdout.write(s)
        self.widget.configure(state='normal')
        self.widget.insert(tk.END, s)
        self.widget.see(tk.END)
        self.widget.configure(state='disabled')
    def flush(self):
        self.stdout.flush()


if __name__ == '__main__':
    try:
        root = tk.Tk()
        app = GCMSApp(root)
        root.update_idletasks()
        w = root.winfo_width()
        h = root.winfo_height()
        sw = root.winfo_screenwidth()
        sh = root.winfo_screenheight()
        root.geometry('+' + str((sw - w) // 2) + '+' + str((sh - h) // 2))
        root.mainloop()
    except Exception as e:
        log_error(f'Fatal error: {e}\n{traceback.format_exc()}')
        # Try to show error in a GUI dialog
        try:
            import tkinter.messagebox as mb
            mb.showerror('Error', f'Failed to start:\n\n{e}')
        except:
            pass
        raise
