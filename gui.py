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

class GCMSApp:
    def __init__(self, root):
        self.root = root
        self.root.title('GC-MS Auto-Processing v2.1')
        self.root.geometry('680x600')
        self.root.resizable(True, True)
        bg = '#f5f5f5'
        ac = '#2F5496'
        self.root.configure(bg=bg)
        tf = tk.Frame(root, bg=ac, height=60)
        tf.pack(fill='x')
        tf.pack_propagate(False)
        tk.Label(tf, text='GC-MS Auto-Processing System v2.1',
                font=('Segoe UI', 16, 'bold'), fg='white', bg=ac).pack(pady=(10, 0))
        tk.Label(tf, text='Copyright (c) 2026 go ho',
                font=('Segoe UI', 8), fg='#b0c4de', bg=ac).pack()

        ff = tk.LabelFrame(root, text='Input Files', font=('Segoe UI', 10, 'bold'),
                           bg=bg, padx=10, pady=10)
        ff.pack(fill='x', padx=20, pady=(10, 8))
        self.lb = tk.Listbox(ff, height=4, font=('Consolas', 9))
        self.lb.pack(fill='x', pady=(0, 5))
        bf = tk.Frame(ff, bg=bg)
        bf.pack(fill='x')
        tk.Button(bf, text='+ Add samples (.RAW / .mzML / .qgd)', command=self.add_sample,
                 bg='#4a90d9', fg='white', font=('Segoe UI', 9),
                 padx=10).pack(side='left', padx=(0, 5))
        tk.Button(bf, text='Clear', command=self.clear,
                 bg='#dc3545', fg='white', font=('Segoe UI', 9),
                 padx=10).pack(side='left')

        sf = tk.LabelFrame(root, text='Settings', font=('Segoe UI', 10, 'bold'),
                           bg=bg, padx=10, pady=10)
        sf.pack(fill='x', padx=20, pady=(0, 8))
        r1 = tk.Frame(sf, bg=bg)
        r1.pack(fill='x', pady=3)
        tk.Label(r1, text='Min S/N:', bg=bg, width=10).pack(side='left')
        self.sn = tk.StringVar(value='10')
        tk.Entry(r1, textvariable=self.sn, width=8).pack(side='left', padx=(0, 25))
        tk.Label(r1, text='Min RMF:', bg=bg, width=10).pack(side='left')
        self.rmf = tk.StringVar(value='700')
        tk.Entry(r1, textvariable=self.rmf, width=8).pack(side='left', padx=(0, 25))
        tk.Label(r1, text='Solvent delay (min):', bg=bg).pack(side='left')
        self.sd = tk.StringVar(value='4.0')
        tk.Entry(r1, textvariable=self.sd, width=6).pack(side='left')
        rs = tk.Frame(sf, bg=bg)
        rs.pack(fill='x', pady=3)
        tk.Label(rs, text='RI standard:', bg=bg, width=10).pack(side='left')
        self.std = tk.StringVar(value='')
        tk.Entry(rs, textvariable=self.std, width=42).pack(side='left', padx=(0, 5))
        tk.Button(rs, text='Browse', command=self.browse_std, bg='#e0e0e0', padx=8).pack(side='left')
        tk.Button(rs, text='X', command=lambda: self.std.set(''), bg='#e0e0e0', padx=4).pack(side='left', padx=(4, 0))
        r2 = tk.Frame(sf, bg=bg)
        r2.pack(fill='x', pady=3)
        tk.Label(r2, text='Output:', bg=bg, width=10).pack(side='left')
        self.out = tk.StringVar(value=str(SCRIPT_DIR / 'output'))
        tk.Entry(r2, textvariable=self.out, width=42).pack(side='left', padx=(0, 5))
        tk.Button(r2, text='Browse', command=self.browse,
                 bg='#e0e0e0', padx=8).pack(side='left')
        r3 = tk.Frame(sf, bg=bg)
        r3.pack(fill='x', pady=5)
        self.nist = tk.BooleanVar(value=True)
        tk.Checkbutton(r3, text='Use NIST engine (recommended)',
                      variable=self.nist, bg=bg).pack(side='left')
        self.dec = tk.BooleanVar(value=False)
        tk.Checkbutton(r3, text='Enable deconvolution',
                      variable=self.dec, bg=bg).pack(side='left', padx=(20, 0))

        self.btn = tk.Button(root, text='START PROCESSING', command=self.run,
                              bg=ac, fg='white', font=('Segoe UI', 13, 'bold'),
                              padx=40, pady=12)
        self.btn.pack(pady=(8, 5))
        self.st = tk.StringVar(value='Ready. Add files and click START.')
        tk.Label(root, textvariable=self.st, font=('Segoe UI', 9),
                fg='#666', bg=bg).pack(pady=(0, 3))
        self.pb = ttk.Progressbar(root, mode='indeterminate', length=450)
        self.pb.pack(pady=(0, 8))
        self.ot = tk.Text(root, height=8, font=('Consolas', 8),
                          bg='#1e1e1e', fg='#d4d4d4', state='disabled')
        self.ot.pack(fill='both', expand=True, padx=20, pady=(0, 5))
        tk.Label(root, text='Copyright (c) 2026 go ho',
                font=('Segoe UI', 7), fg='#aaa', bg=bg).pack(pady=(0, 8))
        sys.stdout = TextRedirector(self.ot)
        self.log('GC-MS Auto-Processing System v2.1')
        self.log('Copyright (c) 2026 go ho')
        self.log('Ready.')

    def log(self, msg):
        print(msg)

    def add_sample(self):
        files = filedialog.askopenfilenames(
            title='Select sample files (Thermo .RAW / .mzML / Shimadzu .qgd)',
            filetypes=[('GC-MS data', '*.raw *.RAW *.mzML *.mzml *.qgd *.QGD'),
                       ('All', '*.*')])
        for f in files:
            if f not in self.lb.get(0, tk.END):
                self.lb.insert(tk.END, f)

    def clear(self):
        self.lb.delete(0, tk.END)

    def browse(self):
        d = filedialog.askdirectory(title='Select output folder')
        if d:
            self.out.set(d)

    def browse_std(self):
        f = filedialog.askopenfilename(
            title='Select n-alkane standard for RI calibration (.qgd/.RAW/.mzML)',
            filetypes=[('GC-MS data', '*.raw *.RAW *.mzML *.mzml *.qgd *.QGD'),
                       ('All', '*.*')])
        if f:
            self.std.set(f)

    def run(self):
        files = list(self.lb.get(0, tk.END))
        if not files:
            messagebox.showwarning('No Files', 'Please add sample files first (.RAW/.mzML/.qgd).')
            return
        # .RAW -> native Thermo reader; .mzML/.qgd -> load_sample (dispatched by extension)
        raw = [f for f in files if f.lower().endswith('.raw')]
        other = [f for f in files if not f.lower().endswith('.raw')]
        self.btn.config(state='disabled', text='Running...')
        self.pb.start(10)
        self.st.set('Processing ' + str(len(files)) + ' sample(s)...')
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
            out_dir = self.out.get() or str(SCRIPT_DIR / 'output')
            result = run_gcms_pipeline(
                raw_files=raw if raw else None,
                mzml_files=mzml if mzml else None,
                output_dir=out_dir,
                config=cfg,
            )
            self.root.after(0, self._done, result)
        except Exception as e:
            import traceback
            self.root.after(0, self._err, str(e) + '\n' + traceback.format_exc())

    def _done(self, result):
        self.btn.config(state='normal', text='START PROCESSING')
        self.pb.stop()
        if result and result.get('output_file'):
            o = result['output_file']
            outs = result.get('output_files', [o])
            self.st.set('Done: ' + str(len(outs)) + ' report(s) in ' + str(Path(o).parent))
            n = result.get('processed', '?')
            t = result.get('total', '?')
            print('\nDone! ' + str(n) + '/' + str(t) + ' samples')
            print('每个样品生成一个独立报告 (' + str(len(outs)) + ' 个 Excel):')
            for f in outs:
                print('  Output: ' + f)
            print('Excel sheets: 汇总 (summary) -> 待复核 (review checklist) -> Results')
            if messagebox.askyesno('Complete',
                                   'Done. ' + str(n) + '/' + str(t) + ' samples → '
                                   + str(len(outs)) + ' 个独立报告(每样品一个).\n\n'
                                   'Excel: 汇总 → 待复核(只看要审的) → Results\n\nOpen output folder?'):
                os.startfile(str(Path(o).parent))
        else:
            self.st.set('Done (no output)')

    def _err(self, msg):
        self.btn.config(state='normal', text='START PROCESSING')
        self.pb.stop()
        self.st.set('Error!')
        print('\n[ERROR]\n' + msg)
        messagebox.showerror('Error', msg[:500])


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
