#!/usr/bin/env python3
"""
Extract retention index from NIST MS Search RI database via DLL calls.

Uses ri64.dll (bundled with NIST MS Search) to read the RI database.
"""
import ctypes, os, sys, json
from ctypes import c_int, c_char_p, c_void_p, c_uint, POINTER, byref

NIST_DIR = r'C:\NIST26-EI-DEMO\MSSEARCH'
HERE = os.path.dirname(os.path.realpath(__file__))


def init_ri_dll():
    """Initialize the NIST RI DLL."""
    dll = ctypes.CDLL(os.path.join(NIST_DIR, 'ri64.dll'))

    # RI_OpenDB(const char *libpath) -> int
    dll.RI_OpenDB.argtypes = [c_char_p]
    dll.RI_OpenDB.restype = c_int

    # RI_GetRIRecords(int spec_loc, int col_type, void **records, int *n_records) -> int
    # spec_loc: compound ID in the main library
    # col_type: 0=non-polar, 1=polar, etc. (-1 for all)
    # records: pointer to array of RI_RECORD structures
    dll.RI_GetRIRecords.argtypes = [c_int, c_int, POINTER(c_void_p), POINTER(c_int)]
    dll.RI_GetRIRecords.restype = c_int

    # RI_CloseDB() -> int
    dll.RI_CloseDB.argtypes = []
    dll.RI_CloseDB.restype = c_int

    # RI_DumpDB(const char *outpath) -> int (exports entire database)
    dll.RI_DumpDB.argtypes = [c_char_p]
    dll.RI_DumpDB.restype = c_int

    return dll


def init_nistretn_dll():
    """Initialize the NIST retention DLL."""
    dll = ctypes.CDLL(os.path.join(NIST_DIR, 'nistretn64.dll'))

    # CheckRIFiles(const char *libpath) -> int
    dll.CheckRIFiles.argtypes = [c_char_p]
    dll.CheckRIFiles.restype = c_int

    # GetRIByCAS(const char *cas, int col_type) -> int
    dll.GetRIByCAS.argtypes = [c_char_p, c_int]
    dll.GetRIByCAS.restype = c_int

    # GetHomologID(int spec_loc) -> int
    dll.GetHomologID.argtypes = [c_int]
    dll.GetHomologID.restype = c_int

    return dll


def try_dump_database():
    """Try to dump the RI database to a text file."""
    print('Attempting RI_DumpDB...')
    dll = init_ri_dll()

    # Open DB first
    result = dll.RI_OpenDB(b'mainlib')
    print(f'RI_OpenDB(mainlib) = {result}')

    if result == 0:
        out_path = os.path.join(HERE, 'nist_ri_dump.txt')
        result2 = dll.RI_DumpDB(out_path.encode())
        print(f'RI_DumpDB = {result2}')
        if result2 == 0:
            if os.path.exists(out_path):
                size = os.path.getsize(out_path)
                print(f'Dump file created: {out_path} ({size:,d} bytes)')
                with open(out_path, 'r') as f:
                    print(f.read()[:2000])
            else:
                print('No output file created')
        dll.RI_CloseDB()
    return False


def try_get_ri_by_cas(cas_list):
    """Try to query RI for specific CAS numbers."""
    print(f'Testing GetRIByCAS for {len(cas_list)} compounds...')
    dll = init_nistretn_dll()

    # Initialize
    result = dll.CheckRIFiles(b'mainlib')
    print(f'CheckRIFiles = {result}')

    if result == 0:
        for cas in cas_list:
            ri = dll.GetRIByCAS(cas.encode(), -1)  # -1 = all column types
            print(f'  CAS={cas}: RI={ri}')
    return False


if __name__ == '__main__':
    print('Approach 1: RI_DumpDB (dump entire database)')
    try:
        try_dump_database()
    except Exception as e:
        print(f'DumpDB failed: {e}')

    print()
    print('Approach 2: GetRIByCAS (query by CAS)')
    try:
        try_get_ri_by_cas(['66-25-1', '124-19-6', '100-52-7', '124-13-0'])
    except Exception as e:
        print(f'GetRIByCAS failed: {e}')

    print()
    print('Note: If both fail with segfault, the DLLs need different initialization')
    print('or the function signatures need adjustment.')
    print()
    print('Alternative: Use NIST MS Search GUI to export RI data:')
    print('  1. Open NIST MS Search')
    print('  2. Tools -> Options -> Export')
    print('  3. Select RI column for export')
    print('  4. Export search results')
