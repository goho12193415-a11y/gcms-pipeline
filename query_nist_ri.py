#!/usr/bin/env python3
"""
NIST Chemistry WebBook — Retention Index query tool.

Queries NIST WebBook for literature RI values and caches results
locally in ri_database.json. Supports batch and on-demand lookup.

Usage:
    python query_nist_ri.py hexanal          # single compound
    python query_nist_ri.py --batch CAS list # from CAS file
    python query_nist_ri.py --rebuild        # rebuild RI DB from food_volatiles.msp
"""
import os, sys, re, json, time, urllib.request, urllib.error, ssl

HERE = os.path.dirname(os.path.realpath(__file__))
RI_DB_PATH = os.path.join(HERE, 'ri_database.json')

# Ignore SSL cert issues (some NIST pages have cert problems)
SSL_CTX = ssl.create_default_context()
SSL_CTX.check_hostname = False
SSL_CTX.verify_mode = ssl.CERT_NONE


def fetch_nist_page_by_cas(cas):
    """Fetch NIST Chemistry WebBook page for a given CAS number."""
    url = f'https://webbook.nist.gov/cgi/cbook.cgi?ID={cas}&Units=SI&cTG=on&cRI=on'
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'GCMS-Pipeline/1.0'})
        with urllib.request.urlopen(req, context=SSL_CTX, timeout=30) as resp:
            return resp.read().decode('utf-8', errors='replace')
    except Exception as e:
        print(f'  [WARN] Failed to fetch {url}: {e}')
        return None


def fetch_nist_page_by_name(name):
    """Fetch NIST page by compound name."""
    from urllib.parse import quote
    url = f'https://webbook.nist.gov/cgi/cbook.cgi?Name={quote(name)}&Units=SI&cTG=on&cRI=on'
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'GCMS-Pipeline/1.0'})
        with urllib.request.urlopen(req, context=SSL_CTX, timeout=30) as resp:
            return resp.read().decode('utf-8', errors='replace')
    except Exception as e:
        print(f'  [WARN] Failed to fetch by name {name}: {e}')
        return None


def parse_ri_from_html(html):
    """Parse retention index values from NIST WebBook HTML.

    Searches for the Gas Chromatography section and extracts:
    - Column type (e.g., DB-5, DB-WAX, HP-5, Carbowax)
    - RI value
    - Temperature program info

    Returns dict: {'DB-5': 1003, 'DB-WAX': 1290, ...}
    """
    if not html:
        return {}

    results = {}
    # Find the Gas Chromatography section
    gc_section_patterns = [
        r'<h2[^>]*>Gas Chromatography</h2>(.*?)(?:<h2|<hr)',
        r'Gas Chromatography.*?<table[^>]*>(.*?)</table>',
    ]

    # More robust: find all tables and look for RI-related headers
    tables = re.findall(r'<table[^>]*>(.*?)</table>', html, re.DOTALL | re.IGNORECASE)

    for table_html in tables:
        # Check if this table contains retention index data
        if 'retention index' not in table_html.lower() and 'kovats' not in table_html.lower():
            # Also check for "I" or "RI" column headers
            if not re.search(r'<th[^>]*>\s*[IR]I?\s*</th>', table_html):
                continue

        rows = re.findall(r'<tr[^>]*>(.*?)</tr>', table_html, re.DOTALL | re.IGNORECASE)
        current_column = ''
        for row in rows:
            cells = re.findall(r'<t[dh][^>]*>(.*?)</t[dh]>', row, re.DOTALL | re.IGNORECASE)
            clean_cells = [re.sub(r'<[^>]+>', '', c).strip() for c in cells]

            # Detect column type from row text
            row_text = ' '.join(clean_cells)
            column_keywords = {
                'DB-5': ['db-5', 'hp-5', 'se-54', 'cp sil 8', 'rtx-5', 'zb-5', 'bp-5', 'spb-5',
                         '007-5', 'ov-5', 'se-52', '5% phenyl'],
                'DB-WAX': ['db-wax', 'hp-20m', 'carbowax', 'cp wax', 'bp-20', 'zb-wax',
                          'innowax', 'hp-innowax', 'supelcowax', 'wax', 'peg', '20m',
                          'polyethylene glycol', 'ffap'],
                'DB-1': ['db-1', 'hp-1', 'se-30', 'ov-1', 'bp-1', 'rtx-1', 'zb-1',
                        '100% dimethyl', 'methyl silicone'],
                'DB-1701': ['db-1701', 'cp sil 19', 'ov-1701', 'bp-10'],
            }

            detected_col = None
            row_lower = row_text.lower()
            for col_name, keywords in column_keywords.items():
                for kw in keywords:
                    if kw in row_lower:
                        detected_col = col_name
                        break
                if detected_col:
                    break

            if detected_col:
                current_column = detected_col

            # Try to extract RI value from cells
            for cell in clean_cells:
                # Match patterns like "1234", "1234.5", "~1234"
                ri_match = re.match(r'^~?(\d{3,4})(?:[.]\d)?$', cell)
                if ri_match:
                    ri_val = int(ri_match.group(1))
                    # Sanity check: RI should be 300-4000 range
                    if 400 <= ri_val <= 4000:
                        col_name = current_column if current_column else 'Unknown'
                        if col_name not in results or ri_val < 2000:  # prefer lower values (standard temp)
                            results[col_name] = ri_val

    return results


def query_compound(name=None, cas=None):
    """Query NIST WebBook for a compound's RI values.

    Returns dict of {column_type: ri_value} or empty dict.
    """
    html = None
    if cas:
        html = fetch_nist_page_by_cas(cas)
    if not html and name:
        html = fetch_nist_page_by_name(name)

    if not html:
        return {}

    return parse_ri_from_html(html)


def load_ri_db():
    """Load existing RI database."""
    if not os.path.exists(RI_DB_PATH):
        return {'compounds': {}}
    with open(RI_DB_PATH, 'r', encoding='utf-8') as f:
        return json.load(f)


def save_ri_db(db):
    """Save RI database."""
    with open(RI_DB_PATH, 'w', encoding='utf-8') as f:
        json.dump(db, f, ensure_ascii=False, indent=2)
    print(f'Saved {len(db.get("compounds", {}))} entries to {RI_DB_PATH}')


def query_and_cache(name, ri_db, cas=None, force=False):
    """Query NIST for a compound and cache results. Rate-limited."""
    key = name.lower()
    existing = ri_db.get('compounds', {}).get(key, {})

    # Skip if already has RI data and not forced
    if not force and existing and any(v is not None for v in existing.values()):
        return existing

    print(f'  Querying: {name}...')
    ri_values = query_compound(name=name, cas=cas)

    if ri_values:
        ri_db.setdefault('compounds', {})[key] = ri_values
        print(f'    -> RI: {ri_values}')
    else:
        print(f'    -> No RI data found on NIST WebBook')
        ri_db.setdefault('compounds', {})[key] = {}

    return ri_values


def main():
    import argparse
    ap = argparse.ArgumentParser(description='Query NIST WebBook for retention indices')
    ap.add_argument('query', nargs='?', help='Compound name or CAS number')
    ap.add_argument('--batch', help='File with CAS numbers (one per line)')
    ap.add_argument('--rebuild', action='store_true',
                    help='Rebuild RI database from food_volatiles.msp')
    ap.add_argument('--force', action='store_true', help='Force re-query even if cached')
    ap.add_argument('--delay', type=float, default=3.0,
                    help='Delay between queries in seconds (default: 3.0)')
    args = ap.parse_args()

    ri_db = load_ri_db()
    print(f'RI database: {len(ri_db.get("compounds", {}))} existing entries')

    if args.rebuild:
        # Load compounds from food_volatiles.msp
        food_path = os.path.join(HERE, 'food_volatiles.msp')
        if not os.path.exists(food_path):
            print(f'food_volatiles.msp not found at {food_path}')
            return

        sys.path.insert(0, HERE)
        from gcms_pipeline import parse_msp
        compounds = parse_msp(food_path)
        print(f'Rebuilding RI database for {len(compounds)} compounds...')

        for i, comp in enumerate(compounds):
            name = comp['name']
            cas = comp.get('cas', '')
            print(f'[{i+1}/{len(compounds)}]', end=' ')
            query_and_cache(name, ri_db, cas=cas if cas else None, force=args.force)
            if i < len(compounds) - 1:
                time.sleep(args.delay)

        save_ri_db(ri_db)

    elif args.batch:
        with open(args.batch, 'r') as f:
            cas_list = [line.strip() for line in f if line.strip()]
        print(f'Processing {len(cas_list)} CAS numbers...')
        for i, cas in enumerate(cas_list):
            print(f'[{i+1}/{len(cas_list)}] CAS={cas}', end=' ')
            query_and_cache(cas, ri_db, cas=cas, force=args.force)
            if i < len(cas_list) - 1:
                time.sleep(args.delay)
        save_ri_db(ri_db)

    elif args.query:
        # Single query
        name = args.query
        query_and_cache(name, ri_db, cas=name if '-' in name else None, force=args.force)
        save_ri_db(ri_db)

    else:
        print("Usage: python query_nist_ri.py <compound_name>")
        print("       python query_nist_ri.py --rebuild")
        print("       python query_nist_ri.py --batch cas_list.txt")


if __name__ == '__main__':
    main()
