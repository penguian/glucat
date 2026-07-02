import os
import re
import sys

base_dir = "/home/leopardi/sync/src/penguian/glucat/benchmarks/AMD-Ryzen-7-8840HS"
cfgs = [
    "armadillo", "armadillo-blas", "armadillo-blas-openmp", "armadillo-flexiblas", "armadillo-flexiblas-openmp",
    "armadillo-openblas", "armadillo-openblas-openmp", "armadillo-openmp",
    "eigen", "eigen-blas", "eigen-blas-openmp", "eigen-flexiblas", "eigen-flexiblas-openmp",
    "eigen-openblas", "eigen-openblas-openmp", "eigen-openmp"
]

print("Comparing gcc vs gcc.sign_helper benchmarks on AMD Ryzen 7 8840HS...")

# We will collect comparison results for:
# 1. products-8: framed_multi<double> Cl(8,8) * under Fill: 0.5
# 2. products-8: matrix_multi<double> Cl(8,8) * under Fill: 0.5
# 3. squaring-11: framed_multi<double> Cl(11,11) squaring * under Fill: 0.5
# 4. squaring-11: matrix_multi<double> Cl(11,11) squaring * under Fill: 0.5
# 5. transforms-8: GFFT old vs new at Cl(8,8) for mm and fm under Fill: 0.5

results = {}

for cfg in cfgs:
    results[cfg] = {}
    for ver in ["gcc", "gcc.sign_helper"]:
        results[cfg][ver] = {}
        
        # 1 & 2: products-8.out
        prod_path = os.path.join(base_dir, ver, cfg, "products-8.out")
        if os.path.exists(prod_path):
            with open(prod_path, 'r') as f:
                content = f.read()
            
            # framed_multi<double>
            fm_idx = content.find("framed_multi<double>")
            if fm_idx != -1:
                fm_sub = content[fm_idx:]
                m = re.search(r"Cl\(\s*8,\s*8\) in Cl\(\s*8,\s*8\) CPU =\s*[\d\.]+ ms \(setup\)\s*([\d\.]+) \(\*\)", fm_sub)
                if m:
                    results[cfg][ver]["prod_fm_mul"] = float(m.group(1))
            
            # matrix_multi<double>
            mm_idx = content.find("matrix_multi<double>")
            if mm_idx != -1:
                mm_sub = content[mm_idx:]
                m = re.search(r"Cl\(\s*8,\s*8\) in Cl\(\s*8,\s*8\) CPU =\s*[\d\.]+ ms \(setup\)\s*([\d\.]+) \(\*\)", mm_sub)
                if m:
                    results[cfg][ver]["prod_mm_mul"] = float(m.group(1))

        # 3 & 4: squaring-11.out
        sq_path = os.path.join(base_dir, ver, cfg, "squaring-11.out")
        if os.path.exists(sq_path):
            with open(sq_path, 'r') as f:
                content = f.read()
            
            # framed_multi<double>
            fm_idx = content.find("framed_multi<double>")
            if fm_idx != -1:
                fm_sub = content[fm_idx:]
                m = re.search(r"Cl\(\s*11,\s*11\) in Cl\(\s*11,\s*11\) CPU =\s*[\d\.]+ ms \(setup\)\s*([\d\.]+) \(\*\)", fm_sub)
                if m:
                    results[cfg][ver]["sq_fm_mul"] = float(m.group(1))
            
            # matrix_multi<double>
            mm_idx = content.find("matrix_multi<double>")
            if mm_idx != -1:
                mm_sub = content[mm_idx:]
                m = re.search(r"Cl\(\s*11,\s*11\) in Cl\(\s*11,\s*11\) CPU =\s*[\d\.]+ ms \(setup\)\s*([\d\.]+) \(\*\)", mm_sub)
                if m:
                    results[cfg][ver]["sq_mm_mul"] = float(m.group(1))

        # 5: transforms-8.out
        trans_path = os.path.join(base_dir, ver, cfg, "transforms-8.out")
        if os.path.exists(trans_path):
            with open(trans_path, 'r') as f:
                content = f.read()
            # find Cl(8,8) line under Fill: 0.5
            m = re.search(r"Cl\(\s*8,\s*8\) in Cl\(\s*8,\s*8\) CPU = mm:\s*([\d\.]+)\s*\(old\)\s*([\d\.]+)\s*\(new\) fm:\s*([\d\.]+)\s*\(old\)\s*([\d\.]+)\s*\(new\)", content)
            if m:
                results[cfg][ver]["trans_mm_old"] = float(m.group(1))
                results[cfg][ver]["trans_mm_new"] = float(m.group(2))
                results[cfg][ver]["trans_fm_old"] = float(m.group(3))
                results[cfg][ver]["trans_fm_new"] = float(m.group(4))

# Print products-8 comparison
print("\n### 1. products-8.out - framed_multi<double> vs matrix_multi<double> Multiplication (Cl(8,8) * at Fill: 0.5)")
print("| Configuration | Backend | gcc framed_multi (ms) | gcc.sign_helper framed_multi (ms) | Ratio | gcc matrix_multi (ms) | gcc.sign_helper matrix_multi (ms) | Ratio |")
print("| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: |")
for cfg in cfgs:
    fm_gcc = results[cfg]["gcc"].get("prod_fm_mul", float('nan'))
    fm_sh = results[cfg]["gcc.sign_helper"].get("prod_fm_mul", float('nan'))
    fm_ratio = fm_sh / fm_gcc if fm_gcc > 0 else float('nan')
    
    mm_gcc = results[cfg]["gcc"].get("prod_mm_mul", float('nan'))
    mm_sh = results[cfg]["gcc.sign_helper"].get("prod_mm_mul", float('nan'))
    mm_ratio = mm_sh / mm_gcc if mm_gcc > 0 else float('nan')
    
    print(f"| {cfg} | {'Eigen' if 'eigen' in cfg else 'Armadillo'} | {fm_gcc:.3f} | {fm_sh:.3f} | {fm_ratio:.2f}x | {mm_gcc:.3f} | {mm_sh:.3f} | {mm_ratio:.2f}x |")

# Print squaring-11 comparison
print("\n### 2. squaring-11.out - Squaring Performance (Cl(11,11) squaring * at Fill: 0.5)")
print("| Configuration | Backend | gcc framed_multi (ms) | gcc.sign_helper framed_multi (ms) | Ratio | gcc matrix_multi (ms) | gcc.sign_helper matrix_multi (ms) | Ratio |")
print("| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: |")
for cfg in cfgs:
    fm_gcc = results[cfg]["gcc"].get("sq_fm_mul", float('nan'))
    fm_sh = results[cfg]["gcc.sign_helper"].get("sq_fm_mul", float('nan'))
    fm_ratio = fm_sh / fm_gcc if fm_gcc > 0 else float('nan')
    
    mm_gcc = results[cfg]["gcc"].get("sq_mm_mul", float('nan'))
    mm_sh = results[cfg]["gcc.sign_helper"].get("sq_mm_mul", float('nan'))
    mm_ratio = mm_sh / mm_gcc if mm_gcc > 0 else float('nan')
    
    print(f"| {cfg} | {'Eigen' if 'eigen' in cfg else 'Armadillo'} | {fm_gcc:.1f} | {fm_sh:.1f} | {fm_ratio:.2f}x | {mm_gcc:.1f} | {mm_sh:.1f} | {mm_ratio:.2f}x |")

# Print transforms-8 comparison
print("\n### 3. transforms-8.out - Transform Performance (Cl(8,8) GFFT at Fill: 0.5)")
print("| Configuration | gcc mm_old (ms) | gcc.sign_helper mm_old (ms) | Ratio | gcc mm_new (ms) | gcc.sign_helper mm_new (ms) | Ratio | gcc fm_old (ms) | gcc.sign_helper fm_old (ms) | Ratio | gcc fm_new (ms) | gcc.sign_helper fm_new (ms) | Ratio |")
print("| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |")
for cfg in cfgs:
    t = results[cfg]
    m_o_g = t["gcc"].get("trans_mm_old", float('nan'))
    m_o_s = t["gcc.sign_helper"].get("trans_mm_old", float('nan'))
    m_o_r = m_o_s / m_o_g if m_o_g > 0 else float('nan')
    
    m_n_g = t["gcc"].get("trans_mm_new", float('nan'))
    m_n_s = t["gcc.sign_helper"].get("trans_mm_new", float('nan'))
    m_n_r = m_n_s / m_n_g if m_n_g > 0 else float('nan')
    
    f_o_g = t["gcc"].get("trans_fm_old", float('nan'))
    f_o_s = t["gcc.sign_helper"].get("trans_fm_old", float('nan'))
    f_o_r = f_o_s / f_o_g if f_o_g > 0 else float('nan')
    
    f_n_g = t["gcc"].get("trans_fm_new", float('nan'))
    f_n_s = t["gcc.sign_helper"].get("trans_fm_new", float('nan'))
    f_n_r = f_n_s / f_n_g if f_n_g > 0 else float('nan')
    
    print(f"| {cfg} | {m_o_g:.2f} | {m_o_s:.2f} | {m_o_r:.2f}x | {m_n_g:.2f} | {m_n_s:.2f} | {m_n_r:.2f}x | {f_o_g:.2f} | {f_o_s:.2f} | {f_o_r:.2f}x | {f_n_g:.2f} | {f_n_s:.2f} | {f_n_r:.2f}x |")
