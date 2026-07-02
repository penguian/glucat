import os
import re

base_dir = "/home/leopardi/sync/src/penguian/glucat/benchmarks/AMD-Ryzen-7-8840HS"
cfgs = [
    "armadillo", "armadillo-blas", "armadillo-blas-openmp", "armadillo-flexiblas", "armadillo-flexiblas-openmp",
    "armadillo-openblas", "armadillo-openblas-openmp", "armadillo-openmp",
    "eigen", "eigen-blas", "eigen-blas-openmp", "eigen-flexiblas", "eigen-flexiblas-openmp",
    "eigen-openblas", "eigen-openblas-openmp", "eigen-openmp"
]

results = {}

for cfg in cfgs:
    results[cfg] = {}
    for ver in ["gcc", "gcc.sign_helper"]:
        results[cfg][ver] = {}
        cfg_dir = os.path.join(base_dir, ver, cfg)
        if not os.path.exists(cfg_dir):
            continue

        # 1. products-8.out
        prod_path = os.path.join(cfg_dir, "products-8.out")
        if os.path.exists(prod_path):
            with open(prod_path, 'r') as f:
                content = f.read()
            # framed_multi<double>
            fm_idx = content.find("framed_multi<double>")
            if fm_idx != -1:
                m = re.search(r"Cl\(\s*8,\s*8\) in Cl\(\s*8,\s*8\) CPU =\s*[\d\.]+ ms \(setup\)\s*([\d\.]+) \(\*\)", content[fm_idx:])
                if m:
                    results[cfg][ver]["prod_fm_mul"] = float(m.group(1))
            # matrix_multi<double>
            mm_idx = content.find("matrix_multi<double>")
            if mm_idx != -1:
                m = re.search(r"Cl\(\s*8,\s*8\) in Cl\(\s*8,\s*8\) CPU =\s*[\d\.]+ ms \(setup\)\s*([\d\.]+) \(\*\)", content[mm_idx:])
                if m:
                    results[cfg][ver]["prod_mm_mul"] = float(m.group(1))

        # 2. squaring-11.out
        sq_path = os.path.join(cfg_dir, "squaring-11.out")
        if os.path.exists(sq_path):
            with open(sq_path, 'r') as f:
                content = f.read()
            # framed_multi<double>
            fm_idx = content.find("framed_multi<double>")
            if fm_idx != -1:
                m = re.search(r"Cl\(\s*11,\s*11\) in Cl\(\s*11,\s*11\) CPU =\s*[\d\.]+ ms \(setup\)\s*([\d\.]+) \(\*\)", content[fm_idx:])
                if m:
                    results[cfg][ver]["sq_fm_mul"] = float(m.group(1))
            # matrix_multi<double>
            mm_idx = content.find("matrix_multi<double>")
            if mm_idx != -1:
                m = re.search(r"Cl\(\s*11,\s*11\) in Cl\(\s*11,\s*11\) CPU =\s*[\d\.]+ ms \(setup\)\s*([\d\.]+) \(\*\)", content[mm_idx:])
                if m:
                    results[cfg][ver]["sq_mm_mul"] = float(m.group(1))

        # 3. transforms-8.out
        trans_path = os.path.join(cfg_dir, "transforms-8.out")
        if os.path.exists(trans_path):
            with open(trans_path, 'r') as f:
                content = f.read()
            m = re.search(r"Cl\(\s*8,\s*8\) in Cl\(\s*8,\s*8\) CPU = mm:\s*([\d\.]+)\s*\(old\)\s*([\d\.]+)\s*\(new\) fm:\s*([\d\.]+)\s*\(old\)\s*([\d\.]+)\s*\(new\)", content)
            if m:
                results[cfg][ver]["trans_mm_old"] = float(m.group(1))
                results[cfg][ver]["trans_mm_new"] = float(m.group(2))
                results[cfg][ver]["trans_fm_old"] = float(m.group(3))
                results[cfg][ver]["trans_fm_new"] = float(m.group(4))

        # 4. gfft_test-11.out
        gfft_path = os.path.join(cfg_dir, "gfft_test-11.out")
        if os.path.exists(gfft_path):
            with open(gfft_path, 'r') as f:
                content = f.read()
            m = re.search(r"R_\{-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8,9,10,11\} in R_\{-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8,9,10,11\}:\s*\n\s*CPU = mm:\s*([\d\.]+)\s*fm:\s*([\d\.]+)", content)
            if m:
                results[cfg][ver]["gfft_mm"] = float(m.group(1))
                results[cfg][ver]["gfft_fm"] = float(m.group(2))

        # 5. versor-16.out
        versor_path = os.path.join(cfg_dir, "versor-16.out")
        if os.path.exists(versor_path):
            with open(versor_path, 'r') as f:
                content = f.read()
            # framed_multi
            fm_idx = content.find("--- Benchmarking framed_multi ---")
            if fm_idx != -1:
                m = re.search(r"Cl\(\s*16,\s*0\) \| naive :\s*([\d\.]+) ms \| operator\| :\s*([\d\.]+) ms \| versor :\s*([\d\.]+) ms \| versor_exp :\s*([\d\.]+) ms", content[fm_idx:])
                if m:
                    results[cfg][ver]["versor_fm_naive"] = float(m.group(1))
                    results[cfg][ver]["versor_fm_or"] = float(m.group(2))
                    results[cfg][ver]["versor_fm_ver"] = float(m.group(3))
                    results[cfg][ver]["versor_fm_ver_exp"] = float(m.group(4))
            # matrix_multi
            mm_idx = content.find("--- Benchmarking matrix_multi ---")
            if mm_idx != -1:
                m = re.search(r"Cl\(\s*16,\s*0\) \| naive :\s*([\d\.]+) ms \| operator\| :\s*([\d\.]+) ms \| versor :\s*([\d\.]+) ms \| versor_exp :\s*([\d\.]+) ms", content[mm_idx:])
                if m:
                    results[cfg][ver]["versor_mm_naive"] = float(m.group(1))
                    results[cfg][ver]["versor_mm_or"] = float(m.group(2))
                    results[cfg][ver]["versor_mm_ver"] = float(m.group(3))
                    results[cfg][ver]["versor_mm_ver_exp"] = float(m.group(4))

        # 6. expressions-8.out
        expr_path = os.path.join(cfg_dir, "expressions-8.out")
        if os.path.exists(expr_path):
            with open(expr_path, 'r') as f:
                content = f.read()
            # framed_multi<double>
            fm_idx = content.find("framed_multi<double>")
            if fm_idx != -1:
                m = re.search(r"Cl\(\s*8,\s*8\) in Cl\(\s*8,\s*8\) CPU =\s*[\d\.]+ ms \(setup\)\s*[\d\.]+ \(comm\)\s*([\d\.]+) \(pade\)\s*([\d\.]+) \(series\)\s*([\d\.]+) \(mix\)\s*([\d\.]+) \(add\)", content[fm_idx:])
                if m:
                    results[cfg][ver]["expr_fm_pade"] = float(m.group(1))
                    results[cfg][ver]["expr_fm_series"] = float(m.group(2))
                    results[cfg][ver]["expr_fm_mix"] = float(m.group(3))
                    results[cfg][ver]["expr_fm_add"] = float(m.group(4))
            # matrix_multi<double>
            mm_idx = content.find("matrix_multi<double>")
            if mm_idx != -1:
                m = re.search(r"Cl\(\s*8,\s*8\) in Cl\(\s*8,\s*8\) CPU =\s*[\d\.]+ ms \(setup\)\s*[\d\.]+ \(comm\)\s*([\d\.]+) \(pade\)\s*([\d\.]+) \(series\)\s*([\d\.]+) \(mix\)\s*([\d\.]+) \(add\)", content[mm_idx:])
                if m:
                    results[cfg][ver]["expr_mm_pade"] = float(m.group(1))
                    results[cfg][ver]["expr_mm_series"] = float(m.group(2))
                    results[cfg][ver]["expr_mm_mix"] = float(m.group(3))
                    results[cfg][ver]["expr_mm_add"] = float(m.group(4))

# Let's print out the tables
print("# Performance Comparison: GCC baseline vs GCC sign_helper on AMD Ryzen 7 8840HS")
print("This report evaluates the performance impact of the `sign_helper` refactoring across 16 different Eigen/Armadillo backend configurations.\n")

print("## 1. Multiplications (products-8.out, Cl(8,8) * under Fill: 0.5)")
print("| Configuration | Backend | gcc fm (ms) | gcc.sign_helper fm (ms) | Ratio fm | gcc mm (ms) | gcc.sign_helper mm (ms) | Ratio mm |")
print("| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: |")
for cfg in cfgs:
    fm_gcc = results[cfg]["gcc"].get("prod_fm_mul", float('nan'))
    fm_sh = results[cfg]["gcc.sign_helper"].get("prod_fm_mul", float('nan'))
    fm_ratio = fm_sh / fm_gcc if fm_gcc > 0 else float('nan')
    mm_gcc = results[cfg]["gcc"].get("prod_mm_mul", float('nan'))
    mm_sh = results[cfg]["gcc.sign_helper"].get("prod_mm_mul", float('nan'))
    mm_ratio = mm_sh / mm_gcc if mm_gcc > 0 else float('nan')
    print(f"| {cfg} | {'Eigen' if 'eigen' in cfg else 'Armadillo'} | {fm_gcc:.3f} | {fm_sh:.3f} | {fm_ratio:.2f}x | {mm_gcc:.3f} | {mm_sh:.3f} | {mm_ratio:.2f}x |")

print("\n## 2. Squaring (squaring-11.out, Cl(11,11) squaring * under Fill: 0.5)")
print("| Configuration | Backend | gcc fm (ms) | gcc.sign_helper fm (ms) | Ratio fm | gcc mm (ms) | gcc.sign_helper mm (ms) | Ratio mm |")
print("| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: |")
for cfg in cfgs:
    fm_gcc = results[cfg]["gcc"].get("sq_fm_mul", float('nan'))
    fm_sh = results[cfg]["gcc.sign_helper"].get("sq_fm_mul", float('nan'))
    fm_ratio = fm_sh / fm_gcc if fm_gcc > 0 else float('nan')
    mm_gcc = results[cfg]["gcc"].get("sq_mm_mul", float('nan'))
    mm_sh = results[cfg]["gcc.sign_helper"].get("sq_mm_mul", float('nan'))
    mm_ratio = mm_sh / mm_gcc if mm_gcc > 0 else float('nan')
    print(f"| {cfg} | {'Eigen' if 'eigen' in cfg else 'Armadillo'} | {fm_gcc:.1f} | {fm_sh:.1f} | {fm_ratio:.2f}x | {mm_gcc:.1f} | {mm_sh:.1f} | {mm_ratio:.2f}x |")

print("\n## 3. GFFT Transforms (transforms-8.out, Cl(8,8) GFFT under Fill: 0.5)")
print("| Configuration | gcc mm_old (ms) | gcc.sign_helper mm_old (ms) | Ratio mm_old | gcc mm_new (ms) | gcc.sign_helper mm_new (ms) | Ratio mm_new | gcc fm_old (ms) | gcc.sign_helper fm_old (ms) | Ratio fm_old | gcc fm_new (ms) | gcc.sign_helper fm_new (ms) | Ratio fm_new |")
print("| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |")
for cfg in cfgs:
    r = results[cfg]
    mo_g = r["gcc"].get("trans_mm_old", float('nan'))
    mo_s = r["gcc.sign_helper"].get("trans_mm_old", float('nan'))
    mo_r = mo_s / mo_g if mo_g > 0 else float('nan')
    mn_g = r["gcc"].get("trans_mm_new", float('nan'))
    mn_s = r["gcc.sign_helper"].get("trans_mm_new", float('nan'))
    mn_r = mn_s / mn_g if mn_g > 0 else float('nan')
    fo_g = r["gcc"].get("trans_fm_old", float('nan'))
    fo_s = r["gcc.sign_helper"].get("trans_fm_old", float('nan'))
    fo_r = fo_s / fo_g if fo_g > 0 else float('nan')
    fn_g = r["gcc"].get("trans_fm_new", float('nan'))
    fn_s = r["gcc.sign_helper"].get("trans_fm_new", float('nan'))
    fn_r = fn_s / fn_g if fn_g > 0 else float('nan')
    print(f"| {cfg} | {mo_g:.2f} | {mo_s:.2f} | {mo_r:.2f}x | {mn_g:.2f} | {mn_s:.2f} | {mn_r:.2f}x | {fo_g:.2f} | {fo_s:.2f} | {fo_r:.2f}x | {fn_g:.2f} | {fn_s:.2f} | {fn_r:.2f}x |")

print("\n## 4. GFFT Test (gfft_test-11.out, R_{-11,11} GFFT under Fill: 0.5)")
print("| Configuration | gcc mm (ms) | gcc.sign_helper mm (ms) | Ratio mm | gcc fm (ms) | gcc.sign_helper fm (ms) | Ratio fm |")
print("| :--- | :---: | :---: | :---: | :---: | :---: | :---: |")
for cfg in cfgs:
    mm_gcc = results[cfg]["gcc"].get("gfft_mm", float('nan'))
    mm_sh = results[cfg]["gcc.sign_helper"].get("gfft_mm", float('nan'))
    mm_ratio = mm_sh / mm_gcc if mm_gcc > 0 else float('nan')
    fm_gcc = results[cfg]["gcc"].get("gfft_fm", float('nan'))
    fm_sh = results[cfg]["gcc.sign_helper"].get("gfft_fm", float('nan'))
    fm_ratio = fm_sh / fm_gcc if fm_gcc > 0 else float('nan')
    print(f"| {cfg} | {mm_gcc:.1f} | {mm_sh:.1f} | {mm_ratio:.2f}x | {fm_gcc:.1f} | {fm_sh:.1f} | {fm_ratio:.2f}x |")

print("\n## 5. Versor Product Benchmarks (versor-16.out, Cl(16,0) under Fill: 0.5)")
print("### Framed Multivector (`framed_multi`) Versors")
print("| Configuration | gcc naive (ms) | gcc.sh naive (ms) | Ratio naive | gcc op| (ms) | gcc.sh op| (ms) | Ratio op| | gcc versor (ms) | gcc.sh versor (ms) | Ratio versor | gcc v_exp (ms) | gcc.sh v_exp (ms) | Ratio v_exp |")
print("| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |")
for cfg in cfgs:
    r = results[cfg]
    n_g = r["gcc"].get("versor_fm_naive", float('nan'))
    n_s = r["gcc.sign_helper"].get("versor_fm_naive", float('nan'))
    n_r = n_s / n_g if n_g > 0 else float('nan')
    o_g = r["gcc"].get("versor_fm_or", float('nan'))
    o_s = r["gcc.sign_helper"].get("versor_fm_or", float('nan'))
    o_r = o_s / o_g if o_g > 0 else float('nan')
    v_g = r["gcc"].get("versor_fm_ver", float('nan'))
    v_s = r["gcc.sign_helper"].get("versor_fm_ver", float('nan'))
    v_r = v_s / v_g if v_g > 0 else float('nan')
    ve_g = r["gcc"].get("versor_fm_ver_exp", float('nan'))
    ve_s = r["gcc.sign_helper"].get("versor_fm_ver_exp", float('nan'))
    ve_r = ve_s / ve_g if ve_g > 0 else float('nan')
    print(f"| {cfg} | {n_g:.1f} | {n_s:.1f} | {n_r:.2f}x | {o_g:.1f} | {o_s:.1f} | {o_r:.2f}x | {v_g:.1f} | {v_s:.1f} | {v_r:.2f}x | {ve_g:.1f} | {ve_s:.1f} | {ve_r:.2f}x |")

print("\n### Dense Matrix Multivector (`matrix_multi`) Versors")
print("| Configuration | gcc naive (ms) | gcc.sh naive (ms) | Ratio naive | gcc op| (ms) | gcc.sh op| (ms) | Ratio op| | gcc versor (ms) | gcc.sh versor (ms) | Ratio versor | gcc v_exp (ms) | gcc.sh v_exp (ms) | Ratio v_exp |")
print("| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |")
for cfg in cfgs:
    r = results[cfg]
    n_g = r["gcc"].get("versor_mm_naive", float('nan'))
    n_s = r["gcc.sign_helper"].get("versor_mm_naive", float('nan'))
    n_r = n_s / n_g if n_g > 0 else float('nan')
    o_g = r["gcc"].get("versor_mm_or", float('nan'))
    o_s = r["gcc.sign_helper"].get("versor_mm_or", float('nan'))
    o_r = o_s / o_g if o_g > 0 else float('nan')
    v_g = r["gcc"].get("versor_mm_ver", float('nan'))
    v_s = r["gcc.sign_helper"].get("versor_mm_ver", float('nan'))
    v_r = v_s / v_g if v_g > 0 else float('nan')
    ve_g = r["gcc"].get("versor_mm_ver_exp", float('nan'))
    ve_s = r["gcc.sign_helper"].get("versor_mm_ver_exp", float('nan'))
    ve_r = ve_s / ve_g if ve_g > 0 else float('nan')
    print(f"| {cfg} | {n_g:.1f} | {n_s:.1f} | {n_r:.2f}x | {o_g:.1f} | {o_s:.1f} | {o_r:.2f}x | {v_g:.1f} | {v_s:.1f} | {v_r:.2f}x | {ve_g:.1f} | {ve_s:.1f} | {ve_r:.2f}x |")

print("\n## 6. Algebra Expressions (expressions-8.out, Cl(8,8) under Fill: 0.5)")
print("### Framed Multivector (`framed_multi<double>`) Expressions")
print("| Configuration | gcc pade (ms) | gcc.sh pade (ms) | Ratio pade | gcc series (ms) | gcc.sh series (ms) | Ratio series | gcc mix (ms) | gcc.sh mix (ms) | Ratio mix | gcc add (ms) | gcc.sh add (ms) | Ratio add |")
print("| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |")
for cfg in cfgs:
    r = results[cfg]
    pa_g = r["gcc"].get("expr_fm_pade", float('nan'))
    pa_s = r["gcc.sign_helper"].get("expr_fm_pade", float('nan'))
    pa_r = pa_s / pa_g if pa_g > 0 else float('nan')
    se_g = r["gcc"].get("expr_fm_series", float('nan'))
    se_s = r["gcc.sign_helper"].get("expr_fm_series", float('nan'))
    se_r = se_s / se_g if se_g > 0 else float('nan')
    mi_g = r["gcc"].get("expr_fm_mix", float('nan'))
    mi_s = r["gcc.sign_helper"].get("expr_fm_mix", float('nan'))
    mi_r = mi_s / mi_g if mi_g > 0 else float('nan')
    ad_g = r["gcc"].get("expr_fm_add", float('nan'))
    ad_s = r["gcc.sign_helper"].get("expr_fm_add", float('nan'))
    ad_r = ad_s / ad_g if ad_g > 0 else float('nan')
    print(f"| {cfg} | {pa_g:.1f} | {pa_s:.1f} | {pa_r:.2f}x | {se_g:.1f} | {se_s:.1f} | {se_r:.2f}x | {mi_g:.1f} | {mi_s:.1f} | {mi_r:.2f}x | {ad_g:.2f} | {ad_s:.2f} | {ad_r:.2f}x |")

print("\n### Dense Matrix Multivector (`matrix_multi<double>`) Expressions")
print("| Configuration | gcc pade (ms) | gcc.sh pade (ms) | Ratio pade | gcc series (ms) | gcc.sh series (ms) | Ratio series | gcc mix (ms) | gcc.sh mix (ms) | Ratio mix | gcc add (ms) | gcc.sh add (ms) | Ratio add |")
print("| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |")
for cfg in cfgs:
    r = results[cfg]
    pa_g = r["gcc"].get("expr_mm_pade", float('nan'))
    pa_s = r["gcc.sign_helper"].get("expr_mm_pade", float('nan'))
    pa_r = pa_s / pa_g if pa_g > 0 else float('nan')
    se_g = r["gcc"].get("expr_mm_series", float('nan'))
    se_s = r["gcc.sign_helper"].get("expr_mm_series", float('nan'))
    se_r = se_s / se_g if se_g > 0 else float('nan')
    mi_g = r["gcc"].get("expr_mm_mix", float('nan'))
    mi_s = r["gcc.sign_helper"].get("expr_mm_mix", float('nan'))
    mi_r = mi_s / mi_g if mi_g > 0 else float('nan')
    ad_g = r["gcc"].get("expr_mm_add", float('nan'))
    ad_s = r["gcc.sign_helper"].get("expr_mm_add", float('nan'))
    ad_r = ad_s / ad_g if ad_g > 0 else float('nan')
    print(f"| {cfg} | {pa_g:.1f} | {pa_s:.1f} | {pa_r:.2f}x | {se_g:.1f} | {se_s:.1f} | {se_r:.2f}x | {mi_g:.1f} | {mi_s:.1f} | {mi_r:.2f}x | {ad_g:.2f} | {ad_s:.2f} | {ad_r:.2f}x |")
