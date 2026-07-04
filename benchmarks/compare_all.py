import os
import re


def parse_prod_output(cfg_dir):
    """Parse products-8.out for framed_multi and matrix_multi times."""
    prod_fm = float("nan")
    prod_mm = float("nan")
    prod_path = os.path.join(cfg_dir, "products-8.out")
    if not os.path.exists(prod_path):
        return prod_fm, prod_mm

    with open(prod_path, 'r', encoding='utf-8') as f:
        content = f.read()

    # framed_multi<double>
    fm_idx = content.find("framed_multi<double>")
    if fm_idx != -1:
        pat = (
            r"Cl\(\s*8,\s*8\) in Cl\(\s*8,\s*8\) CPU =\s*[\d\.]+"
            r" ms \(setup\)\s*([\d\.]+) \(\*\)"
        )
        m = re.search(pat, content[fm_idx:])
        if m:
            prod_fm = float(m.group(1))

    # matrix_multi<double>
    mm_idx = content.find("matrix_multi<double>")
    if mm_idx != -1:
        pat = (
            r"Cl\(\s*8,\s*8\) in Cl\(\s*8,\s*8\) CPU =\s*[\d\.]+"
            r" ms \(setup\)\s*([\d\.]+) \(\*\)"
        )
        m = re.search(pat, content[mm_idx:])
        if m:
            prod_mm = float(m.group(1))

    return prod_fm, prod_mm


def parse_sq_output(cfg_dir):
    """Parse squaring-11.out for framed_multi and matrix_multi times."""
    sq_fm = float("nan")
    sq_mm = float("nan")
    sq_path = os.path.join(cfg_dir, "squaring-11.out")
    if not os.path.exists(sq_path):
        return sq_fm, sq_mm

    with open(sq_path, 'r', encoding='utf-8') as f:
        content = f.read()

    # framed_multi<double>
    fm_idx = content.find("framed_multi<double>")
    if fm_idx != -1:
        pat = (
            r"Cl\(\s*11,\s*11\) in Cl\(\s*11,\s*11\) CPU =\s*[\d\.]+"
            r" ms \(setup\)\s*([\d\.]+) \(\*\)"
        )
        m = re.search(pat, content[fm_idx:])
        if m:
            sq_fm = float(m.group(1))

    # matrix_multi<double>
    mm_idx = content.find("matrix_multi<double>")
    if mm_idx != -1:
        pat = (
            r"Cl\(\s*11,\s*11\) in Cl\(\s*11,\s*11\) CPU =\s*[\d\.]+"
            r" ms \(setup\)\s*([\d\.]+) \(\*\)"
        )
        m = re.search(pat, content[mm_idx:])
        if m:
            sq_mm = float(m.group(1))

    return sq_fm, sq_mm


def parse_trans_output(cfg_dir):
    """Parse transforms-8.out for mm and fm old/new GFFT times."""
    trans_mm_old = float("nan")
    trans_mm_new = float("nan")
    trans_fm_old = float("nan")
    trans_fm_new = float("nan")
    trans_path = os.path.join(cfg_dir, "transforms-8.out")
    if not os.path.exists(trans_path):
        return trans_mm_old, trans_mm_new, trans_fm_old, trans_fm_new

    with open(trans_path, 'r', encoding='utf-8') as f:
        content = f.read()

    pat = (
        r"Cl\(\s*8,\s*8\) in Cl\(\s*8,\s*8\) CPU = mm:\s*([\d\.]+)"
        r"\s*\(old\)\s*([\d\.]+)\s*\(new\) fm:\s*([\d\.]+)"
        r"\s*\(old\)\s*([\d\.]+)\s*\(new\)"
    )
    m = re.search(pat, content)
    if m:
        trans_mm_old = float(m.group(1))
        trans_mm_new = float(m.group(2))
        trans_fm_old = float(m.group(3))
        trans_fm_new = float(m.group(4))

    return trans_mm_old, trans_mm_new, trans_fm_old, trans_fm_new


def parse_gfft_test_output(cfg_dir):
    """Parse gfft_test-11.out for mm and fm GFFT times."""
    gfft_mm = float("nan")
    gfft_fm = float("nan")
    gfft_path = os.path.join(cfg_dir, "gfft_test-11.out")
    if not os.path.exists(gfft_path):
        return gfft_mm, gfft_fm

    with open(gfft_path, 'r', encoding='utf-8') as f:
        content = f.read()

    pat = (
        r"R_\{-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8,"
        r"9,10,11\} in R_\{-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,1,2"
        r",3,4,5,6,7,8,9,10,11\}:\s*\n\s*CPU = mm:\s*([\d\.]+)"
        r"\s*fm:\s*([\d\.]+)"
    )
    m = re.search(pat, content)
    if m:
        gfft_mm = float(m.group(1))
        gfft_fm = float(m.group(2))

    return gfft_mm, gfft_fm


def parse_versor_output(cfg_dir):
    """Parse versor-16.out for fm and mm naive, op|, and versor times."""
    results = {}
    versor_path = os.path.join(cfg_dir, "versor-16.out")
    if not os.path.exists(versor_path):
        return results

    with open(versor_path, 'r', encoding='utf-8') as f:
        content = f.read()

    # framed_multi
    fm_idx = content.find("--- Benchmarking framed_multi ---")
    if fm_idx != -1:
        pat = (
            r"Cl\(\s*16,\s*0\) \| naive :\s*([\d\.]+) ms \| "
            r"operator\| :\s*([\d\.]+) ms \| versor :\s*([\d\.]+)"
            r" ms \| versor_exp :\s*([\d\.]+) ms"
        )
        m = re.search(pat, content[fm_idx:])
        if m:
            results["versor_fm_naive"] = float(m.group(1))
            results["versor_fm_or"] = float(m.group(2))
            results["versor_fm_ver"] = float(m.group(3))
            results["versor_fm_ver_exp"] = float(m.group(4))

    # matrix_multi
    mm_idx = content.find("--- Benchmarking matrix_multi ---")
    if mm_idx != -1:
        pat = (
            r"Cl\(\s*16,\s*0\) \| naive :\s*([\d\.]+) ms \| "
            r"operator\| :\s*([\d\.]+) ms \| versor :\s*([\d\.]+)"
            r" ms \| versor_exp :\s*([\d\.]+) ms"
        )
        m = re.search(pat, content[mm_idx:])
        if m:
            results["versor_mm_naive"] = float(m.group(1))
            results["versor_mm_or"] = float(m.group(2))
            results["versor_mm_ver"] = float(m.group(3))
            results["versor_mm_ver_exp"] = float(m.group(4))

    return results


def parse_expr_output(cfg_dir):
    """Parse expressions-8.out for fm/mm algebra expression times."""
    results = {}
    expr_path = os.path.join(cfg_dir, "expressions-8.out")
    if not os.path.exists(expr_path):
        return results

    with open(expr_path, 'r', encoding='utf-8') as f:
        content = f.read()

    # framed_multi<double>
    fm_idx = content.find("framed_multi<double>")
    if fm_idx != -1:
        pat = (
            r"Cl\(\s*8,\s*8\) in Cl\(\s*8,\s*8\) CPU =\s*[\d\.]+"
            r" ms \(setup\)\s*[\d\.]+ \(comm\)\s*([\d\.]+) \(pade\)"
            r"\s*([\d\.]+) \(series\)\s*([\d\.]+) \(mix\)"
            r"\s*([\d\.]+) \(add\)"
        )
        m = re.search(pat, content[fm_idx:])
        if m:
            results["expr_fm_pade"] = float(m.group(1))
            results["expr_fm_series"] = float(m.group(2))
            results["expr_fm_mix"] = float(m.group(3))
            results["expr_fm_add"] = float(m.group(4))

    # matrix_multi<double>
    mm_idx = content.find("matrix_multi<double>")
    if mm_idx != -1:
        pat = (
            r"Cl\(\s*8,\s*8\) in Cl\(\s*8,\s*8\) CPU =\s*[\d\.]+"
            r" ms \(setup\)\s*[\d\.]+ \(comm\)\s*([\d\.]+) \(pade\)"
            r"\s*([\d\.]+) \(series\)\s*([\d\.]+) \(mix\)"
            r"\s*([\d\.]+) \(add\)"
        )
        m = re.search(pat, content[mm_idx:])
        if m:
            results["expr_mm_pade"] = float(m.group(1))
            results["expr_mm_series"] = float(m.group(2))
            results["expr_mm_mix"] = float(m.group(3))
            results["expr_mm_add"] = float(m.group(4))

    return results


def load_results(base_dir, cfgs):
    """Parse all benchmark output files into a results dictionary."""
    results = {}
    for cfg in cfgs:
        results[cfg] = {}
        for ver in ["gcc", "gcc.sign_helper"]:
            results[cfg][ver] = {}
            cfg_dir = os.path.join(base_dir, ver, cfg)
            if not os.path.exists(cfg_dir):
                continue

            # 1. products-8
            prod_fm, prod_mm = parse_prod_output(cfg_dir)
            results[cfg][ver]["prod_fm_mul"] = prod_fm
            results[cfg][ver]["prod_mm_mul"] = prod_mm

            # 2. squaring-11
            sq_fm, sq_mm = parse_sq_output(cfg_dir)
            results[cfg][ver]["sq_fm_mul"] = sq_fm
            results[cfg][ver]["sq_mm_mul"] = sq_mm

            # 3. transforms-8
            mm_old, mm_new, fm_old, fm_new = parse_trans_output(cfg_dir)
            results[cfg][ver]["trans_mm_old"] = mm_old
            results[cfg][ver]["trans_mm_new"] = mm_new
            results[cfg][ver]["trans_fm_old"] = fm_old
            results[cfg][ver]["trans_fm_new"] = fm_new

            # 4. gfft_test-11
            gfft_mm, gfft_fm = parse_gfft_test_output(cfg_dir)
            results[cfg][ver]["gfft_mm"] = gfft_mm
            results[cfg][ver]["gfft_fm"] = gfft_fm

            # 5. versor-16
            v_res = parse_versor_output(cfg_dir)
            for k, val in v_res.items():
                results[cfg][ver][k] = val

            # 6. expressions-8
            e_res = parse_expr_output(cfg_dir)
            for k, val in e_res.items():
                results[cfg][ver][k] = val

    return results


def print_row(*cols):
    """Print a list of columns as a markdown table row."""
    print("| " + " | ".join(str(c) for c in cols) + " |")


def print_multiplications_table(results, cfgs):
    """Print products benchmark results."""
    print("## 1. Multiplications (products-8.out, Cl(8,8) * under Fill: 0.5)")
    print_row(
        "Configuration", "Backend",
        "gcc fm (ms)", "gcc.sign_helper fm (ms)", "Ratio fm",
        "gcc mm (ms)", "gcc.sign_helper mm (ms)", "Ratio mm"
    )
    print_row(*[":---"] * 2 + [":---:"] * 6)

    for cfg in cfgs:
        fm_gcc = results[cfg]["gcc"].get("prod_fm_mul", float("nan"))
        fm_sh = results[cfg]["gcc.sign_helper"].get("prod_fm_mul", float("nan"))
        fm_ratio = fm_sh / fm_gcc if fm_gcc > 0 else float("nan")
        mm_gcc = results[cfg]["gcc"].get("prod_mm_mul", float("nan"))
        mm_sh = results[cfg]["gcc.sign_helper"].get("prod_mm_mul", float("nan"))
        mm_ratio = mm_sh / mm_gcc if mm_gcc > 0 else float("nan")
        print_row(
            cfg,
            "Eigen" if "eigen" in cfg else "Armadillo",
            f"{fm_gcc:.3f}", f"{fm_sh:.3f}", f"{fm_ratio:.2f}x",
            f"{mm_gcc:.3f}", f"{mm_sh:.3f}", f"{mm_ratio:.2f}x"
        )


def print_squaring_table(results, cfgs):
    """Print squaring benchmark results."""
    print(
        "\n## 2. Squaring (squaring-11.out, Cl(11,11) squaring * "
        "under Fill: 0.5)"
    )
    print_row(
        "Configuration", "Backend",
        "gcc fm (ms)", "gcc.sign_helper fm (ms)", "Ratio fm",
        "gcc mm (ms)", "gcc.sign_helper mm (ms)", "Ratio mm"
    )
    print_row(*[":---"] * 2 + [":---:"] * 6)

    for cfg in cfgs:
        fm_gcc = results[cfg]["gcc"].get("sq_fm_mul", float("nan"))
        fm_sh = results[cfg]["gcc.sign_helper"].get("sq_fm_mul", float("nan"))
        fm_ratio = fm_sh / fm_gcc if fm_gcc > 0 else float("nan")
        mm_gcc = results[cfg]["gcc"].get("sq_mm_mul", float("nan"))
        mm_sh = results[cfg]["gcc.sign_helper"].get("sq_mm_mul", float("nan"))
        mm_ratio = mm_sh / mm_gcc if mm_gcc > 0 else float("nan")
        print_row(
            cfg,
            "Eigen" if "eigen" in cfg else "Armadillo",
            f"{fm_gcc:.1f}", f"{fm_sh:.1f}", f"{fm_ratio:.2f}x",
            f"{mm_gcc:.1f}", f"{mm_sh:.1f}", f"{mm_ratio:.2f}x"
        )


def print_transforms_table(results, cfgs):
    """Print GFFT transform benchmark results."""
    print(
        "\n## 3. GFFT Transforms (transforms-8.out, Cl(8,8) GFFT "
        "under Fill: 0.5)"
    )
    print_row(
        "Configuration",
        "gcc mm_old (ms)", "gcc.sign_helper mm_old (ms)", "Ratio mm_old",
        "gcc mm_new (ms)", "gcc.sign_helper mm_new (ms)", "Ratio mm_new",
        "gcc fm_old (ms)", "gcc.sign_helper fm_old (ms)", "Ratio fm_old",
        "gcc fm_new (ms)", "gcc.sign_helper fm_new (ms)", "Ratio fm_new"
    )
    print_row(*[":---"] + [":---:"] * 12)

    for cfg in cfgs:
        r = results[cfg]
        mo_g = r["gcc"].get("trans_mm_old", float("nan"))
        mo_s = r["gcc.sign_helper"].get("trans_mm_old", float("nan"))
        mo_r = mo_s / mo_g if mo_g > 0 else float("nan")
        mn_g = r["gcc"].get("trans_mm_new", float("nan"))
        mn_s = r["gcc.sign_helper"].get("trans_mm_new", float("nan"))
        mn_r = mn_s / mn_g if mn_g > 0 else float("nan")
        fo_g = r["gcc"].get("trans_fm_old", float("nan"))
        fo_s = r["gcc.sign_helper"].get("trans_fm_old", float("nan"))
        fo_r = fo_s / fo_g if fo_g > 0 else float("nan")
        fn_g = r["gcc"].get("trans_fm_new", float("nan"))
        fn_s = r["gcc.sign_helper"].get("trans_fm_new", float("nan"))
        fn_r = fn_s / fn_g if fn_g > 0 else float("nan")
        print_row(
            cfg,
            f"{mo_g:.2f}", f"{mo_s:.2f}", f"{mo_r:.2f}x",
            f"{mn_g:.2f}", f"{mn_s:.2f}", f"{mn_r:.2f}x",
            f"{fo_g:.2f}", f"{fo_s:.2f}", f"{fo_r:.2f}x",
            f"{fn_g:.2f}", f"{fn_s:.2f}", f"{fn_r:.2f}x"
        )


def print_gfft_test_table(results, cfgs):
    """Print GFFT high-dimension benchmark results."""
    print(
        "\n## 4. GFFT Test (gfft_test-11.out, R_{-11,11} GFFT "
        "under Fill: 0.5)"
    )
    print_row(
        "Configuration",
        "gcc mm (ms)", "gcc.sign_helper mm (ms)", "Ratio mm",
        "gcc fm (ms)", "gcc.sign_helper fm (ms)", "Ratio fm"
    )
    print_row(*[":---"] + [":---:"] * 6)

    for cfg in cfgs:
        mm_gcc = results[cfg]["gcc"].get("gfft_mm", float("nan"))
        mm_sh = results[cfg]["gcc.sign_helper"].get("gfft_mm", float("nan"))
        mm_ratio = mm_sh / mm_gcc if mm_gcc > 0 else float("nan")
        fm_gcc = results[cfg]["gcc"].get("gfft_fm", float("nan"))
        fm_sh = results[cfg]["gcc.sign_helper"].get("gfft_fm", float("nan"))
        fm_ratio = fm_sh / fm_gcc if fm_gcc > 0 else float("nan")
        print_row(
            cfg,
            f"{mm_gcc:.1f}", f"{mm_sh:.1f}", f"{mm_ratio:.2f}x",
            f"{fm_gcc:.1f}", f"{fm_sh:.1f}", f"{fm_ratio:.2f}x"
        )


def print_versors_table(results, cfgs):
    """Print versor product benchmark results."""
    print(
        "\n## 5. Versor Product Benchmarks (versor-16.out, Cl(16,0) under "
        "Fill: 0.5)"
    )
    print("### Framed Multivector (`framed_multi`) Versors")
    print_row(
        "Configuration",
        "gcc naive (ms)", "gcc.sh naive (ms)", "Ratio naive",
        "gcc op| (ms)", "gcc.sh op| (ms)", "Ratio op|",
        "gcc versor (ms)", "gcc.sh versor (ms)", "Ratio versor",
        "gcc v_exp (ms)", "gcc.sh v_exp (ms)", "Ratio v_exp"
    )
    print_row(*[":---"] + [":---:"] * 12)

    for cfg in cfgs:
        r = results[cfg]
        n_g = r["gcc"].get("versor_fm_naive", float("nan"))
        n_s = r["gcc.sign_helper"].get("versor_fm_naive", float("nan"))
        n_r = n_s / n_g if n_g > 0 else float("nan")
        o_g = r["gcc"].get("versor_fm_or", float("nan"))
        o_s = r["gcc.sign_helper"].get("versor_fm_or", float("nan"))
        o_r = o_s / o_g if o_g > 0 else float("nan")
        v_g = r["gcc"].get("versor_fm_ver", float("nan"))
        v_s = r["gcc.sign_helper"].get("versor_fm_ver", float("nan"))
        v_r = v_s / v_g if v_g > 0 else float("nan")
        ve_g = r["gcc"].get("versor_fm_ver_exp", float("nan"))
        ve_s = r["gcc.sign_helper"].get("versor_fm_ver_exp", float("nan"))
        ve_r = ve_s / ve_g if ve_g > 0 else float("nan")
        print_row(
            cfg,
            f"{n_g:.1f}", f"{n_s:.1f}", f"{n_r:.2f}x",
            f"{o_g:.1f}", f"{o_s:.1f}", f"{o_r:.2f}x",
            f"{v_g:.1f}", f"{v_s:.1f}", f"{v_r:.2f}x",
            f"{ve_g:.1f}", f"{ve_s:.1f}", f"{ve_r:.2f}x"
        )

    print("\n### Dense Matrix Multivector (`matrix_multi`) Versors")
    print_row(
        "Configuration",
        "gcc naive (ms)", "gcc.sh naive (ms)", "Ratio naive",
        "gcc op| (ms)", "gcc.sh op| (ms)", "Ratio op|",
        "gcc versor (ms)", "gcc.sh versor (ms)", "Ratio versor",
        "gcc v_exp (ms)", "gcc.sh v_exp (ms)", "Ratio v_exp"
    )
    print_row(*[":---"] + [":---:"] * 12)

    for cfg in cfgs:
        r = results[cfg]
        n_g = r["gcc"].get("versor_mm_naive", float("nan"))
        n_s = r["gcc.sign_helper"].get("versor_mm_naive", float("nan"))
        n_r = n_s / n_g if n_g > 0 else float("nan")
        o_g = r["gcc"].get("versor_mm_or", float("nan"))
        o_s = r["gcc.sign_helper"].get("versor_mm_or", float("nan"))
        o_r = o_s / o_g if o_g > 0 else float("nan")
        v_g = r["gcc"].get("versor_mm_ver", float("nan"))
        v_s = r["gcc.sign_helper"].get("versor_mm_ver", float("nan"))
        v_r = v_s / v_g if v_g > 0 else float("nan")
        ve_g = r["gcc"].get("versor_mm_ver_exp", float("nan"))
        ve_s = r["gcc.sign_helper"].get("versor_mm_ver_exp", float("nan"))
        ve_r = ve_s / ve_g if ve_g > 0 else float("nan")
        print_row(
            cfg,
            f"{n_g:.1f}", f"{n_s:.1f}", f"{n_r:.2f}x",
            f"{o_g:.1f}", f"{o_s:.1f}", f"{o_r:.2f}x",
            f"{v_g:.1f}", f"{v_s:.1f}", f"{v_r:.2f}x",
            f"{ve_g:.1f}", f"{ve_s:.1f}", f"{ve_r:.2f}x"
        )


def print_expressions_table(results, cfgs):
    """Print expression algebra benchmark results."""
    print(
        "\n## 6. Algebra Expressions (expressions-8.out, Cl(8,8) "
        "under Fill: 0.5)"
    )
    print("### Framed Multivector (`framed_multi<double>`) Expressions")
    print_row(
        "Configuration",
        "gcc pade (ms)", "gcc.sh pade (ms)", "Ratio pade",
        "gcc series (ms)", "gcc.sh series (ms)", "Ratio series",
        "gcc mix (ms)", "gcc.sh mix (ms)", "Ratio mix",
        "gcc add (ms)", "gcc.sh add (ms)", "Ratio add"
    )
    print_row(*[":---"] + [":---:"] * 12)

    for cfg in cfgs:
        r = results[cfg]
        pa_g = r["gcc"].get("expr_fm_pade", float("nan"))
        pa_s = r["gcc.sign_helper"].get("expr_fm_pade", float("nan"))
        pa_r = pa_s / pa_g if pa_g > 0 else float("nan")
        se_g = r["gcc"].get("expr_fm_series", float("nan"))
        se_s = r["gcc.sign_helper"].get("expr_fm_series", float("nan"))
        se_r = se_s / se_g if se_g > 0 else float("nan")
        mi_g = r["gcc"].get("expr_fm_mix", float("nan"))
        mi_s = r["gcc.sign_helper"].get("expr_fm_mix", float("nan"))
        mi_r = mi_s / mi_g if mi_g > 0 else float("nan")
        ad_g = r["gcc"].get("expr_fm_add", float("nan"))
        ad_s = r["gcc.sign_helper"].get("expr_fm_add", float("nan"))
        ad_r = ad_s / ad_g if ad_g > 0 else float("nan")
        print_row(
            cfg,
            f"{pa_g:.1f}", f"{pa_s:.1f}", f"{pa_r:.2f}x",
            f"{se_g:.1f}", f"{se_s:.1f}", f"{se_r:.2f}x",
            f"{mi_g:.1f}", f"{mi_s:.1f}", f"{mi_r:.2f}x",
            f"{ad_g:.2f}", f"{ad_s:.2f}", f"{ad_r:.2f}x"
        )

    print("\n### Dense Matrix Multivector (`matrix_multi<double>`) Expressions")
    print_row(
        "Configuration",
        "gcc pade (ms)", "gcc.sh pade (ms)", "Ratio pade",
        "gcc series (ms)", "gcc.sh series (ms)", "Ratio series",
        "gcc mix (ms)", "gcc.sh mix (ms)", "Ratio mix",
        "gcc add (ms)", "gcc.sh add (ms)", "Ratio add"
    )
    print_row(*[":---"] + [":---:"] * 12)

    for cfg in cfgs:
        r = results[cfg]
        pa_g = r["gcc"].get("expr_mm_pade", float("nan"))
        pa_s = r["gcc.sign_helper"].get("expr_mm_pade", float("nan"))
        pa_r = pa_s / pa_g if pa_g > 0 else float("nan")
        se_g = r["gcc"].get("expr_mm_series", float("nan"))
        se_s = r["gcc.sign_helper"].get("expr_mm_series", float("nan"))
        se_r = se_s / se_g if se_g > 0 else float("nan")
        mi_g = r["gcc"].get("expr_mm_mix", float("nan"))
        mi_s = r["gcc.sign_helper"].get("expr_mm_mix", float("nan"))
        mi_r = mi_s / mi_g if mi_g > 0 else float("nan")
        ad_g = r["gcc"].get("expr_mm_add", float("nan"))
        ad_s = r["gcc.sign_helper"].get("expr_mm_add", float("nan"))
        ad_r = ad_s / ad_g if ad_g > 0 else float("nan")
        print_row(
            cfg,
            f"{pa_g:.1f}", f"{pa_s:.1f}", f"{pa_r:.2f}x",
            f"{se_g:.1f}", f"{se_s:.1f}", f"{se_r:.2f}x",
            f"{mi_g:.1f}", f"{mi_s:.1f}", f"{mi_r:.2f}x",
            f"{ad_g:.2f}", f"{ad_s:.2f}", f"{ad_r:.2f}x"
        )


def main():
    """Main entrypoint for comprehensive benchmark comparison."""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    base_dir = os.path.join(script_dir, "AMD-Ryzen-7-8840HS")
    cfgs = [
        "armadillo", "armadillo-blas", "armadillo-blas-openmp",
        "armadillo-flexiblas", "armadillo-flexiblas-openmp",
        "armadillo-openblas", "armadillo-openblas-openmp",
        "armadillo-openmp",
        "eigen", "eigen-blas", "eigen-blas-openmp",
        "eigen-flexiblas", "eigen-flexiblas-openmp",
        "eigen-openblas", "eigen-openblas-openmp",
        "eigen-openmp"
    ]

    print(
        "# Performance Comparison: GCC baseline vs GCC sign_helper on "
        "AMD Ryzen 7 8840HS"
    )
    print(
        "This report evaluates the performance impact of the `sign_helper` "
        "refactoring across 16 different Eigen/Armadillo backend "
        "configurations.\n"
    )

    results = load_results(base_dir, cfgs)
    print_multiplications_table(results, cfgs)
    print_squaring_table(results, cfgs)
    print_transforms_table(results, cfgs)
    print_gfft_test_table(results, cfgs)
    print_versors_table(results, cfgs)
    print_expressions_table(results, cfgs)


if __name__ == "__main__":
    main()
