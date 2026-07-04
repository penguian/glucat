import os
import re


def parse_prod_output(ver_cfg_dir):
    """Parse products-8.out for framed_multi and matrix_multi times."""
    prod_fm = float("nan")
    prod_mm = float("nan")
    prod_path = os.path.join(ver_cfg_dir, "products-8.out")
    if not os.path.exists(prod_path):
        return prod_fm, prod_mm

    with open(prod_path, "r", encoding="utf-8") as f:
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


def parse_sq_output(ver_cfg_dir):
    """Parse squaring-11.out for framed_multi and matrix_multi times."""
    sq_fm = float("nan")
    sq_mm = float("nan")
    sq_path = os.path.join(ver_cfg_dir, "squaring-11.out")
    if not os.path.exists(sq_path):
        return sq_fm, sq_mm

    with open(sq_path, "r", encoding="utf-8") as f:
        content = f.read()

    # framed_multi<double>
    fm_idx = content.find("framed_multi<double>")
    if fm_idx != -1:
        pat = (
            r"Cl\(\s*11,\s*11\) in "
            r"Cl\(\s*11,\s*11\) CPU =\s*[\d\.]+"
            r" ms \(setup\)\s*([\d\.]+) \(\*\)"
        )
        m = re.search(pat, content[fm_idx:])
        if m:
            sq_fm = float(m.group(1))

    # matrix_multi<double>
    mm_idx = content.find("matrix_multi<double>")
    if mm_idx != -1:
        pat = (
            r"Cl\(\s*11,\s*11\) in "
            r"Cl\(\s*11,\s*11\) CPU =\s*[\d\.]+"
            r" ms \(setup\)\s*([\d\.]+) \(\*\)"
        )
        m = re.search(pat, content[mm_idx:])
        if m:
            sq_mm = float(m.group(1))

    return sq_fm, sq_mm


def parse_trans_output(ver_cfg_dir):
    """Parse transforms-8.out for mm and fm old/new GFFT times."""
    trans_mm_old = float("nan")
    trans_mm_new = float("nan")
    trans_fm_old = float("nan")
    trans_fm_new = float("nan")
    trans_path = os.path.join(ver_cfg_dir, "transforms-8.out")
    if not os.path.exists(trans_path):
        return trans_mm_old, trans_mm_new, trans_fm_old, trans_fm_new

    with open(trans_path, "r", encoding="utf-8") as f:
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


def load_results(base_dir, cfgs):
    """Parse products, squaring, and transforms benchmark output files."""
    results = {}
    for cfg in cfgs:
        results[cfg] = {}
        for ver in ["gcc", "gcc.sign_helper"]:
            results[cfg][ver] = {}
            ver_cfg_dir = os.path.join(base_dir, ver, cfg)

            prod_fm, prod_mm = parse_prod_output(ver_cfg_dir)
            results[cfg][ver]["prod_fm_mul"] = prod_fm
            results[cfg][ver]["prod_mm_mul"] = prod_mm

            sq_fm, sq_mm = parse_sq_output(ver_cfg_dir)
            results[cfg][ver]["sq_fm_mul"] = sq_fm
            results[cfg][ver]["sq_mm_mul"] = sq_mm

            mm_old, mm_new, fm_old, fm_new = parse_trans_output(ver_cfg_dir)
            results[cfg][ver]["trans_mm_old"] = mm_old
            results[cfg][ver]["trans_mm_new"] = mm_new
            results[cfg][ver]["trans_fm_old"] = fm_old
            results[cfg][ver]["trans_fm_new"] = fm_new

    return results


def print_row(*cols):
    """Print a list of columns as a markdown table row."""
    print("| " + " | ".join(str(c) for c in cols) + " |")


def print_products_table(results, cfgs):
    """Print multiplication benchmark results."""
    print(
        "\n### 1. products-8.out - framed_multi<double> vs "
        "matrix_multi<double> Multiplication (Cl(8,8) * at Fill: 0.5)"
    )
    print_row(
        "Configuration", "Backend",
        "gcc framed_multi (ms)", "gcc.sign_helper framed_multi (ms)", "Ratio",
        "gcc matrix_multi (ms)", "gcc.sign_helper matrix_multi (ms)", "Ratio"
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
        "\n### 2. squaring-11.out - Squaring Performance "
        "(Cl(11,11) squaring * at Fill: 0.5)"
    )
    print_row(
        "Configuration", "Backend",
        "gcc framed_multi (ms)", "gcc.sign_helper framed_multi (ms)", "Ratio",
        "gcc matrix_multi (ms)", "gcc.sign_helper matrix_multi (ms)", "Ratio"
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
        "\n### 3. transforms-8.out - Transform Performance "
        "(Cl(8,8) GFFT at Fill: 0.5)"
    )
    print_row(
        "Configuration",
        "gcc mm_old (ms)", "gcc.sign_helper mm_old (ms)", "Ratio",
        "gcc mm_new (ms)", "gcc.sign_helper mm_new (ms)", "Ratio",
        "gcc fm_old (ms)", "gcc.sign_helper fm_old (ms)", "Ratio",
        "gcc fm_new (ms)", "gcc.sign_helper fm_new (ms)", "Ratio"
    )
    print_row(*[":---"] + [":---:"] * 12)

    for cfg in cfgs:
        t = results[cfg]
        m_o_g = t["gcc"].get("trans_mm_old", float("nan"))
        m_o_s = t["gcc.sign_helper"].get("trans_mm_old", float("nan"))
        m_o_r = m_o_s / m_o_g if m_o_g > 0 else float("nan")

        m_n_g = t["gcc"].get("trans_mm_new", float("nan"))
        m_n_s = t["gcc.sign_helper"].get("trans_mm_new", float("nan"))
        m_n_r = m_n_s / m_n_g if m_n_g > 0 else float("nan")

        f_o_g = t["gcc"].get("trans_fm_old", float("nan"))
        f_o_s = t["gcc.sign_helper"].get("trans_fm_old", float("nan"))
        f_o_r = f_o_s / f_o_g if f_o_g > 0 else float("nan")

        f_n_g = t["gcc"].get("trans_fm_new", float("nan"))
        f_n_s = t["gcc.sign_helper"].get("trans_fm_new", float("nan"))
        f_n_r = f_n_s / f_n_g if f_n_g > 0 else float("nan")

        print_row(
            cfg,
            f"{m_o_g:.2f}", f"{m_o_s:.2f}", f"{m_o_r:.2f}x",
            f"{m_n_g:.2f}", f"{m_n_s:.2f}", f"{m_n_r:.2f}x",
            f"{f_o_g:.2f}", f"{f_o_s:.2f}", f"{f_o_r:.2f}x",
            f"{f_n_g:.2f}", f"{f_n_s:.2f}", f"{f_n_r:.2f}x"
        )


def main():
    """Main entrypoint for benchmark comparison."""
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
        "Comparing gcc vs gcc.sign_helper benchmarks on "
        "AMD Ryzen 7 8840HS..."
    )
    results = load_results(base_dir, cfgs)
    print_products_table(results, cfgs)
    print_squaring_table(results, cfgs)
    print_transforms_table(results, cfgs)


if __name__ == "__main__":
    main()
