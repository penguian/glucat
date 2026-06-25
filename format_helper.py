#!/usr/bin/env python3
import sys
import re
import subprocess

def preprocess(content):
    # Transforms:
    # auto blah()
    # { return foo; }
    # into:
    # auto blah()
    # {
    #   return foo;
    # }
    pattern = re.compile(r"\n(\s*)\{\s*(return\s+[^;]+;)\s*\}")
    return pattern.sub(r"\n\1{\n\1  \2\n\1}", content)

def postprocess(content):
    # Transforms:
    # auto blah()
    # {
    #   return foo;
    # }
    # back into:
    # auto blah()
    # { return foo; }
    pattern = re.compile(r"\n(\s*)\{\n\s*(return\s+[^;]+;)\n\s*\}")
    return pattern.sub(r"\n\1{ \2 }", content)


def check_file(filepath):
    if "glucat_config.h" in filepath:
        return True
    try:
        with open(filepath, "r") as f:
            content = f.read()
    except Exception as e:
        print(f"Error reading {filepath}: {e}", file=sys.stderr)
        return False

    preprocessed = preprocess(content)
    
    process = subprocess.Popen(
        ["clang-format", "-dry-run", "-Werror", f"--assume-filename={filepath}"],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )
    stdout, stderr = process.communicate(input=preprocessed)
    
    if process.returncode != 0:
        # Print compiler-like diagnostics for the file
        print(f"Formatting violations in {filepath}:", file=sys.stderr)
        sys.stderr.write(stderr)
        return False
    return True

def fix_file(filepath):
    if "glucat_config.h" in filepath:
        return True
    try:
        with open(filepath, "r") as f:
            content = f.read()
    except Exception as e:
        print(f"Error reading {filepath}: {e}", file=sys.stderr)
        return False


    preprocessed = preprocess(content)
    
    process = subprocess.Popen(
        ["clang-format", f"--assume-filename={filepath}"],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )
    formatted, stderr = process.communicate(input=preprocessed)
    
    if process.returncode != 0:
        print(f"clang-format error on {filepath}:", file=sys.stderr)
        sys.stderr.write(stderr)
        return False
        
    postprocessed = postprocess(formatted)
    
    if postprocessed != content:
        with open(filepath, "w") as f:
            f.write(postprocessed)
        print(f"Formatted {filepath}")
    return True

def main():
    if len(sys.argv) < 3:
        print("Usage: format_helper.py [check|fix] [files...]", file=sys.stderr)
        sys.exit(1)
        
    mode = sys.argv[1]
    files = sys.argv[2:]
    
    success = True
    if mode == "check":
        for filepath in files:
            if not check_file(filepath):
                success = False
    elif mode == "fix":
        for filepath in files:
            if not fix_file(filepath):
                success = False
    else:
        print(f"Unknown mode: {mode}", file=sys.stderr)
        sys.exit(1)
        
    if not success:
        sys.exit(1)

if __name__ == "__main__":
    main()
