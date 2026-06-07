# Common environment setup for GluCat benchmark configurations
# This script is designed to be sourced by the env-*.sh wrapper scripts.

# 1. Operating System and Core/Thread Architecture Resolution
OS_TYPE=$(uname -s)
PHYSICAL_CORES=""

if [ "$OS_TYPE" = "Linux" ]; then
    P_CORE_LIST=""
    P_CORE_PLACES=""

    # Helper function to filter out SMT hyperthreads and keep only the primary thread of each core
    filter_smt_and_add() {
        local cpu_id=$1
        local cpu_dir="/sys/devices/system/cpu/cpu$cpu_id"
        if [ -d "$cpu_dir" ]; then
            local siblings
            siblings=$(cat "$cpu_dir/topology/thread_siblings_list" 2>/dev/null)
            local first_sibling
            first_sibling=$(echo "$siblings" | cut -d',' -f1 | cut -d'-' -f1)
            if [ "$cpu_id" = "$first_sibling" ]; then
                P_CORE_LIST="${P_CORE_LIST:+$P_CORE_LIST,}${cpu_id}"
                P_CORE_PLACES="${P_CORE_PLACES:+$P_CORE_PLACES,}{${cpu_id}}"
            fi
        fi
    }

    # 1. Option A: Intel Hybrid CPU Types API (recent Linux kernels)
    if [ -d "/sys/devices/system/cpu/types" ]; then
        for type_dir in /sys/devices/system/cpu/types/*core*; do
            if [ -f "$type_dir/cpulist" ]; then
                cpulist=$(cat "$type_dir/cpulist")
                expanded_cpus=$(python3 -c "
import sys
res = []
for part in sys.argv[1].split(','):
    if '-' in part:
        s, e = map(int, part.split('-'))
        res.extend(range(s, e+1))
    else:
        res.append(int(part))
print(' '.join(map(str, res)))
" "$cpulist" 2>/dev/null)
                for cpu_id in $expanded_cpus; do
                    filter_smt_and_add "$cpu_id"
                done
            fi
        done
    fi

    # 2. Option B: cpu_capacity Fallback (Apple Silicon, ARM64 big.LITTLE)
    if [ -z "$P_CORE_LIST" ] && ls /sys/devices/system/cpu/cpu*/cpu_capacity >/dev/null 2>&1; then
        MAX_CAP=$(cat /sys/devices/system/cpu/cpu*/cpu_capacity 2>/dev/null | sort -nr | head -n 1)
        MIN_CAP=$(cat /sys/devices/system/cpu/cpu*/cpu_capacity 2>/dev/null | sort -n | head -n 1)
        if [ -n "$MAX_CAP" ] && [ -n "$MIN_CAP" ] && [ "$MIN_CAP" -lt "$MAX_CAP" ]; then
            for cpu_dir in /sys/devices/system/cpu/cpu[0-9]*; do
                cap=$(cat "$cpu_dir/cpu_capacity" 2>/dev/null)
                if [ "$cap" = "$MAX_CAP" ]; then
                    cpu_id=$(basename "$cpu_dir" | sed 's/cpu//')
                    filter_smt_and_add "$cpu_id"
                fi
            done
        fi
    fi

    # Set PHYSICAL_CORES if asymmetric P-cores were successfully isolated
    if [ -n "$P_CORE_LIST" ]; then
        PHYSICAL_CORES=$(echo "$P_CORE_LIST" | tr ',' '\n' | wc -l)
        export GLUCAT_P_CORE_PLACES="$P_CORE_PLACES"
        export GLUCAT_P_CORE_MASK="$P_CORE_LIST"
    fi

    # Fallback to standard homogeneous x86-64 socket/core extraction
    if [ -z "$PHYSICAL_CORES" ] || [ "$PHYSICAL_CORES" -eq 0 ]; then
        SOCKETS=$(lscpu 2>/dev/null | grep -i "^Socket(s):" | awk '{print $2}')
        CORES_PER_SOCKET=$(lscpu 2>/dev/null | grep -i "^Core(s) per socket:" | awk '{print $4}')
        if [ -n "$SOCKETS" ] && [ -n "$CORES_PER_SOCKET" ]; then
            PHYSICAL_CORES=$(( SOCKETS * CORES_PER_SOCKET ))
        else
            PHYSICAL_CORES=$(nproc 2>/dev/null)
        fi
    fi
else
    # Fallback for other platforms
    PHYSICAL_CORES=8
fi

# Ensure PHYSICAL_CORES is a valid positive integer
if ! [[ "$PHYSICAL_CORES" =~ ^[0-9]+$ ]] || [ "$PHYSICAL_CORES" -le 0 ]; then
    PHYSICAL_CORES=8
fi

# 2. Pre-existing OMP_NUM_THREADS Precedence Checks
# Respect pre-existing OMP_NUM_THREADS if set to a valid positive integer, else default to PHYSICAL_CORES
if ! [[ "$OMP_NUM_THREADS" =~ ^[1-9][0-9]*$ ]]; then
    RESOLVED_THREADS=$PHYSICAL_CORES
else
    RESOLVED_THREADS=$OMP_NUM_THREADS
fi

# 3. Environment Variable Dispatcher by Configuration
if [ -z "$GLUCAT_BENCHMARKS_CONFIG" ]; then
    if [ "$GLUCAT_BENCHMARKS_VERBOSE" = "1" ] || [ "$GLUCAT_BENCHMARKS_VERBOSE" = "true" ]; then
        echo "Warning: GLUCAT_BENCHMARKS_CONFIG is not defined. Defaulting to sequential profile."
    fi
    GLUCAT_BENCHMARKS_CONFIG="sequential"
fi

# Core classification of profiles
IS_OPENMP=0
IS_BLAS=0
IS_FLEXIBLAS=0
IS_OPENBLAS=0

case "$GLUCAT_BENCHMARKS_CONFIG" in
    *openmp*)
        IS_OPENMP=1
        ;;
esac

case "$GLUCAT_BENCHMARKS_CONFIG" in
    *blas*)
        IS_BLAS=1
        ;;
esac

case "$GLUCAT_BENCHMARKS_CONFIG" in
    *flexiblas*)
        IS_FLEXIBLAS=1
        ;;
esac

case "$GLUCAT_BENCHMARKS_CONFIG" in
    *openblas*)
        IS_OPENBLAS=1
        ;;
esac

# Set OPENBLAS alternatives for FlexiBLAS
if [ -x "$(which flexiblas)" ]; then
    if flexiblas list | grep -q "OPENBLAS-SERIAL"; then
        # Fedora / Red Hat
        OPENBLAS_SERIAL="OPENBLAS-SERIAL"
    elif flexiblas list | grep -q "OPENBLASSERIAL"; then
        # Generic upstream Linux
        OPENBLAS_SERIAL="OPENBLASSERIAL"
    else
        # Reference BLAS/LAPACK
        OPENBLAS_SERIAL="NETLIB"
    fi
    if flexiblas list | grep -q "OPENBLAS-OPENMP"; then
        # Fedora / Red Hat
        OPENBLAS_OPENMP="OPENBLAS-OPENMP"
    elif flexiblas list | grep -q "OPENBLASOPENMP"; then
        # Generic upstream Linux
        OPENBLAS_OPENMP="OPENBLASOPENMP"
    else
        # Reference BLAS/LAPACK
        OPENBLAS_SERIAL="NETLIB"
    fi
else
    OPENBLAS_SERIAL="NETLIB"
    OPENBLAS_OPENMP="NETLIB"
fi

# Execute variable mapping
if [ "$IS_OPENMP" -eq 1 ]; then
    # OpenMP application profile: Armadillo or Eigen runs multithreaded.
    # We must restrict BLAS/OpenBLAS to serial single-threaded mode to prevent context-switching storms.
    export OMP_NUM_THREADS=$RESOLVED_THREADS
    export OPENBLAS_NUM_THREADS=1

    if [ "$IS_FLEXIBLAS" -eq 1 ] || [ "$IS_BLAS" -eq 1 ]; then
        export FLEXIBLAS="${OPENBLAS_SERIAL}"
    else
        unset FLEXIBLAS
    fi

    # Pin threads to hardware cores (preferring performance cores on hybrid systems)
    if [ -n "$GLUCAT_P_CORE_PLACES" ]; then
        export OMP_PLACES="$GLUCAT_P_CORE_PLACES"
    else
        export OMP_PLACES=cores
    fi
    export OMP_PROC_BIND=close

    # Serial math environment overrides for nested loops
    export MKL_NUM_THREADS=1
    export VECLIB_MAXIMUM_THREADS=1
    export NUMEXPR_NUM_THREADS=1

else
    # Sequential application profile: Armadillo or Eigen runs strictly single-threaded.
    if [ "$IS_BLAS" -eq 1 ] || [ "$IS_FLEXIBLAS" -eq 1 ] || [ "$IS_OPENBLAS" -eq 1 ]; then
        # Backends are parallelized up to the resolved thread count.
        export OMP_NUM_THREADS=$RESOLVED_THREADS
        export OPENBLAS_NUM_THREADS=$RESOLVED_THREADS

        if [ "$IS_FLEXIBLAS" -eq 1 ] || [ "$IS_BLAS" -eq 1 ]; then
            export FLEXIBLAS="${OPENBLAS_OPENMP}"
        else
            unset FLEXIBLAS
        fi

        # Pin threads if utilizing multiple threads
        if [ "$RESOLVED_THREADS" -gt 1 ]; then
            if [ -n "$GLUCAT_P_CORE_PLACES" ]; then
                export OMP_PLACES="$GLUCAT_P_CORE_PLACES"
            else
                export OMP_PLACES=cores
            fi
            export OMP_PROC_BIND=close
        else
            unset OMP_PLACES
            unset OMP_PROC_BIND
        fi
    else
        # Pure homogeneous sequential run (e.g. armadillo or eigen without BLAS or OpenMP)
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        unset OMP_PLACES
        unset OMP_PROC_BIND
    fi
fi

# Print feedback to console (optional, off by default. Enable with GLUCAT_BENCHMARKS_VERBOSE=1)
if [ "$GLUCAT_BENCHMARKS_VERBOSE" = "1" ] || [ "$GLUCAT_BENCHMARKS_VERBOSE" = "true" ]; then
    echo "=========================================================="
    echo "GluCat Environment Sourced: $GLUCAT_BENCHMARKS_CONFIG"
    echo "=========================================================="
    echo "Resolved Physical CPU Cores : $PHYSICAL_CORES"
    echo "Environment Variables Set:"
    echo "  OMP_NUM_THREADS         : $OMP_NUM_THREADS"
    echo "  OPENBLAS_NUM_THREADS    : $OPENBLAS_NUM_THREADS"
    [ -n "$FLEXIBLAS" ] && echo "  FLEXIBLAS               : $FLEXIBLAS"
    [ -n "$OMP_PLACES" ] && echo "  OMP_PLACES              : $OMP_PLACES"
    [ -n "$OMP_PROC_BIND" ] && echo "  OMP_PROC_BIND           : $OMP_PROC_BIND"
    echo "=========================================================="
fi
