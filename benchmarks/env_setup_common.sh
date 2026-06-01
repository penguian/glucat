# Common environment setup for GluCat benchmark configurations
# This script is designed to be sourced by the env-*.sh wrapper scripts.

# 1. Operating System and Core/Thread Architecture Resolution
OS_TYPE=$(uname -s)
PHYSICAL_CORES=""

if [ "$OS_TYPE" = "Linux" ]; then
    # Dynamic detection for hybrid architectures (e.g. big.LITTLE or Apple Silicon under Asahi Linux)
    if ls /sys/devices/system/cpu/cpu*/cpu_capacity >/dev/null 2>&1; then
        # Find maximum core capacity reported by the kernel (Performance cores)
        MAX_CAP=$(cat /sys/devices/system/cpu/cpu*/cpu_capacity 2>/dev/null | sort -nr | head -n 1)
        if [ -n "$MAX_CAP" ]; then
            PHYSICAL_CORES=$(cat /sys/devices/system/cpu/cpu*/cpu_capacity 2>/dev/null | grep -c "^$MAX_CAP$")
        fi
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

# Execute variable mapping
if [ "$IS_OPENMP" -eq 1 ]; then
    # OpenMP application profile: Armadillo or Eigen runs multithreaded.
    # We must restrict BLAS/OpenBLAS to serial single-threaded mode to prevent context-switching storms.
    export OMP_NUM_THREADS=$RESOLVED_THREADS
    export OPENBLAS_NUM_THREADS=1
    
    if [ "$IS_FLEXIBLAS" -eq 1 ]; then
        export FLEXIBLAS=OPENBLASSERIAL
    else
        unset FLEXIBLAS
    fi

    # Pin threads to hardware cores
    export OMP_PLACES=cores
    export OMP_PROC_BIND=close

    # Serial math environment overrides for nested loops
    export MKL_NUM_THREADS=1
    export VECLIB_MAXIMUM_THREADS=1
    export NUMEXPR_NUM_THREADS=1

else
    # Sequential application profile: Armadillo or Eigen runs strictly single-threaded.
    unset FLEXIBLAS
    
    if [ "$IS_BLAS" -eq 1 ] || [ "$IS_FLEXIBLAS" -eq 1 ] || [ "$IS_OPENBLAS" -eq 1 ]; then
        # Backends are parallelized up to the resolved thread count.
        export OMP_NUM_THREADS=$RESOLVED_THREADS
        export OPENBLAS_NUM_THREADS=$RESOLVED_THREADS
        
        if [ "$IS_FLEXIBLAS" -eq 1 ]; then
            export FLEXIBLAS=OPENBLASOPENMP
        fi

        # Pin threads if utilizing multiple threads
        if [ "$RESOLVED_THREADS" -gt 1 ]; then
            export OMP_PLACES=cores
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
