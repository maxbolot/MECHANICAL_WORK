#!/usr/bin/env bash
# build_intel_netcdf.sh - Automates building NetCDF-C and NetCDF-Fortran with Intel oneAPI
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

PREFIX="${1:-$HOME/local}"
BUILD_DIR="${TMPDIR:-/tmp}/build_netcdf_oneapi_$$"

echo "========================================================================"
echo " Building NetCDF-C & NetCDF-Fortran with Intel oneAPI (icx / ifx)"
echo " Install Destination (PREFIX): $PREFIX"
echo " Temporary Build Dir:          $BUILD_DIR"
echo "========================================================================"

# 1. Initialize Intel oneAPI compilers if needed
if ! command -v ifx >/dev/null 2>&1 || ! command -v icx >/dev/null 2>&1; then
    ONEAPI_SETVARS="${ONEAPI_SETVARS:-/opt/intel/oneapi/setvars.sh}"
    if [[ -f "$ONEAPI_SETVARS" ]]; then
        echo "Sourcing Intel oneAPI from $ONEAPI_SETVARS..."
        # Intel setvars.sh relies on unbound variables and internal non-zero returns.
        # Temporarily disable -e and -u while sourcing.
        set +eu
        # shellcheck disable=SC1091
        source "$ONEAPI_SETVARS" --force >/dev/null 2>&1 || true
        set -eu
    fi
fi

if ! command -v ifx >/dev/null 2>&1 || ! command -v icx >/dev/null 2>&1; then
    echo "Error: Intel compilers (icx, ifx) not found in PATH." >&2
    echo "Please ensure Intel oneAPI is installed and source /opt/intel/oneapi/setvars.sh." >&2
    exit 1
fi

echo "Using C Compiler:       $(command -v icx) ($(icx --version | head -n 1))"
echo "Using Fortran Compiler: $(command -v ifx) ($(ifx --version | head -n 1))"

# 2. Locate HDF5 paths (Debian/Ubuntu multi-arch puts hdf5.h under /usr/include/hdf5/serial)
HDF5_INC=""
HDF5_LIB=""

if [[ -f "/usr/include/hdf5/serial/hdf5.h" ]]; then
    HDF5_INC="/usr/include/hdf5/serial"
    HDF5_LIB="/usr/lib/x86_64-linux-gnu/hdf5/serial"
elif [[ -f "/usr/include/hdf5.h" ]]; then
    HDF5_INC="/usr/include"
    HDF5_LIB="/usr/lib/x86_64-linux-gnu"
elif pkg-config --exists hdf5-serial 2>/dev/null; then
    HDF5_INC=$(pkg-config --cflags-only-I hdf5-serial 2>/dev/null | sed -e 's/-I//g' -e 's/^[ \t]*//' -e 's/[ \t]*$//' || echo "")
    HDF5_LIB=$(pkg-config --libs-only-L hdf5-serial 2>/dev/null | sed -e 's/-L//g' -e 's/^[ \t]*//' -e 's/[ \t]*$//' || echo "")
fi

echo "Detected HDF5 Include:  ${HDF5_INC:-Not found}"
echo "Detected HDF5 Lib:      ${HDF5_LIB:-Not found}"

# 3. Check essential system build dependencies
EXTRA_INC_FLAGS=()
if [[ -n "$HDF5_INC" ]]; then
    EXTRA_INC_FLAGS+=("-I$HDF5_INC")
fi

MISSING_PKGS=()
for header in "curl/curl.h" "zlib.h" "hdf5.h"; do
    if ! echo "#include <$header>" | icx "${EXTRA_INC_FLAGS[@]}" -E - >/dev/null 2>&1; then
        case "$header" in
            "curl/curl.h") MISSING_PKGS+=("libcurl4-openssl-dev") ;;
            "zlib.h")      MISSING_PKGS+=("zlib1g-dev") ;;
            "hdf5.h")      MISSING_PKGS+=("libhdf5-dev") ;;
        esac
    fi
done

if [[ ${#MISSING_PKGS[@]} -gt 0 ]]; then
    echo ""
    echo "Warning: Missing system development libraries for NetCDF: ${MISSING_PKGS[*]}"
    echo "Please install them via APT before building:"
    echo "  sudo apt update && sudo apt install -y ${MISSING_PKGS[*]} m4 cmake"
    echo ""
    read -r -p "Continue anyway? [y/N] " response
    if [[ ! "$response" =~ ^[Yy]$ ]]; then
        exit 1
    fi
fi

mkdir -p "$PREFIX/include" "$PREFIX/lib" "$PREFIX/bin" "$BUILD_DIR"
cd "$BUILD_DIR"

export PATH="$PREFIX/bin:$PATH"
export LD_LIBRARY_PATH="$PREFIX/lib:${HDF5_LIB:+$HDF5_LIB:}${LD_LIBRARY_PATH:-}"
export PKG_CONFIG_PATH="$PREFIX/lib/pkgconfig:${PKG_CONFIG_PATH:-}"
export CPPFLAGS="-I$PREFIX/include ${HDF5_INC:+-I$HDF5_INC} ${CPPFLAGS:-}"
export LDFLAGS="-L$PREFIX/lib ${HDF5_LIB:+-L$HDF5_LIB} ${LDFLAGS:-}"

NETCDF_C_VER="4.9.2"
NETCDF_FORTRAN_VER="4.6.1"

# 4. Build NetCDF-C
echo -e "\n--- Step 1/2: Downloading and building NetCDF-C v${NETCDF_C_VER} ---"
if [[ ! -f "netcdf-c-${NETCDF_C_VER}.tar.gz" ]]; then
    curl -fsSL "https://downloads.unidata.ucar.edu/netcdf-c/${NETCDF_C_VER}/netcdf-c-${NETCDF_C_VER}.tar.gz" -o "netcdf-c-${NETCDF_C_VER}.tar.gz"
fi
tar -xzf "netcdf-c-${NETCDF_C_VER}.tar.gz"

cd "netcdf-c-${NETCDF_C_VER}"
export CC=icx
export CFLAGS="-O3 -fPIC ${HDF5_INC:+-I$HDF5_INC}"

if command -v cmake >/dev/null 2>&1; then
    cmake -B build \
        -DCMAKE_INSTALL_PREFIX="$PREFIX" \
        -DCMAKE_C_COMPILER=icx \
        -DCMAKE_C_FLAGS="-O3 -fPIC ${HDF5_INC:+-I$HDF5_INC}" \
        -DCMAKE_INCLUDE_PATH="${HDF5_INC:+$HDF5_INC;}$PREFIX/include" \
        -DCMAKE_LIBRARY_PATH="${HDF5_LIB:+$HDF5_LIB;}$PREFIX/lib" \
        -DHDF5_ROOT=/usr \
        -DBUILD_SHARED_LIBS=ON \
        -DENABLE_TESTS=OFF \
        -DENABLE_DAP=ON
    cmake --build build -j"$(nproc)"
    cmake --install build
else
    ./configure --prefix="$PREFIX" \
        --enable-shared \
        --disable-dap-remote-tests \
        CPPFLAGS="-I$PREFIX/include ${HDF5_INC:+-I$HDF5_INC}" \
        LDFLAGS="-L$PREFIX/lib ${HDF5_LIB:+-L$HDF5_LIB}"
    make -j"$(nproc)"
    make install
fi

echo "NetCDF-C installed successfully in $PREFIX."
cd "$BUILD_DIR"

# 5. Build NetCDF-Fortran
echo -e "\n--- Step 2/2: Downloading and building NetCDF-Fortran v${NETCDF_FORTRAN_VER} ---"
if [[ ! -f "netcdf-fortran-${NETCDF_FORTRAN_VER}.tar.gz" ]]; then
    curl -fsSL "https://downloads.unidata.ucar.edu/netcdf-fortran/${NETCDF_FORTRAN_VER}/netcdf-fortran-${NETCDF_FORTRAN_VER}.tar.gz" -o "netcdf-fortran-${NETCDF_FORTRAN_VER}.tar.gz"
fi
tar -xzf "netcdf-fortran-${NETCDF_FORTRAN_VER}.tar.gz"

cd "netcdf-fortran-${NETCDF_FORTRAN_VER}"
export FC=ifx
export FFLAGS="-O3 -fPIC"

if command -v cmake >/dev/null 2>&1; then
    cmake -B build \
        -DCMAKE_INSTALL_PREFIX="$PREFIX" \
        -DCMAKE_C_COMPILER=icx \
        -DCMAKE_Fortran_COMPILER=ifx \
        -DCMAKE_Fortran_FLAGS="-O3 -fPIC" \
        -DCMAKE_INCLUDE_PATH="$PREFIX/include;${HDF5_INC:+$HDF5_INC}" \
        -DCMAKE_LIBRARY_PATH="$PREFIX/lib;${HDF5_LIB:+$HDF5_LIB}" \
        -DNETCDF_C_LIBRARY="$PREFIX/lib/libnetcdf.so" \
        -DNETCDF_C_INCLUDE_DIR="$PREFIX/include" \
        -DBUILD_SHARED_LIBS=ON \
        -DENABLE_TESTS=OFF
    cmake --build build -j"$(nproc)"
    cmake --install build
else
    export CPPFLAGS="-I$PREFIX/include ${HDF5_INC:+-I$HDF5_INC}"
    export LDFLAGS="-L$PREFIX/lib ${HDF5_LIB:+-L$HDF5_LIB}"
    export LIBS="-lnetcdf"
    ./configure --prefix="$PREFIX" --enable-shared
    make -j"$(nproc)"
    make install
fi

# 6. Verification
echo -e "\n========================================================================"
echo " NetCDF Build Verification"
echo "========================================================================"
if [[ -x "$PREFIX/bin/nf-config" ]]; then
    echo "nf-config location: $PREFIX/bin/nf-config"
    echo "Fortran Flags:      $("$PREFIX/bin/nf-config" --fflags)"
    echo "Fortran Libs:       $("$PREFIX/bin/nf-config" --flibs)"
    echo ""
    echo "netcdf.mod check:   $(ls -la "$PREFIX/include/netcdf.mod" 2>/dev/null || echo 'Not found')"
    echo "libnetcdff check:   $(ls -la "$PREFIX/lib/libnetcdff"* 2>/dev/null | head -n 2)"
    echo ""
    echo "SUCCESS: NetCDF with Intel oneAPI is installed and ready in $PREFIX."
else
    echo "Warning: nf-config was not found in $PREFIX/bin."
fi

# Cleanup
rm -rf "$BUILD_DIR"
