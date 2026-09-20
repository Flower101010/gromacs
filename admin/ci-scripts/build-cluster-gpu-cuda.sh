#!/usr/bin/env bash
# Build a cluster-compatible CUDA/MPI GROMACS package in a CentOS 7 userspace.
#
# GitHub-hosted runners do not have the cluster's glibc, CUDA driver, or MPI
# stack. This script therefore uses an NVIDIA CUDA 12.4 CentOS 7 image for
# the userspace/toolchain and leaves libcuda.so.1 and the cluster MPI module
# to the target machine at runtime.

set -euo pipefail

repo_root="${GITHUB_WORKSPACE:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
artifact_dir="${GROMACS_ARTIFACT_DIR:-${repo_root}/build/artifacts}"
build_jobs="${BUILD_JOBS:-4}"
cuda_container="${CUDA_CONTAINER:-nvidia/cuda:12.4.0-devel-centos7}"

if ! command -v docker >/dev/null 2>&1; then
    echo "Docker is required on the GitHub-hosted runner." >&2
    exit 1
fi

mkdir -p "${artifact_dir}"
docker pull "${cuda_container}"

docker run --rm \
    --user 0:0 \
    -e BUILD_JOBS="${build_jobs}" \
    -e GITHUB_SHA="${GITHUB_SHA:-local}" \
    -v "${repo_root}:/src:ro" \
    -v "${repo_root}/build:/out" \
    -w /src \
    "${cuda_container}" \
    bash -s -- <<'CONTAINER_SCRIPT'
set -euo pipefail

build_jobs="${BUILD_JOBS:-4}"
commit="${GITHUB_SHA:-local}"
short_commit="${commit:0:12}"
build_root="/out/cluster-gpu-cuda-build"
install_prefix="/out/cluster-gpu-cuda-package"
artifact="/out/artifacts/gromacs-2026.3-cluster-gpu-cuda-${short_commit}.tar.gz"

# CentOS 7 is EOL, so point its repositories at the immutable vault before
# installing the build prerequisites and Software Collections GCC 11.
for repo_file in /etc/yum.repos.d/CentOS-*.repo; do
    [[ -f "${repo_file}" ]] || continue
    sed -i \
        -e 's|^mirrorlist=|#mirrorlist=|' \
        -e 's|^#baseurl=http://mirror.centos.org|baseurl=http://vault.centos.org|' \
        "${repo_file}"
done
yum clean all >/dev/null
yum install -y centos-release-scl || true
yum install -y \
    bzip2 curl file git gzip make perl tar wget xz zlib-devel \
    devtoolset-11-gcc devtoolset-11-gcc-c++

# Keep the toolchain selection explicit. CUDA 12.4 accepts GCC 11, and GCC
# 11 supplies the C++17 implementation required by current GROMACS.
# shellcheck disable=SC1091
source /opt/rh/devtoolset-11/enable

deps=/tmp/gromacs-deps
rm -rf "${deps}" "${build_root}" "${install_prefix}"
mkdir -p "${deps}" /out/artifacts

cmake_version=3.31.8
curl -fsSL \
    "https://github.com/Kitware/CMake/releases/download/v${cmake_version}/cmake-${cmake_version}-linux-x86_64.sh" \
    -o "${deps}/cmake.sh"
chmod +x "${deps}/cmake.sh"
"${deps}/cmake.sh" --skip-license --prefix=/opt/cmake
cmake_bin=/opt/cmake/bin/cmake

fftw_version=3.3.10
curl -fsSL "https://www.fftw.org/fftw-${fftw_version}.tar.gz" -o "${deps}/fftw.tar.gz"
tar -xzf "${deps}/fftw.tar.gz" -C "${deps}"
(
    cd "${deps}/fftw-${fftw_version}"
    CC=gcc CFLAGS='-O3 -fPIC' ./configure \
        --prefix=/opt/fftw \
        --enable-shared \
        --disable-static \
        --enable-float \
        --enable-openmp
    make -j"${build_jobs}"
    make install
)

# Build the same OpenMPI major/minor line provided by the cluster. The MPI
# libraries are intentionally not copied into the package: the target job
# must load its site OpenMPI/4.0.4 module so that launcher and fabric
# integration remain under site control.
openmpi_version=4.0.4
curl -fsSL \
    "https://download.open-mpi.org/release/open-mpi/v4.0/openmpi-${openmpi_version}.tar.gz" \
    -o "${deps}/openmpi.tar.gz"
tar -xzf "${deps}/openmpi.tar.gz" -C "${deps}"
(
    cd "${deps}/openmpi-${openmpi_version}"
    ./configure --prefix=/opt/openmpi --disable-static --enable-shared --without-verbs
    make -j"${build_jobs}"
    make install
)

cuda_root=/usr/local/cuda
export PATH="/opt/openmpi/bin:${cuda_root}/bin:${PATH}"
export LD_LIBRARY_PATH="/opt/openmpi/lib:/opt/fftw/lib:${cuda_root}/lib64:${LD_LIBRARY_PATH:-}"
export CMAKE_PREFIX_PATH=/opt/fftw

"${cmake_bin}" -S /src -B "${build_root}" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="${install_prefix}" \
    -DCMAKE_C_COMPILER=/opt/openmpi/bin/mpicc \
    -DCMAKE_CXX_COMPILER=/opt/openmpi/bin/mpicxx \
    -DMPI_C_COMPILER=/opt/openmpi/bin/mpicc \
    -DMPI_CXX_COMPILER=/opt/openmpi/bin/mpicxx \
    -DCUDAToolkit_ROOT="${cuda_root}" \
    -DCMAKE_INSTALL_RPATH='$ORIGIN;$ORIGIN/host;$ORIGIN/cuda;$ORIGIN/../lib;$ORIGIN/../lib/host;$ORIGIN/../lib/cuda' \
    -DGMX_MPI=ON \
    -DGMX_THREAD_MPI=OFF \
    -DGMX_GPU=CUDA \
    -DGMX_CUDA_ARCHITECTURES='70-real;80-real' \
    -DGMX_SIMD=AVX2_256 \
    -DGMX_BUILD_OWN_FFTW=OFF \
    -DGMX_FFT_LIBRARY=fftw3 \
    -DGMX_BUILD_MANUAL=OFF \
    -DGMX_BUILD_UNITTESTS=OFF \
    -DGMXAPI=OFF \
    -DBUILD_TESTING=OFF

"${cmake_bin}" --build "${build_root}" --parallel "${build_jobs}"
"${cmake_bin}" --install "${build_root}"

mkdir -p "${install_prefix}/lib/cuda" "${install_prefix}/lib/host"

# Copy the CUDA user-space libraries needed by the installed GROMACS shared
# library. libcuda.so.1 is deliberately excluded: it is supplied by the GPU
# driver on each compute node and must never be replaced by the toolkit stub.
copy_cuda_closure() {
    local pending=() seen_file dep real base
    pending=("${install_prefix}/lib/libgromacs_mpi.so.11")
    declare -A seen=()
    while ((${#pending[@]})); do
        dep="${pending[0]}"
        pending=("${pending[@]:1}")
        [[ -f "${dep}" ]] || continue
        [[ -n "${seen[${dep}]+yes}" ]] && continue
        seen["${dep}"]=1
        while read -r real; do
            [[ -n "${real}" ]] || continue
            case "${real}" in
                ${cuda_root}/*)
                    base="$(basename "${real}")"
                    case "${base}" in
                        libcuda.so*|libnvidia-ml.so*) continue ;;
                    esac
                    cp -L "${real}" "${install_prefix}/lib/cuda/${base}"
                    pending+=("${real}")
                    ;;
            esac
        done < <(
            ldd "${dep}" 2>/dev/null | awk \
                '/=> \/.*\(/ {print $3} /^[[:space:]]*\/.*\(/ {print $1}' \
                | grep '^/' || true
        )
    done
}
copy_cuda_closure

# Bundle the GCC runtime, avoiding the GCC 10.2 libstdc++ mismatch seen in
# the old build when a shell selected the wrong compiler module.
for runtime_name in libstdc++.so.6 libgcc_s.so.1 libgomp.so.1; do
    runtime_lib="$(find /opt/rh/devtoolset-11/root/usr -name "${runtime_name}" -print -quit)"
    [[ -f "${runtime_lib}" ]] || { echo "Missing ${runtime_name}" >&2; exit 1; }
    cp -L "${runtime_lib}" "${install_prefix}/lib/host/"
done

# Include the FFTW runtime as well; the target only needs its site OpenMPI
# module in addition to this package.
for fftw_lib in /opt/fftw/lib/libfftw3f.so*; do
    [[ -e "${fftw_lib}" ]] || { echo "Missing FFTW runtime" >&2; exit 1; }
    cp -L "${fftw_lib}" "${install_prefix}/lib/"
done

cat > "${install_prefix}/bin/activate" <<'ACTIVATE'
# Source this file after unpacking the package.
_gmx_cluster_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
if command -v module >/dev/null 2>&1; then
    module load OpenMPI/4.0.4 >/dev/null 2>&1 || true
fi
export GMXBIN="${_gmx_cluster_root}/bin"
export GMXLDLIB="${_gmx_cluster_root}/lib"
export GMXDATA="${_gmx_cluster_root}/share/gromacs"
export PATH="${GMXBIN}:${PATH}"
export LD_LIBRARY_PATH="${GMXLDLIB}/host:${GMXLDLIB}/cuda:${GMXLDLIB}${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
unset _gmx_cluster_root
ACTIVATE
chmod +x "${install_prefix}/bin/activate"

env -i \
    PATH="/opt/openmpi/bin:/opt/cmake/bin:/usr/bin:/bin" \
    LD_LIBRARY_PATH="${install_prefix}/lib/host:${install_prefix}/lib/cuda:${install_prefix}/lib:/opt/openmpi/lib" \
    GMXDATA="${install_prefix}/share/gromacs" \
    "${install_prefix}/bin/gmx_mpi" --version | tee /out/gromacs-gpu-version.txt

if ldd "${install_prefix}/bin/gmx_mpi" 2>&1 | grep -E 'not found' | grep -v 'libcuda.so.1' >/dev/null; then
    echo "The packaged executable has unresolved non-driver dependencies." >&2
    ldd "${install_prefix}/bin/gmx_mpi" >&2
    exit 1
fi

tar -C "${install_prefix}" -czf "${artifact}" .
sha256sum "${artifact}" | tee "${artifact}.sha256"
echo "Cluster GPU artifact: ${artifact}"
CONTAINER_SCRIPT

echo "Cluster GPU artifact(s):"
ls -lh "${artifact_dir}"/gromacs-2026.3-cluster-gpu-cuda-*.tar.gz "${artifact_dir}"/*.sha256
