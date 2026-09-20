# GitHub Actions cluster build

`.github/workflows/cluster-build.yml` builds a CUDA/MPI GROMACS package on a
GitHub-hosted runner inside `nvidia/cuda:12.4.0-devel-centos7`. The CentOS 7
userspace keeps the package compatible with the cluster's glibc 2.17, while
CUDA 12.4 satisfies GROMACS 2026.3's CUDA requirement. The cluster does not
need outbound network access and does not need a self-hosted runner.

## Cluster requirements

The cluster only needs the following runtime pieces, all of which are already
present in the inspected environment:

- NVIDIA driver new enough to run CUDA 12.x user-space libraries;
- `OpenMPI/4.0.4`, loaded through the site's module system;
- GPU jobs submitted to a compute node, not the login node.

The package includes the CUDA user-space libraries, FFTW, and GCC runtime. It
does not include `libcuda.so.1`, because that library must come from the GPU
driver, and it does not replace the cluster's MPI library or launcher.

## Result

The workflow uploads an MPI-enabled CUDA GROMACS installation. Download the
artifact from the completed Actions run, copy it to the cluster, and unpack it
wherever you want. The package includes a relocatable `activate` script:

```bash
mkdir -p "$HOME/opt/gromacs-2026.3-cluster-gpu-cuda"
tar -xzf gromacs-2026.3-cluster-gpu-cuda-*.tar.gz \
    -C "$HOME/opt/gromacs-2026.3-cluster-gpu-cuda"
source "$HOME/opt/gromacs-2026.3-cluster-gpu-cuda/bin/activate"
gmx_mpi --version

# Example GPU job steps (run on a GPU allocation):
module load OpenMPI/4.0.4
srun --mpi=pmix -n 1 gmx_mpi mdrun -nb gpu -pme gpu ...
```

The old user build confirmed the correct approach: it used
`~/opt/cuda-12.6.0`, `GMX_GPU=CUDA`, `GMX_MPI=ON`, and `AVX2_256`. The new
workflow reproduces that toolchain without requiring the cluster to download
anything. The CUDA code is explicitly compiled for `sm_70` (V100) and `sm_80`
(A100), so one package supports both GPU types. The package still assumes the
cluster's GPU driver supports CUDA 12.x; check `nvidia-smi` inside a GPU
allocation if the first job reports a driver/runtime incompatibility.
