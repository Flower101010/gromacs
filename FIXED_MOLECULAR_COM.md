# Experimental fixed molecular-COM constrained MD

This branch provides a narrowly scoped CPU fast path for constrained MD in
which every whole molecule has a fixed, mass-weighted center of mass (COM).
It is intended for the TIP4P/2005 water and TraPPE-UA methanol mean-force
labelling workflow in this repository. It is experimental; it is not an
upstream GROMACS input feature.

## What it constrains

For molecule `I`, the target is the mass-weighted molecular COM in the input
configuration:

```
R_I = sum_a(m_a r_Ia) / sum_a(m_a)
R_I(t) = R_I(0)
```

After GROMACS has applied SETTLE/LINCS, the fast path calculates the COM
displacement and translates every site in that molecule by `-dR_I`. It applies
the consistent velocity correction `v_Ia -= dR_I / dt` to massive particles.

Massless virtual sites do not contribute to the COM or COM velocity, but are
translated with their parent molecule. Thus the TIP4P M site is handled
correctly. A uniform translation preserves intramolecular distances and
orientation; the molecule remains free to rotate.

## Leap-frog initialization and time labels

`integrator=md` stores positions at `t_n` and velocities at `t_n-dt/2`.
The update produces `v_(n+1/2)` and then `x_(n+1) = x_n + dt*v_(n+1/2)`.
The COM correction therefore remains `-dR/dt`, not `-dR/(dt/2)`.
At trajectory time `t_n`, the force stream corresponds to `x_n`; the stored
velocity belongs to the preceding half step. The reported leap-frog kinetic
energy averages adjacent half-step kinetic energies, so it need not equal
the kinetic energy calculated from that frame's velocity alone.

For a fresh run, targets are captured before initial internal-constraint
corrections, and the initial mass-weighted molecular translational velocity
is removed before kinetic-energy history and thermostat initialization.
Both `continuation=yes` and `continuation=no` MDP settings are supported for
fresh runs; this does not enable checkpoint restart. With `continuation=no`,
the position-only constraint call and reverse initialization at `t0-dt`
are handled explicitly. A final initial velocity projection removes the
roundoff residue from that reverse step without changing the saved targets.
With `continuation=yes`, GROMACS still trusts the supplied internal geometry
and velocities; the new molecular COM constraints are initialized regardless.

The local regression can be repeated on the RTX 3060 with:

```bash
python3 admin/validate-fixed-com-initialization.py \
    --gmx /home/flos/gromacs-fixed-com-build/bin/gmx_mpi \
    --gro /home/flos/MeOH-H2O/systems/cluster_grid/xm_050/production.gro \
    --top /home/flos/CGWorkflow/benchmarks/pull-constraint-20260901/inputs/topol.local.top \
    --output /tmp/fixed-com-startup-test
```

Choose a new output directory. The script runs both continuation settings
with the feature on and off, reads high-precision TRR dumps, and checks COM
positions, half-step COM velocities, the discrete position/velocity relation,
rigid distances, initial kinetic energy and runtime Ndf. Fixture masses and
atom order are specific to the 500 SOL + 500 MEOH system above.

On 2026-09-05, the 20-step tests passed with initial maximum COM velocities
below `4e-8 nm/ps`, maximum COM position errors below `7.5e-7 nm`, and maximum
rigid-distance errors below `4.6e-7 nm`. Initial kinetic energy for
`continuation=yes` was `3697.078920 kJ/mol`, matching the independent value
`7531.697010 - 3834.618088 = 3697.078921 kJ/mol` after removing translation.
During integration, maximum COM velocity residue was `3.87e-4 nm/ps` and the
maximum discrete velocity/position discrepancy was `3.93e-4 nm/ps`, consistent
with single-precision coordinate roundoff divided by `dt=0.002 ps`. These
are finite-precision tolerances, not exact zeros or long-run ensemble tests.

## Running

Use a fresh run from the target configuration and enable the experimental
path at runtime:

```bash
export GMX_FIXED_MOLECULAR_COM=1
gmx_mpi mdrun -s input.tpr -deffnm fixed-com \
    -ntomp 8 -pin on -update cpu -nb gpu -pme gpu
```

The run log reports both the fixed-COM path and the corrected temperature
coupling degrees of freedom. For the 500 water + 500 methanol benchmark,
`nrdf` is reduced from 6000 to 3000: three translational degrees of freedom
are removed for each of 1000 molecular COM constraints. This correction is
required for the reported temperature and for v-rescale sampling.

## Online force and virial averages (local build)

Enable a final small text file without storing instantaneous trajectories:

```bash
export GMX_FIXED_MOLECULAR_COM=1
export GMX_FIXED_MOLECULAR_COM_AVERAGE_FILE="$PWD/rmf-average.txt"
export GMX_FIXED_MOLECULAR_COM_AVERAGE_SKIP_STEPS=5000
export GMX_FIXED_MOLECULAR_COM_AVERAGE_STRIDE=10
unset GMX_FIXED_MOLECULAR_COM_FORCE_FILE
gmx_mpi mdrun -s input.tpr -deffnm fixed-com \
    -ntomp 8 -update cpu -nb gpu -pme gpu
```

The average switch defaults to off. Skip defaults to 0 and stride to 1.
Both are integer MD steps: sample when `s >= skip` and
`(s - skip) % stride == 0`, where `s = step - init_step`.
Thus skip=5000 at dt=0.002 ps excludes the first 10 ps and includes the
configuration at 10 ps. The initial and final configurations are included
if they meet this rule; nsteps=N can give N+1 samples with default settings.
This schedule is independent of `nstfout`, `nstenergy`, and `nstcalcenergy`.
Set `nstfout=0`, `nstxout=0`, `nstvout=0`, and `nstxout-compressed=0` in the
TPR to avoid large trajectories. Normal small GROMACS energy/log output
can remain enabled.

The force estimator is exactly the existing molecular force stream estimator:
the same float per-molecule summation over all sites, followed by double
accumulation across samples. To keep the hot path inexpensive while reducing
long-run summation error, selected samples first accumulate into a fixed
256-sample double block and the completed block is then added to the global
double sum. This is not Kahan/Neumaier compensation and does not change the
force calculation or GPU reduction. Force rows use 1-based topology molecule
order and units kJ mol^-1 nm^-1. No constraint-force estimator is changed.

Each mapped force component and virial component is checked for a finite value,
and both block and global additions are checked for a finite result. A
non-finite value or an exhausted sample counter aborts the run with the
relative step, time, and component rather than silently writing a NaN average.

Virial rows are the nine components of the physical `force_vir` from
`do_force()` at the same configuration as the forces, before integration,
with GROMACS's analytic long-range dispersion correction removed. Thus this
is the requested microscopic configurational **no-DispCorr** virial: it
includes short-range and PME contributions plus virtual-site force handling,
but excludes only the `DispCorr` tail term. It uses GROMACS Xi convention
(the usual -1/2 force-position convention) and units kJ mol^-1.

It is **not** the later `total_vir = force_vir + shake_vir`: it contains no
kinetic stress, SETTLE/LINCS reaction virial, or fixed-COM projection reaction.
It is not reconstructed from pressure, and it is not a rerun quantity. The
fixed-COM projection does not modify either sampled force buffer or this
virial. Do not compare it directly with constraint-inclusive EDR `Vir-XX`
etc.

The implementation requests virial computation on every selected sample;
otherwise the MD loop may not have a valid virial on that step. This can
increase GPU/CPU force-evaluation cost. Stride controls that cost as well as
sample density. Accumulation uses O(number of molecules) memory; output is
written once on normal loop exit, including a graceful stop. The block
accumulator uses one additional double vector of molecule-force size and nine
additional doubles for the virial; it does not write per-sample data. Abrupt
process termination loses the averages. Checkpoint restart remains
unsupported, and multiple time stepping is rejected for this output.

The text file contains `samples`, `skip_steps`, `stride`, `beads`, one
`force bead Fx Fy Fz` row per bead and three `virial` rows (x,y,z).
With zero samples it reports `samples 0` and omits undefined force/virial
averages. Output paths are opened at startup, so use new filenames for new runs.

For validation only, `GMX_FIXED_MOLECULAR_COM_AVERAGE_VIRIAL_FILE` enables
instantaneous text rows `relative_step time_ps XX XY XZ YX YY YZ ZX ZY ZZ`
of this no-DispCorr tensor at selected samples. Leave it unset for production.
All output paths must be distinct. Repeat the local online/offline regression
with:

```bash
python3 admin/validate-fixed-com-average.py \
    --gmx /home/flos/gromacs-fixed-com-build/bin/gmx_mpi \
    --gro /home/flos/MeOH-H2O/systems/cluster_grid/xm_050/production.gro \
    --top /home/flos/CGWorkflow/benchmarks/pull-constraint-20260901/inputs/topol.local.top \
    --output /tmp/fixed-com-average-test
```

The dispersion-correction-specific regression uses normal MD only (never
`mdrun -rerun`):

```bash
python3 admin/validate-fixed-com-no-dispcorr-virial.py \
    --gmx /home/flos/gromacs-fixed-com-build/bin/gmx_mpi \
    --gro /home/flos/MeOH-H2O/systems/cluster_grid/xm_050/production.gro \
    --top /home/flos/CGWorkflow/benchmarks/pull-constraint-20260901/inputs/topol.local.top \
    --output /home/flos/.cache/fixed-com-validation/no-dispcorr-virial \
    --steps 0
```

At zero steps all runs evaluate identical initial coordinates. It checks that
`DispCorr=EnerPres` after online subtraction agrees with `DispCorr=no` within
the independent GPU repeat floor. A short multi-step run also checks that the
accumulated tensor equals the arithmetic mean of every diagnostic sample
exactly at text precision; cross-run multi-step differences include normal GPU
trajectory roundoff divergence.

Local validation on 2026-09-08 used 500 water + 500 methanol molecules on
RTX 3060, GPU nonbonded/PME and CPU update. Over 8 integration steps,
same-run online/offline force and all nine virial averages differed by
exactly 0 at the recorded precision: 9 samples with skip=0/stride=1,
6 with skip=3/stride=1, and 3 with skip=2/stride=3. Nonzero init-step=10
verified run-relative skipping. The sparse case used nstcalcenergy=100,
so its samples explicitly exercised additional virial requests.
Skip=9 produced zero samples. A run with nstfout=0 produced averages and
no TRR. Invalid skip, stride and missing fixed-COM activation failed loudly.
Separate GPU trajectories are not bitwise reproducible: average on/off
maximum instantaneous force difference was 0.02759 kJ mol^-1 nm^-1, versus
0.02728 between two average-off repeats. These short tests validate the
accumulator and integration hook, not statistical convergence of RMF labels.

## Molecular physical-force stream

Set `GMX_FIXED_MOLECULAR_COM_FORCE_FILE` to request a binary stream of
instantaneous mapped molecular forces:

```bash
export GMX_FIXED_MOLECULAR_COM=1
export GMX_FIXED_MOLECULAR_COM_FORCE_FILE="$PWD/molecular-forces.bin"
gmx_mpi mdrun -s input-with-nstfout.tpr -deffnm fixed-com \
    -ntomp 8 -pin on -update cpu -nb gpu -pme gpu
```

The TPR must have `nstfout > 0`; one record is written at each force-output
step. A record contains

```
F_map,I = sum_{a in molecule I} F_a
```

over every site, including virtual sites. It is sampled after virtual-site
force spreading and before the next coordinate update. The analytic COM
projection changes positions and velocities only; it does not add a reaction
force to this buffer. This is the physical mapped force suitable for later
conditional averaging, not `pullf.xvg` and not a constraint reaction force.

### Saved force versus `mdrun -rerun`

For this supported fixed-COM workflow, a rerun is not required to make the
mapped force physical. The saved stream is calculated by the same `do_force()`
call as `mdrun -rerun`, at the same saved coordinates, after virtual-site force
spreading and before the coordinate update. Rerun reads positions only and
does not apply the COM projection or integrate a timestep.

The following local test generates a complete per-frame TRR (coordinates,
box and forces), reruns it twice with the same TPR, and compares the saved
whole-molecule force sums against independently summed rerun atom forces:

```bash
python3 admin/validate-fixed-com-saved-vs-rerun.py \
    --gmx /home/flos/gromacs-fixed-com-build/bin/gmx_mpi \
    --gro /home/flos/MeOH-H2O/systems/cluster_grid/xm_050/production.gro \
    --top /home/flos/CGWorkflow/benchmarks/pull-constraint-20260901/inputs/topol.local.top \
    --output /tmp/fixed-com-saved-vs-rerun
```

The test records three independent links. It compares the binary stream with
atom forces in the production TRR written by `do_md_trajectory_writing`
immediately before the stream writer. It also runs the same fresh TPR twice
with `-nsteps 0`, fixed-COM OFF and ON. The ON run executes fixed-COM
initialization and internal constraints but no integration step, so it is a
direct test for force-array contamination. Finally it compares time averages
to a contiguous-block SEM from the live saved samples.

On 2026-09-08, the 500 SOL + 500 MEOH fixture was run on RTX 3060 with GPU
nonbonded/PME and CPU update. In the 101-frame, 10 ps test (sampled every
0.1 ps), stream versus same-live-buffer TRR force sums differed by at most
`1.39e-4 kJ mol^-1 nm^-1` (RMS `1.48e-5`), limited by float trajectory text
decoding. Production and rerun coordinates agreed within minimum-image
`1.21e-7 nm`; the rerun input is therefore the `do_force()` state.

At a real live `-nsteps 0` fixed-COM ON/OFF single point, mapped forces
differed by at most `9.46e-4` (RMS `1.26e-4`) kJ mol^-1 nm^-1. The ON/OFF
rerun comparison across 101 frames was `1.45e-3` maximum and `1.20e-4` RMS.
These are GPU roundoff/order differences and rule out a COM reaction force in
the sampled array.

Saved-versus-rerun instantaneous forces in the 101-frame mixed test had
maximum difference `0.02877` and RMS `0.001336 kJ mol^-1 nm^-1`, against a
largest mapped component of `1806.31` (relative maximum `1.59e-5`). Their
mean-force difference was `0.000866` maximum component and `0.000174` RMS.
The signed global `(x,y,z)` difference was
`(-1.19e-7, -1.04e-6, 4.85e-7) kJ mol^-1 nm^-1`; it has no detectable bias.
Using ten contiguous blocks of ten saved samples, mean-difference/block-SEM
had median `9.12e-6`, 95th percentile `3.89e-5`, and maximum `1.92e-4`.
This is far below the RMF sampling error; this short test estimates numerical
comparison uncertainty, not production RMF convergence.

The virtual-site and no-virtual-site paths were separately covered with 51
frames each: pure TIP4P/2005 water (`x_MeOH=0`) and pure TraPPE-UA methanol
(`x_MeOH=1`). Stream-versus-live-buffer maxima were `1.27e-4` and `1.29e-4`.
Their saved-versus-rerun mean/SEM maximum ratios were `3.56e-4` and
`1.01e-4`, respectively. All force values use kJ mol^-1 nm^-1.

Consequently, use the online average or saved force stream directly for the
same force field and run settings. Use rerun only when you intentionally need
forces under a different TPR/force field, or when an older trajectory omitted
the forces and must be evaluated afterwards. A rerun does not improve the
scientific definition or remove this normal GPU rounding difference.

The little-endian file protocol is:

| Offset | Type | Content |
| --- | --- | --- |
| 0 | 8 bytes | ASCII `GMXCOMF1` |
| 8 | `uint32` | version, currently `1` |
| 12 | `uint64` | molecule count |
| 20 | repeated | `float64 time_ps`, then `3*moleculeCount` `float32` values in topology molecule order (`Fx,Fy,Fz`) |

`admin/fixed-molecular-com-force-validation.mdp` is a five-step local input
for validating this stream. In a local RTX 3060 test, all 6000 molecular force
vectors agreed with independent sums of the corresponding TRR atom forces to
a maximum `9.69e-3 kJ mol^-1 nm^-1` (RMS `1.44e-3`), consistent with
single-precision output and text-dump rounding.

## Supported first-prototype scope

- `md` with temperature coupling, fixed volume/NVT, no FEP;
- orthorhombic `xyz` PBC and small molecules that fit within the minimum-image
  convention around one massive reference atom;
- one PP rank with all molecule atoms local;
- CPU update (`-update cpu`), with GPU nonbonded and GPU PME supported;
- whole molecules whose massive atoms are all in one temperature-coupling
  group.

The program exits rather than silently falling back when these requirements
are not met. Generic pull constraints cannot be enabled simultaneously.

## Not supported

- checkpoint continuation (initial COM targets are not checkpointed);
- multi-rank/domain-decomposition production runs;
- triclinic boxes, pressure coupling, FEP, or GPU update;
- a molecule split across temperature-coupling groups.

## Performance reference

On the V100 benchmark (one rank, eight OpenMP threads, 50,000 steps, GPU
nonbonded + GPU PME), the pre-Ndf-correction fast path measured 0.206 ms/step
and 839.8 ns/day, compared with 0.189 ms/step and 915.5 ns/day for the
feature-off build. The Ndf correction changes thermostat bookkeeping only;
nevertheless rerun the target production settings if an exact performance
number is needed.

The scientifically correct generic-pull reference used 1000 complete
molecular groups and 3000 Cartesian constraints. It measured 0.957 ms/step,
showing why this specialized path is needed. Earlier O/C single-atom pull
benchmarks are not a valid comparison for molecular-COM mapping.
