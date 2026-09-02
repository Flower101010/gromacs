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
