# TLMC — Tiled Lattice Monte Carlo

**Tiled Lattice Monte Carlo (TLMC)** enables **large-scale canonical MC (CMC) and kinetic MC (KMC)** on lattice models by removing the *per-site topology/neighbour-list memory bottleneck*.  
It is designed to scale to **O(10⁸+) lattice sites (≈ 100+ million atoms)** on HPC by storing neighbour/topology information **once** for a small reference cell and reusing it across a tiled supercell.

<p align="center">
  <img src="docs/images/prim_supercell_to_tiled_supercell.png" alt="TLMC framework" width="800">
</p>

---

## Why TLMC?

Conventional lattice MC/KMC stores the species/vacancy on every lattice site **and** often stores additional *per-site metadata* (e.g., neighbour lists used for ΔE and barrier evaluation).  
For very large supercells, storing neighbour/topology information for *every* site becomes the dominant memory cost.

**TLMC avoids replicating topology data.** It stores:
- A single **global occupation array** for the full tiled system (unavoidable in any lattice method)
- **Neighbour/topology tables only for a small reference configuration** (the *Prim Supercell*)
- Compact binary dumps of the occupation state for restart/post-processing

This makes very large MC/KMC feasible without the neighbour-list memory blow-up.

---

## Core Idea

### 1) Prim Supercell (PrimSC) → Tiled Supercell
Start with a **Prim Supercell** containing `Nprim` lattice sites (fixed geometry + site coordinates).  
Build the large system by replicating PrimSC on a 3D tiling `(nx, ny, nz)`:

- `Ntiles = nx * ny * nz`
- `N = Ntiles * Nprim`

Every tile has **identical** lattice geometry and neighbour topology; **only occupations differ**.

### 2) Global indexing + atom vector
Each site in the tiled system is labeled by:
- `p` = PrimSC replica (tile) id
- `i` = site id inside PrimSC

Global id:
`g = p * Nprim + i`

The simulation state is stored as one contiguous atom vector:
`{ s_g } (g = 0..N-1)`  
where `s_g` is the element/vacancy at site `g`.

### 3) Neighbour lookup tables (stored once)
For each PrimSC site `i`, we precompute a neighbour list where each neighbour is encoded as:

`(i_nbr, κ)`

- `i_nbr` = neighbour site id inside PrimSC  
- `κ` = neighbour-tile pointer
  - `κ = -1` if the neighbour is inside the *same* tile
  - otherwise `κ` indexes into a fixed, sorted list of neighbouring tile ids for tile `p`

During MC/KMC:
1. Read PrimSC neighbour list for `i`
2. Recover neighbour tile id `p'` from `κ` (either `p` or lookup in the tile-neighbour list)
3. Convert `(i_nbr, p') → g_nbr` using `g = p' * Nprim + i_nbr`

✅ Exact connectivity across tile boundaries  
✅ No per-site neighbour list replication across `N`

---

## Physics models supported

### Cluster expansion Hamiltonian
TLMC uses a **cluster expansion (CE) effective Hamiltonian** for fast energy evaluation:
- Fit CE coefficients using **CASM**, **icet**, or any CE workflow (including orthogonal polynomial bases such as **Chebyshev** for appropriate DOFs), then export coefficients in the format TLMC expects.
- Works naturally for multi-component alloys with occupational degrees of freedom.

### KMC with kinetically-resolved activation (KRA) barriers
For KMC, TLMC supports barrier models based on a **kinetically-resolved activation (KRA)** framework, where a direction-independent resolved barrier is used as the basis for hop barriers (classic Van der Ven / Ceder formulation).

> If you use a different barrier model, you can implement it as long as you can compute event rates consistently from local environments.

---

## Installation

See the full installation guide here: [docs/installation_guide.md](docs/installation_guide.md)


---

## Lineage

## Lineage

TLMC extends https://github.com/pravendra12/LatticeMonteCarlo-eigen (BCC/FCC CE-LMC; HCP untested)  

which is based on https://github.com/zhucongx/LatticeMonteCarlo (original framework).

---



---

## Citations / References

### Cluster expansion tools

* CASM: [https://prisms-center.github.io/CASMcode_docs/](https://prisms-center.github.io/CASMcode_docs/)
* icet: [https://icet.materialsmodeling.org/](https://icet.materialsmodeling.org/)
* Ångqvist et al., *icet — A Python library for constructing and sampling alloy cluster expansions*, arXiv:1901.08790 (2019)

### KMC / KRA barrier framework

* Van der Ven et al., *First-principles theory of ionic diffusion with nondilute carriers*, Phys. Rev. B **64**, 184307 (2001)

### Upstream Framework Reference (Original LatticeMonteCarlo)

- Z. Xi *et al.*, “Kinetic Monte Carlo simulations of solute clustering during quenching and aging of Al–Mg–Zn alloys,” *Acta Materialia* **269** (2024) 119795. https://doi.org/10.1016/j.actamat.2024.119795


---

## Contact / Contributing

Issues and PRs are welcome. If you use TLMC in your work, please cite the references above (and add your TLMC paper/thesis citation here once public).

