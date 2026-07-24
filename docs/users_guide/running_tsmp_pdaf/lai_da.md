(laida)=
# LAI Data Assimilation #

Leaf Area Index (LAI) data assimilation allows observed LAI fields to
be used to update the vegetation state and, optionally, the specific
leaf area parameters of eCLM.

**eCLM only.** Requires BGC (Biogeochemical Mode) to be active.
Not available for CLM3.5 or CLM5.0 without BGC.

## Configuration ##

LAI DA is controlled by three parameters in `enkfpf.par`:

- [`CLM:update_lai`](enkfpf:clm:update_lai) — selects the state vector
  layout and the location of the forward observation operator.
- [`CLM:update_lai_params`](enkfpf:clm:update_lai_params) — enables
  joint estimation of specific leaf area parameters alongside LAI.
- [`CLM:update_lai_incr_w`](enkfpf:clm:update_lai_incr_w) — blend between
  additive and multiplicative patch LAI increments (options 1 and 2 only).

(laida:background)=
## Background ##

eCLM with BGC does not carry LAI as an independent prognostic
variable.  Instead, LAI is a diagnostic quantity derived from leaf
carbon (`leafc`) and plant-functional-type (PFT) constants for
specific leaf area:

$$
\text{LAI}(p) =
\begin{cases}
  \dfrac{s_{\text{top}}(p)\,\bigl(\exp\!\bigl(C_\ell(p)\,s'(p)\bigr) - 1\bigr)}{s'(p)}
    & s'(p) > 0 \\[6pt]
  s_{\text{top}}(p)\, C_\ell(p) & s'(p) = 0
\end{cases}
$$

where
- $C_\ell$ is `leafc` (leaf carbon, gC m$^{-2}$)
- $s_{\text{top}}$ is `slatop` (specific leaf area at the canopy top,
m$^2$ gC$^{-1}$)
- $s'$ is `dsladlai` (slope of specific leaf
area with respect to LAI, m$^2$ gC$^{-1}$)

Both $s_{\text{top}}$ and $s'$ are PFT constants read from the CLM
surface dataset. The formula follows Eq. 3 of Thornton and Zimmermann
(2007, J. Clim., 20, 3902–3923).

Because LAI is not a direct model state, assimilation updates `leafc`
(and, consistently, `leafn = leafc / leafcn`) rather than LAI
itself.

## State Vector ##

The choice of [`CLM:update_lai`](enkfpf:clm:update_lai) determines the
state vector layout:

| `update_lai` | `update_lai_params` | State vector contents                                                           | Size        |
|:------------:|:-------------------:|---------------------------------------------------------------------------------|-------------|
|      1       |          0          | Gridcell-averaged LAI                                                           | $N_g$       |
|      1       |          1          | Gridcell-averaged LAI + `slatop` per patch                                      | $N_g + N_p$ |
|      2       |          0          | `leafc`, `slatop`, `dsladlai` per patch                                         | $3\,N_p$    |
|      2       |          2          | `leafc`, `slatop`, `dsladlai` per patch (parameters updated after assimilation) | $3\,N_p$    |
|      3       |          0          | `leafc`, `livestemc`, `deadstemc` per vegetated patch                           | $3\,N_v$    |
|      3       |          1          | `leafc`, `livestemc`, `deadstemc` + `slatop`, `medlynslope` per vegetated patch | $5\,N_v$    |

$N_g$ = number of local grid cells, $N_p$ = number of local patches,
$N_v$ = number of local vegetated patches (non-bare-ground PFTs).

## Option 1 — Gridcell LAI State Vector ##

With `CLM:update_lai=1`, the state vector holds one LAI value per grid
cell (a weighted average over all patches in that cell):

**Set phase** (`set_clm_statevec`):

1. Per-patch LAI is computed from `leafc` using the formula above and
   clipped to zero.
2. The gridcell LAI is the patch-weight-averaged sum:
   $\text{LAI}_g = \sum_p w_p\,\text{LAI}(p)$, where $w_p$ is
   `patch%wtgcell`.
3. Per-patch forecast LAI $\text{LAI}_\text{f}(p)$ is stored for the
   update phase;

**PDAF update**: PDAF operates on the gridcell-level LAI values.

**Update phase** (`update_clm`):

1. The gridcell analysis LAI $\text{LAI}'_g$ is mapped to patches using
   a blend of additive and multiplicative increments controlled by
   [`CLM:update_lai_incr_w`](enkfpf:clm:update_lai_incr_w)
   (default `1.0` = purely additive).
2. The updated per-patch LAI is inverted to obtain the new `leafc`:

$$
     C'_\ell(p) =
     \begin{cases}
       \dfrac{\ln\!\bigl((\text{LAI}'(p)\,s'(p)\,/\,s_{\text{top}}(p)) + 1\bigr)}{s'(p)}
         & s'(p) > 0 \\[6pt]
       \text{LAI}'(p)\,/\,s_{\text{top}}(p) & s'(p) = 0
     \end{cases}
$$

3. `leafc` is clipped to zero (no negative leaf carbon).
4. `leafn` is set consistently: $N'_\ell = C'_\ell\,/\,\text{leafcn}$.

**Joint parameter estimation** (`update_lai_params=1`): `slatop` for
each patch is appended to the state vector and updated by PDAF. The
updated `slatop` values are written back to the PFT constants before
the LAI inversion in step 2.

## Option 2 — Patch-Level State Vector with Observation Operator ##

:::{warning}
This option is not tested. Use [Option 3](laida:option3) for patch-level state vectors instead.
:::

With `CLM:update_lai=2`, the state vector holds `leafc`, `slatop`,
and `dsladlai` for every patch in three consecutive blocks:

$$
\mathbf{x} = \bigl[\underbrace{C_\ell(1),\ldots,C_\ell(N_p)}_{\text{block 1}},
              \underbrace{s_{\text{top}}(1),\ldots,s_{\text{top}}(N_p)}_{\text{block 2}},
              \underbrace{s'(1),\ldots,s'(N_p)}_{\text{block 3}}\bigr]
$$

The observation operator (`obs_op_pdaf`, active when
`CLM:update_lai=2`) maps this state to observed gridcell LAI:

$$
\mathcal{H}(\mathbf{x})_i = \sum_{p:\,g(p)=g_i} w_p\,\text{LAI}(p)
$$

where $g_i$ is the grid cell corresponding to observation $i$.

After the PDAF update:

- Patch LAI is updated with the same additive/multiplicative blend as in
  option 1 (using weighted gridcell-mean LAI for the ratio and increment).
- `leafc` is obtained by inverting the blended patch LAI and clipped to
  zero; `leafn` is updated consistently.
- If `CLM:update_lai_params=2`, `slatop` and `dsladlai` are also
  read back from blocks 2 and 3 and written to the PFT constants.
  The diagnostic LAI update is then left to eCLM's own phenology
  routines.

(laida:option3)=
## Option 3 — Carbon-Pool State Vector ##

With `CLM:update_lai=3`, the state vector holds three carbon pools for
every vegetated patch (bare-ground PFT excluded):

$$
\mathbf{x} = \bigl[\underbrace{C_\ell(1),\ldots,C_\ell(N_v)}_{\text{block 1: leafc}},
              \underbrace{C_s(1),\ldots,C_s(N_v)}_{\text{block 2: livestemc}},
              \underbrace{C_d(1),\ldots,C_d(N_v)}_{\text{block 3: deadstemc}}\bigr]
$$

where $C_\ell$ = `leafc`, $C_s$ = `livestemc`, $C_d$ = `deadstemc`.

The observation operator (`obs_op_pdaf`) maps this state to observed
gridcell LAI using fixed `dsladlai` from the PFT constants and `slatop`
either from the state vector (when `update_lai_params=1`) or from the PFT
constants:

$$
\mathcal{H}(\mathbf{x})_i = \sum_{p:\,g(p)=g_i} w_p\,\text{LAI}(p)
$$

with $\text{LAI}(p)$ computed from $C_\ell(p)$ using the Thornton and
Zimmermann (2007) formula described in the [Background](laida:background)
section.

**Set phase** (`set_clm_statevec`):

1. `leafc`, `livestemc`, and `deadstemc` for each vegetated patch are
   copied into the three consecutive state vector blocks.
2. If `update_lai_params=1`, two additional blocks are appended:
   `slatop` (block 4) and `medlynslope` (block 5) per patch, read from
   the PFT constants for each patch's PFT type.

**PDAF update**: PDAF operates on the full state vector.

**Update phase** (`update_clm`):

1. `leafc`, `livestemc`, and `deadstemc` are read from their respective
   state vector blocks and clipped to zero (no negative carbon pools).
2. Nitrogen pools are updated consistently using PFT C:N ratios from
   `pftcon`:
   - $N'_\ell = C'_\ell / \text{leafcn}$
   - $N'_s = C'_s / \text{livewdcn}$
   - $N'_d = C'_d / \text{deadwdcn}$
3. If `update_lai_params=1`, `slatop` and `medlynslope` are read from
   blocks 4 and 5 and written back to the PFT constants
   (`pftcon%slatop`, `pftcon%medlynslope`).

**Note**: Unlike options 1 and 2, the `CLM:update_lai_incr_w` blending
parameter has no effect in option 3. LAI is not part of the state vector;
only the underlying carbon pools are updated directly.

## Configuration Examples ##

**Example 1** — gridcell LAI assimilation, no parameter estimation:

```ini
[CLM]
update_lai        = 1
update_lai_params = 0
update_lai_incr_w = 1.0
```

**Example 2** — gridcell LAI assimilation with joint `slatop` estimation:

```ini
[CLM]
update_lai        = 1
update_lai_params = 1
```

**Example 3** — patch-level state vector with custom observation operator:

```ini
[CLM]
update_lai        = 2
update_lai_params = 0
```

**Example 4** — patch-level state vector with joint `slatop`/`dsladlai` estimation:

```ini
[CLM]
update_lai        = 2
update_lai_params = 2
```

**Example 5** — carbon-pool state vector, no parameter estimation:

```ini
[CLM]
update_lai        = 3
update_lai_params = 0
```

**Example 6** — carbon-pool state vector with joint `slatop`/`medlynslope` estimation:

```ini
[CLM]
update_lai        = 3
update_lai_params = 1
```
