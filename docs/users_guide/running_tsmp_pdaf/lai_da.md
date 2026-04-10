(laida)=
# LAI Data Assimilation #

Leaf Area Index (LAI) data assimilation allows observed LAI fields to
be used to update the vegetation state and, optionally, the specific
leaf area parameters of eCLM.

**eCLM only.** Requires BGC (Biogeochemical Mode) to be active.
Not available for CLM3.5 or CLM5.0 without BGC.

## Configuration ##

LAI DA is controlled by two parameters in `enkfpf.par`:

- [`CLM:update_lai`](enkfpf:clm:update_lai) — selects the state vector
  layout and the location of the forward observation operator.
- [`CLM:update_lai_params`](enkfpf:clm:update_lai_params) — enables
  joint estimation of specific leaf area parameters alongside LAI.

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

| `update_lai` | `update_lai_params` | State vector contents | Size |
|:---:|:---:|---|---|
| 1 | 0 | Gridcell-averaged LAI | $N_g$ |
| 1 | 1 | Gridcell-averaged LAI + `slatop` per patch | $N_g + N_p$ |
| 2 | 0 | `leafc`, `slatop`, `dsladlai` per patch | $3\,N_p$ |
| 2 | 2 | `leafc`, `slatop`, `dsladlai` per patch (parameters updated after assimilation) | $3\,N_p$ |

$N_g$ = number of local grid cells, $N_p$ = number of local patches.

## Option 1 — Gridcell LAI State Vector ##

With `CLM:update_lai=1`, the state vector holds one LAI value per grid
cell (a weighted average over all patches in that cell):

**Set phase** (`set_clm_statevec`):

1. Per-patch LAI is computed from `leafc` using the formula above and
   clipped to zero.
2. The gridcell LAI is the patch-weight-averaged sum:
   $\text{LAI}_g = \sum_p w_p\,\text{LAI}(p)$, where $w_p$ is
   `patch%wtgcell`.
3. The relative contribution of each patch is stored as
   $f_p = w_p\,\text{LAI}(p)\,/\,\text{LAI}_g$, which is used in the
   update phase.

**PDAF update**: PDAF operates on the gridcell-level LAI values.

**Update phase** (`update_clm`):

1. The updated gridcell LAI from the state vector is distributed back
   to patches: $\text{LAI}'(p) = \text{LAI}'_g\,f_p\,/\,w_p$.
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

- `leafc` is read back from block 1 and clipped to zero;
  `leafn` is updated consistently.
- If `CLM:update_lai_params=2`, `slatop` and `dsladlai` are also
  read back from blocks 2 and 3 and written to the PFT constants.
  The diagnostic LAI update is then left to eCLM's own phenology
  routines.

## Configuration Examples ##

**Example 1** — gridcell LAI assimilation, no parameter estimation:

```ini
[CLM]
update_lai        = 1
update_lai_params = 0
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
