import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D
import warnings
warnings.filterwarnings('ignore')

# 0. Parameter settings
N_SIM = 10000
RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)

COMPOUNDS = ['IMI', 'THM', 'DIN', 'ACE', 'CLO', 'NIP']

# 1. Load input data
print("=" * 60)
print("Loading input data...")

rq_global = pd.read_csv('RQ_Global.txt', sep='\t')
rq_global.columns = rq_global.columns.str.strip()
print(f"Global RQ data: {rq_global.shape[0]} rows × {rq_global.shape[1]} columns")
print(f"Site counts by continent: \n{rq_global['Continent'].value_counts().to_string()}")

mag = pd.read_csv('magnification.txt', sep='\t')
mag.columns = mag.columns.str.strip()
print(f"\nRisk increment factor data: {mag.shape[0]} rows × {mag.shape[1]} columns")
print(f"Column names: {list(mag.columns)}")

# 2. Extract compound-specific risk increment factors and exclude invalid values (≤1 indicates no risk increment)
print("\n" + "=" * 60)
print("Descriptive statistics for risk increment factors (valid values >1 only):")

mag_data = {}
for cpd in COMPOUNDS:
    col = f't-{cpd}_magnification'
    vals = mag[col].values
    vals_valid = vals[vals > 1.0]
    mag_data[cpd] = vals_valid
    print(f"  {cpd}: n={len(vals_valid)}, "
          f"min={vals_valid.min():.2f}, "
          f"median={np.median(vals_valid):.2f}, "
          f"max={vals_valid.max():.1f}, "
          f"mean={vals_valid.mean():.2f}")

aneos_vals = mag['a-NEOs_magnification'].values
aneos_valid = aneos_vals[aneos_vals > 1.0]
print(f"  a-NEOs: n={len(aneos_valid)}, "
      f"min={aneos_valid.min():.2f}, "
      f"median={np.median(aneos_valid):.2f}, "
      f"max={aneos_valid.max():.1f}, "
      f"mean={aneos_valid.mean():.2f}")

# 3. Fit log-normal distributions
print("\n" + "=" * 60)
print("Log-normal distribution parameters (mu=mean of log-transformed values, sigma=standard deviation of log-transformed values):")

fit_params = {}

def fit_lognormal(data, name):
    log_data = np.log(data)
    mu = log_data.mean()
    sigma = log_data.std()
    if len(log_data) <= 5000:
        stat, p = stats.shapiro(log_data)
        norm_ok = "Passed" if p > 0.05 else f"Failed: p={p:.4f}"
    else:
        norm_ok = "Skipped: sample size exceeds 5000"
    print(f"  {name}: mu={mu:.3f}, sigma={sigma:.3f}, "
          f"Normality test: {norm_ok}")
    return mu, sigma

for cpd in COMPOUNDS:
    fit_params[cpd] = fit_lognormal(mag_data[cpd], f't-{cpd}')

fit_params['a-NEOs'] = fit_lognormal(aneos_valid, 'a-NEOs')

# 4. Extract global parent-NEO RQs (excluding zero values)
print("\n" + "=" * 60)
print("Descriptive statistics for global parent-NEO RQs (positive values only):")

rq_data = {}
rq_fit = {}

for cpd in COMPOUNDS:
    col = f'RQ-{cpd}'
    vals = rq_global[col].values
    vals_nonzero = vals[vals > 0]
    rq_data[cpd] = vals_nonzero
    if len(vals_nonzero) > 3:
        mu_rq, sigma_rq = fit_lognormal(vals_nonzero, f'RQ-{cpd}')
        rq_fit[cpd] = (mu_rq, sigma_rq)
        print(f"  {cpd}: n={len(vals_nonzero)}, "
              f"median={np.median(vals_nonzero):.4f}, "
              f"max={vals_nonzero.max():.1f}")
    else:
        rq_fit[cpd] = None
        print(f"  {cpd}: Insufficient valid data (n={len(vals_nonzero)}); skipped")

# 5. Monte Carlo simulation for individual compounds
print("\n" + "=" * 60)
print(f"Starting Monte Carlo simulation (N={N_SIM} iterations)...")

results = {}

for cpd in COMPOUNDS:
    if rq_fit[cpd] is None:
        print(f"  Skipping {cpd} (insufficient data)")
        continue

    mu_mag, sigma_mag = fit_params[cpd]
    mu_rq, sigma_rq = rq_fit[cpd]
    n_sites = len(rq_data[cpd])

    rq_orig_2d = rq_data[cpd][:, np.newaxis]

    sampled_mag = np.random.lognormal(
        mean=mu_mag, sigma=sigma_mag,
        size=(n_sites, N_SIM)
    )
    sampled_mag = np.maximum(sampled_mag, 1.0)
    rq_adjusted_2d = rq_orig_2d * sampled_mag

    # Calculate the median across sites for each simulation iteration
    # to derive the 95% simulation interval of the risk increment ratio
    col_medians = np.median(rq_adjusted_2d, axis=0)   # shape: (N_SIM,)
    orig_median = np.median(rq_data[cpd])
    col_ratios  = col_medians / orig_median if orig_median > 0 else np.full(N_SIM, np.nan)

    results[cpd] = {
        'rq_original':  rq_data[cpd],
        'rq_adjusted':  rq_adjusted_2d,
        'sim_medians':  col_medians,
        'sim_means':    np.mean(rq_adjusted_2d, axis=0),
        'ratio_p2_5':     np.percentile(col_ratios, 2.5),
        'ratio_p50':    np.percentile(col_ratios, 50),
        'ratio_p97_5':    np.percentile(col_ratios, 97.5),
    }
    print(f"  {cpd}: Completed {n_sites} sites × {N_SIM} simulation iterations")

# Monte Carlo simulation for a-NEOs
print("  Simulating a-NEOs (total p-NEO RQ × analogue increment factor)...")
mu_aneos, sigma_aneos = fit_params['a-NEOs']

rq_pNEOs_list = []
for i in range(len(rq_global)):
    total = sum(rq_global[f'RQ-{c}'].iloc[i] for c in COMPOUNDS)
    rq_pNEOs_list.append(total)
rq_pNEOs = np.array(rq_pNEOs_list)
rq_pNEOs_nonzero = rq_pNEOs[rq_pNEOs > 0]

n_sites_aneos = len(rq_pNEOs_nonzero)
sampled_aneos = np.random.lognormal(
    mean=mu_aneos, sigma=sigma_aneos,
    size=(n_sites_aneos, N_SIM)
)
sampled_aneos = np.maximum(sampled_aneos, 1.0)
rq_aneos_adjusted = rq_pNEOs_nonzero[:, np.newaxis] * sampled_aneos

col_medians_a = np.median(rq_aneos_adjusted, axis=0)
orig_median_a = np.median(rq_pNEOs_nonzero)
col_ratios_a  = col_medians_a / orig_median_a if orig_median_a > 0 else np.full(N_SIM, np.nan)

results['a-NEOs'] = {
    'rq_original':  rq_pNEOs_nonzero,
    'rq_adjusted':  rq_aneos_adjusted,
    'sim_medians':  col_medians_a,
    'sim_means':    np.mean(rq_aneos_adjusted, axis=0),
    'ratio_p2_5':     np.percentile(col_ratios_a, 2.5),
    'ratio_p50':    np.percentile(col_ratios_a, 50),
    'ratio_p97_5':    np.percentile(col_ratios_a, 97.5),
}
print(f"  a-NEOs: Completed {n_sites_aneos} sites × {N_SIM} simulation iterations")

# 5b. Site-aligned simulation for combined risk scenarios
print("\n  Simulating combined risk scenarios (aligned by site)...")

n_sites_global = len(rq_global)

rq_pNEOs_all = np.array([
    sum(rq_global[f'RQ-{c}'].iloc[i] for c in COMPOUNDS)
    for i in range(n_sites_global)
])

rq_pt_NEOs_2d = np.zeros((n_sites_global, N_SIM))

for cpd in COMPOUNDS:
    if rq_fit[cpd] is None:
        rq_cpd_all = rq_global[f'RQ-{cpd}'].values
        rq_pt_NEOs_2d += rq_cpd_all[:, np.newaxis]
        continue

    mu_mag, sigma_mag = fit_params[cpd]
    rq_cpd_all = rq_global[f'RQ-{cpd}'].values

    sampled_mag_cpd = np.random.lognormal(
        mean=mu_mag, sigma=sigma_mag,
        size=(n_sites_global, N_SIM)
    )
    sampled_mag_cpd = np.maximum(sampled_mag_cpd, 1.0)
    rq_pt_NEOs_2d += rq_cpd_all[:, np.newaxis] * sampled_mag_cpd

valid_mask = rq_pNEOs_all > 0
rq_pNEOs_valid   = rq_pNEOs_all[valid_mask]
rq_pt_NEOs_valid = rq_pt_NEOs_2d[valid_mask, :]
n_valid = valid_mask.sum()

sampled_aneos_aligned = np.random.lognormal(
    mean=mu_aneos, sigma=sigma_aneos,
    size=(n_valid, N_SIM)
)
sampled_aneos_aligned = np.maximum(sampled_aneos_aligned, 1.0)

rq_pa_NEOs_2d  = rq_pNEOs_valid[:, np.newaxis] * sampled_aneos_aligned
rq_pta_NEOs_2d = rq_pt_NEOs_valid * sampled_aneos_aligned

# Calculate the 95% simulation interval of the risk increment ratio for combined scenarios using the median p-NEO RQ as the reference
ref_median = np.median(rq_pNEOs_valid)

def calc_ratio_ci(arr_2d, ref_med, ci_lo=2.5, ci_hi=97.5):
    col_meds = np.median(arr_2d, axis=0)
    ratios = col_meds / ref_med if ref_med > 0 else np.full(arr_2d.shape[1], np.nan)
    return np.percentile(ratios, ci_lo), np.percentile(ratios, 50), np.percentile(ratios, ci_hi)

pt_r2_5,  pt_r50,  pt_r97_5  = calc_ratio_ci(rq_pt_NEOs_valid, ref_median)
pa_r2_5,  pa_r50,  pa_r97_5  = calc_ratio_ci(rq_pa_NEOs_2d,    ref_median)
pta_r2_5, pta_r50, pta_r97_5 = calc_ratio_ci(rq_pta_NEOs_2d,   ref_median)

results['p-NEOs']     = {'rq_all': rq_pNEOs_valid}
results['p+t-NEOs']   = {'rq_all': rq_pt_NEOs_valid,
                          'ratio_p2_5': pt_r2_5,  'ratio_p50': pt_r50,  'ratio_p97_5': pt_r97_5}
results['p+a-NEOs']   = {'rq_all': rq_pa_NEOs_2d,
                          'ratio_p2_5': pa_r2_5,  'ratio_p50': pa_r50,  'ratio_p97_5': pa_r97_5}
results['p+t+a-NEOs'] = {'rq_all': rq_pta_NEOs_2d,
                          'ratio_p2_5': pta_r2_5, 'ratio_p50': pta_r50, 'ratio_p97_5': pta_r97_5}

print(f"  Combined simulation completed: {n_valid} valid sites × {N_SIM} iterations")

# 6. Statistical summary
print("\n" + "=" * 60)
print("Statistical summary of Monte Carlo results:")

summary_rows = []
for key in COMPOUNDS + ['a-NEOs']:
    if key not in results:
        continue
    res = results[key]
    orig_median = np.median(res['rq_original'])
    all_adjusted = res['rq_adjusted'].flatten()
    adj_p5   = np.percentile(all_adjusted, 5)
    adj_p50  = np.percentile(all_adjusted, 50)
    adj_p95  = np.percentile(all_adjusted, 95)

    summary_rows.append({
        'Compound': key,
        'Original_RQ_P5': np.percentile(res['rq_original'], 5),
        'Original_RQ_median': orig_median,
        'Original_RQ_P95': np.percentile(res['rq_original'], 95),
        'Adjusted_RQ_P5': adj_p5,
        'Adjusted_RQ_median': adj_p50,
        'Adjusted_RQ_P95': adj_p95,
        'Risk_increment_ratio_median': res['ratio_p50'],
    })

summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv('MC_summary.csv', index=False, encoding='utf-8-sig')
print("Statistical summary saved: MC_summary.csv")

# Helper functions
def add_pdf_shadow(ax, data, color, alpha_fill=0.22, scale_to=0.28):
    d = data[data > 0]
    if len(d) < 5:
        return
    log_d = np.log10(d)
    kde = gaussian_kde(log_d, bw_method='scott')
    x_log = np.linspace(log_d.min() - 0.5, log_d.max() + 0.5, 500)
    density = kde(x_log)
    density_scaled = density / density.max() * scale_to
    ax.fill_between(10 ** x_log, 0, density_scaled,
                    color=color, alpha=alpha_fill, linewidth=0)


def smooth_cdf_logspace(data, n_points=800, bw=0.15):
    d = data[data > 0]
    log_d = np.log10(d)
    kde = gaussian_kde(log_d, bw_method=bw)
    x_log = np.linspace(log_d.min() - 0.3, log_d.max() + 0.3, n_points)
    density = kde(x_log)
    dx = x_log[1] - x_log[0]
    cdf = np.cumsum(density) * dx
    cdf /= cdf[-1]
    return 10 ** x_log, cdf


PALETTE = {
    'IMI': '#C0392B',
    'THM': '#1A5E9B',
    'DIN': '#27AE60',
    'ACE': '#D4870A',
    'CLO': '#6C3483',
    'NIP': '#E67E22',
}
PALETTE_LIST = [PALETTE[c] for c in COMPOUNDS]

# Figure 1: Risk increments from transformation products for individual compounds (six panels)
print("\nGenerating Figure 1: Risk increments from transformation products (six panels)...")
fig1, axes = plt.subplots(2, 3, figsize=(14, 9))
axes = axes.flatten()

for idx, cpd in enumerate(COMPOUNDS):
    if cpd not in results:
        axes[idx].set_visible(False)
        continue

    ax = axes[idx]
    res = results[cpd]
    color = PALETTE[cpd]

    add_pdf_shadow(ax, res['rq_original'], color=color,
                   alpha_fill=0.40, scale_to=0.28)

    x_o_sm, y_o_sm = smooth_cdf_logspace(res['rq_original'])
    ax.semilogx(x_o_sm, y_o_sm, color=color, lw=1.8,
                linestyle='--', alpha=0.85, label='p-NEOs only')

    n_flat = res['rq_adjusted'].size
    if n_flat > 50000:
        idx_sample = np.random.choice(n_flat, 50000, replace=False)
        adj_flat = res['rq_adjusted'].flatten()[idx_sample]
    else:
        adj_flat = res['rq_adjusted'].flatten()

    x_adj_sm, y_adj_sm = smooth_cdf_logspace(adj_flat)
    ax.semilogx(x_adj_sm, y_adj_sm, color=color, lw=2.2,
                label='p+t-NEOs (MC)')

    ax.fill_between(x_o_sm, 0, y_o_sm, color=color, alpha=0.08)

    med_orig = np.median(res['rq_original'])
    med_adj  = np.median(adj_flat)
    ax.axvline(med_orig, color=color, lw=1.0, linestyle=':', alpha=0.7)
    ax.axvline(med_adj,  color=color, lw=1.5, linestyle='-', alpha=0.9)
    ax.axvline(1.0, color='gray', lw=0.8, linestyle='--', alpha=0.4)

    r_p50 = res['ratio_p50']
    r_p2_5  = res['ratio_p2_5']
    r_p97_5 = res['ratio_p97_5']

    ax.text(0.97, 0.05,
            f'Median ratio: {r_p50:.2f}×\n(95% simulation interval: {r_p2_5:.2f}–{r_p97_5:.2f}×)',
            transform=ax.transAxes,
            ha='right', va='bottom',
            fontsize=9, color=color,
            bbox=dict(boxstyle='round,pad=0.3',
                      facecolor='white', alpha=0.7, edgecolor=color))

    ax.set_xlabel('RQ', fontsize=10)
    ax.set_ylabel('Cumulative probability', fontsize=10)
    ax.set_title(f'{cpd}', fontsize=12, fontweight='bold', color=color)
    ax.set_ylim(0, 1.05)
    ax.set_xlim(left=1e-5)
    ax.legend(fontsize=8, loc='upper left')
    ax.grid(True, which='both', alpha=0.2)
    ax.tick_params(labelsize=9)

fig1.suptitle(
    'Monte Carlo Simulation: Risk Increment from Transformation Products\n'
    '(p-NEOs → p+t-NEOs, N=10,000)',
    fontsize=13, fontweight='bold', y=1.01
)
plt.tight_layout()
plt.savefig('MC_CDF_transformation_products.pdf',
            dpi=300, bbox_inches='tight')
plt.close()
print("  Saved: MC_CDF_transformation_products.pdf")

# Figure 2: Risk increment from structural analogues (single panel)
print("Generating Figure 2: Risk increment from structural analogues...")
fig2, ax2 = plt.subplots(1, 1, figsize=(7, 5.5))

res_a = results['a-NEOs']
color_orig = '#A93226'
color_adj  = '#1A5276'

add_pdf_shadow(ax2, res_a['rq_original'],
               color=color_orig, alpha_fill=0.40, scale_to=0.28)
add_pdf_shadow(ax2, res_a['rq_adjusted'].flatten(),
               color=color_adj,  alpha_fill=0.30, scale_to=0.28)

x_o_sm, y_o_sm = smooth_cdf_logspace(res_a['rq_original'])
ax2.semilogx(x_o_sm, y_o_sm, color=color_orig, lw=2.0,
             linestyle='--', label='p-NEOs (original)')

n_flat = res_a['rq_adjusted'].size
if n_flat > 50000:
    idx_sample = np.random.choice(n_flat, 50000, replace=False)
    adj_flat = res_a['rq_adjusted'].flatten()[idx_sample]
else:
    adj_flat = res_a['rq_adjusted'].flatten()

x_adj_sm, y_adj_sm = smooth_cdf_logspace(adj_flat)
ax2.semilogx(x_adj_sm, y_adj_sm, color=color_adj, lw=2.5,
             label='p+a-NEOs (MC adjusted)')

ax2.fill_between(x_o_sm,   0, y_o_sm,   color=color_orig, alpha=0.08)
ax2.fill_between(x_adj_sm, 0, y_adj_sm, color=color_adj,  alpha=0.08)

med_orig = np.median(res_a['rq_original'])
med_adj  = np.median(adj_flat)
ax2.axvline(med_orig, color=color_orig, lw=1.5, linestyle='--',
            alpha=0.85, label=f'Median (original): {med_orig:.3f}')
ax2.axvline(med_adj,  color=color_adj,  lw=1.5, linestyle='-',
            alpha=0.85, label=f'Median (adjusted): {med_adj:.3f}')
ax2.axvline(1.0, color='gray', lw=1.2, linestyle=':', alpha=0.5,
            label='RQ = 1')

r_p50_a = res_a['ratio_p50']
r_p2_5_a  = res_a['ratio_p2_5']
r_p97_5_a = res_a['ratio_p97_5']

ax2.text(0.97, 0.05,
         f'Median risk increment: {r_p50_a:.2f}×\n(95% simulation interval: {r_p2_5_a:.2f}–{r_p97_5_a:.2f}×)',
         transform=ax2.transAxes,
         ha='right', va='bottom', fontsize=11,
         bbox=dict(boxstyle='round,pad=0.4',
                   facecolor='lightyellow', alpha=0.9,
                   edgecolor=color_adj))

ax2.set_xlabel('RQ (Σp-NEOs or Σp+a-NEOs)', fontsize=11)
ax2.set_ylabel('Cumulative probability', fontsize=11)
ax2.set_title(
    'Monte Carlo Simulation: Risk Increment from Structural Analogues\n'
    '(p-NEOs → p+a-NEOs, N=10,000)',
    fontsize=11, fontweight='bold'
)
ax2.set_ylim(0, 1.05)
ax2.legend(fontsize=9, loc='upper left')
ax2.grid(True, which='both', alpha=0.2)
plt.tight_layout()
plt.savefig('MC_CDF_analogues.pdf', dpi=300, bbox_inches='tight')
plt.close()
print("  Saved: MC_CDF_analogues.pdf")

# Figure 3: Comparison across all compounds
print("Generating Figure 3: Comparison across all compounds...")
fig3, ax3 = plt.subplots(1, 1, figsize=(9, 6.5))

for i, cpd in enumerate(COMPOUNDS):
    if cpd not in results:
        continue
    res = results[cpd]
    color = PALETTE_LIST[i]

    add_pdf_shadow(ax3, res['rq_original'], color=color,
                   alpha_fill=0.40, scale_to=0.22)

    n_flat = res['rq_adjusted'].size
    if n_flat > 30000:
        idx_s = np.random.choice(n_flat, 30000, replace=False)
        adj_flat = res['rq_adjusted'].flatten()[idx_s]
    else:
        adj_flat = res['rq_adjusted'].flatten()

    x_adj_sm, y_adj_sm = smooth_cdf_logspace(adj_flat)
    ax3.semilogx(x_adj_sm, y_adj_sm, color=color, lw=2.2,
                 label=f'{cpd}')

    x_o_sm, y_o_sm = smooth_cdf_logspace(res['rq_original'])
    ax3.semilogx(x_o_sm, y_o_sm, color=color, lw=1.4,
                 linestyle='--', alpha=0.65)

ax3.axvline(1.0, color='gray', lw=1.0, linestyle=':', alpha=0.6)

# Add a summary box at the lower right of Figure 3 showing the 95% simulation interval for each compound
si_text_lines = ['95% simulation interval of median ratio:']
for cpd in COMPOUNDS:
    if cpd not in results:
        continue
    res = results[cpd]
    si_text_lines.append(
        f'  {cpd}: {res["ratio_p50"]:.2f}× ({res["ratio_p2_5"]:.2f}–{res["ratio_p97_5"]:.2f}×)'
    )
ax3.text(0.97, 0.03, '\n'.join(si_text_lines),
         transform=ax3.transAxes,
         ha='right', va='bottom', fontsize=7.5,
         fontfamily='monospace',
         bbox=dict(boxstyle='round,pad=0.4',
                   facecolor='white', alpha=0.85,
                   edgecolor='#AAAAAA'))

legend_lines = [
    Line2D([0], [0], color=PALETTE_LIST[i], lw=2.2, label=COMPOUNDS[i])
    for i in range(len(COMPOUNDS))
]
legend_style = [
    Line2D([0], [0], color='gray', lw=2.2, label='p+t-NEOs (solid)'),
    Line2D([0], [0], color='gray', lw=1.4, linestyle='--',
           label='p-NEOs only (dashed)'),
]
ax3.legend(handles=legend_lines + legend_style,
           fontsize=9, ncol=2, loc='upper left')
ax3.set_xlabel('RQ', fontsize=12)
ax3.set_ylabel('Cumulative probability', fontsize=12)
ax3.set_title(
    'Monte Carlo Simulation: Global Surface Water Neonicotinoid Risk\n'
    'Solid = p+t-NEOs adjusted  |  Dashed = p-NEOs only',
    fontsize=11, fontweight='bold'
)
ax3.set_ylim(0, 1.05)
ax3.grid(True, which='both', alpha=0.2)
plt.tight_layout()
plt.savefig('MC_CDF_all_compounds.pdf', dpi=300, bbox_inches='tight')
plt.close()
print("  Saved: MC_CDF_all_compounds.pdf")

# Figure 4: Comparison of four combined risk scenarios
print("\nGenerating Figure 4: Comparison of four combined risk scenarios...")

np.random.seed(42)

C_pNEOs = '#C0392B'
C_pt    = '#1A5276'
C_pa    = '#27AE60'
C_pta   = '#6C3483'

fig4, ax4 = plt.subplots(1, 1, figsize=(8, 6))

def get_flat_sample(arr_2d, max_pts=80000):
    flat = arr_2d.flatten()
    if len(flat) > max_pts:
        flat = flat[np.random.choice(len(flat), max_pts, replace=False)]
    return flat

# 1. p-NEOs (dashed line)
pNEOs_data = results['p-NEOs']['rq_all']
add_pdf_shadow(ax4, pNEOs_data, color=C_pNEOs, alpha_fill=0.35, scale_to=0.25)
x_pn, y_pn = smooth_cdf_logspace(pNEOs_data)
ax4.semilogx(x_pn, y_pn, color=C_pNEOs, lw=2.0, linestyle='--',
             label='p-NEOs', zorder=4)
med_pn = np.median(pNEOs_data)

# 2. p+t-NEOs
pt_flat = get_flat_sample(results['p+t-NEOs']['rq_all'])
add_pdf_shadow(ax4, pt_flat, color=C_pt, alpha_fill=0.30, scale_to=0.25)
x_pt, y_pt = smooth_cdf_logspace(pt_flat)
ax4.semilogx(x_pt, y_pt, color=C_pt, lw=2.3, linestyle='-',
             label='p+t-NEOs', zorder=3)
med_pt = np.median(pt_flat)

# 3. p+a-NEOs
pa_flat = get_flat_sample(results['p+a-NEOs']['rq_all'])
add_pdf_shadow(ax4, pa_flat, color=C_pa, alpha_fill=0.30, scale_to=0.25)
x_pa, y_pa = smooth_cdf_logspace(pa_flat)
ax4.semilogx(x_pa, y_pa, color=C_pa, lw=2.3, linestyle='-',
             label='p+a-NEOs', zorder=3)
med_pa = np.median(pa_flat)

# 4. p+t+a-NEOs
pta_flat = get_flat_sample(results['p+t+a-NEOs']['rq_all'])
add_pdf_shadow(ax4, pta_flat, color=C_pta, alpha_fill=0.35, scale_to=0.25)
x_pta, y_pta = smooth_cdf_logspace(pta_flat)
ax4.semilogx(x_pta, y_pta, color=C_pta, lw=2.8, linestyle='-',
             label='p+t+a-NEOs', zorder=5)
med_pta = np.median(pta_flat)

ax4.axvline(1.0, color='gray', lw=1.0, linestyle=':', alpha=0.55,
            label='RQ = 1')

# Summary box: median RQ, median risk increment ratio, and 95% simulation interval
pt_r2_5,  pt_r50,  pt_r97_5  = results['p+t-NEOs']['ratio_p2_5'],  results['p+t-NEOs']['ratio_p50'],  results['p+t-NEOs']['ratio_p97_5']
pa_r2_5,  pa_r50,  pa_r97_5  = results['p+a-NEOs']['ratio_p2_5'],  results['p+a-NEOs']['ratio_p50'],  results['p+a-NEOs']['ratio_p97_5']
pta_r2_5, pta_r50, pta_r97_5 = results['p+t+a-NEOs']['ratio_p2_5'], results['p+t+a-NEOs']['ratio_p50'], results['p+t+a-NEOs']['ratio_p97_5']

info_text = (
    f"Median (p-NEOs):      {med_pn:.2f}\n"
    f"p+t-NEOs:    {med_pt:.2f}  (×{pt_r50:.2f}, 95% simulation interval: {pt_r2_5:.2f}–{pt_r97_5:.2f})\n"
    f"p+a-NEOs:    {med_pa:.2f}  (×{pa_r50:.2f}, 95% simulation interval: {pa_r2_5:.2f}–{pa_r97_5:.2f})\n"
    f"p+t+a-NEOs:  {med_pta:.2f}  (×{pta_r50:.2f}, 95% simulation interval: {pta_r2_5:.2f}–{pta_r97_5:.2f})"
)
ax4.text(0.02, 0.97, info_text,
         transform=ax4.transAxes,
         ha='left', va='top', fontsize=8.5,
         fontfamily='monospace',
         bbox=dict(boxstyle='round,pad=0.5',
                   facecolor='white', alpha=0.85,
                   edgecolor='#AAAAAA'))

ax4.set_xlabel('Risk quotients', fontsize=12)
ax4.set_ylabel('Cumulative probability', fontsize=12)
ax4.set_title(
    'Monte Carlo Simulation: Cumulative Risk Increment\n'
    'p-NEOs → p+t-NEOs / p+a-NEOs / p+t+a-NEOs  (N=10,000)',
    fontsize=11, fontweight='bold'
)
ax4.set_ylim(0, 1.05)
ax4.legend(fontsize=10, loc='upper left',
           framealpha=0.9, edgecolor='#CCCCCC')
ax4.grid(True, which='both', alpha=0.2)
ax4.tick_params(labelsize=10)

plt.tight_layout()
plt.savefig('MC_CDF_combined_4curves.pdf', dpi=300, bbox_inches='tight')
plt.close()
print("  Saved: MC_CDF_combined_4curves.pdf")

# Table 4: Percentile summary of the four combined risk scenarios
print("\nGenerating Table 4: Percentile summary of the four combined risk scenarios...")

percentiles_list = [5, 10, 25, 50, 75, 90, 95]

curve_data = {
    'p-NEOs':     pNEOs_data,
    'p+t-NEOs':   pt_flat,
    'p+a-NEOs':   pa_flat,
    'p+t+a-NEOs': pta_flat,
}

table4_rows = []
for curve_name, data in curve_data.items():
    d = data[data > 0]
    row = {'Curve': curve_name}
    for p in percentiles_list:
        row[f'P{p}'] = np.percentile(d, p)
    row['Mean'] = np.mean(d)
    ref_median = np.median(pNEOs_data[pNEOs_data > 0])
    row['Median_ratio_vs_pNEOs'] = row['P50'] / ref_median if ref_median > 0 else np.nan
    table4_rows.append(row)

table4_df = pd.DataFrame(table4_rows)
table4_df.to_csv('MC_combined_4curves_summary.csv',
                 index=False, encoding='utf-8-sig')
print("  Saved: MC_combined_4curves_summary.csv")

# Generate MC_ratio_SI.csv
print("\nGenerating MC_ratio_SI.csv (95% simulation interval: P2.5–P97.5)...")

ratio_si_rows = []

# Individual compounds (t-NEOs)
for cpd in COMPOUNDS:
    if cpd not in results:
        continue
    res = results[cpd]
    ratio_si_rows.append({
        'Compound_or_Curve':       cpd,
        'Category':                't-NEOs (transformation products)',
        'Original_RQ_median':      np.median(res['rq_original']),
        'Adjusted_RQ_median':      np.median(res['rq_adjusted'].flatten()),
        'Ratio_P50 (median)':      res['ratio_p50'],
        'Ratio_P2.5  (95% simulation interval_lower)': res['ratio_p2_5'],
        'Ratio_P97.5 (95% simulation interval_upper)': res['ratio_p97_5'],
    })

# a-NEOs
res_a = results['a-NEOs']
ratio_si_rows.append({
    'Compound_or_Curve':       'a-NEOs',
    'Category':                'a-NEOs (structural analogues)',
    'Original_RQ_median':      np.median(res_a['rq_original']),
    'Adjusted_RQ_median':      np.median(res_a['rq_adjusted'].flatten()),
    'Ratio_P50 (median)':      res_a['ratio_p50'],
    'Ratio_P2.5  (95% simulation interval_lower)': res_a['ratio_p2_5'],
    'Ratio_P97.5 (95% simulation interval_upper)': res_a['ratio_p97_5'],
})

# Combined risk scenarios
for curve_key, label in [
    ('p+t-NEOs',   'p+t-NEOs (combined)'),
    ('p+a-NEOs',   'p+a-NEOs (combined)'),
    ('p+t+a-NEOs', 'p+t+a-NEOs (combined)'),
]:
    res_c = results[curve_key]
    ratio_si_rows.append({
        'Compound_or_Curve':       curve_key,
        'Category':                label,
        'Original_RQ_median':      ref_median,
        'Adjusted_RQ_median':      np.median(res_c['rq_all'].flatten()),
        'Ratio_P50 (median)':      res_c['ratio_p50'],
        'Ratio_P2.5  (95% simulation interval_lower)': res_c['ratio_p2_5'],
        'Ratio_P97.5 (95% simulation interval_upper)': res_c['ratio_p97_5'],
    })

ratio_si_df = pd.DataFrame(ratio_si_rows)
ratio_si_df.to_csv('MC_ratio_SI.csv', index=False, encoding='utf-8-sig')

print("\n" + "=" * 85)
print(f"{'Compound/Curve':<18} {'Orig_med':>10} {'Adj_med':>10} "
      f"{'Ratio_P2.5':>10} {'Ratio_P50':>10} {'Ratio_P97.5':>10}")
print("-" * 85)
for _, row in ratio_si_df.iterrows():
    print(f"  {row['Compound_or_Curve']:<16} "
          f"{row['Original_RQ_median']:>10.4f} "
          f"{row['Adjusted_RQ_median']:>10.4f} "
          f"{row['Ratio_P2.5  (95% simulation interval_lower)']:>10.3f} "
          f"{row['Ratio_P50 (median)']:>10.3f} "
          f"{row['Ratio_P97.5 (95% simulation interval_upper)']:>10.3f}")
print("=" * 85)
print("  Saved: MC_ratio_SI.csv")

# 8. Detailed percentile statistics
percentiles = [5, 10, 25, 50, 75, 90, 95]
detail_rows = []

for key in COMPOUNDS + ['a-NEOs']:
    if key not in results:
        continue
    res = results[key]
    adj_all = res['rq_adjusted'].flatten()
    row = {'Compound': key}
    for p in percentiles:
        row[f'Original_RQ_P{p}'] = np.percentile(res['rq_original'], p)
        row[f'Adjusted_RQ_P{p}'] = np.percentile(adj_all, p)
    detail_rows.append(row)

detail_df = pd.DataFrame(detail_rows)
detail_df.to_csv('MC_percentiles.csv', index=False, encoding='utf-8-sig')
print("  Saved: MC_percentiles.csv")

print("\n" + "=" * 60)
print("Monte Carlo simulation completed.")
print("Output files:")
print("  MC_CDF_transformation_products.pdf")
print("  MC_CDF_analogues.pdf")
print("  MC_CDF_all_compounds.pdf")
print("  MC_CDF_combined_4curves.pdf")
print("  MC_summary.csv")
print("  MC_percentiles.csv")
print("  MC_combined_4curves_summary.csv")
print("  MC_ratio_SI.csv")
print("=" * 60)