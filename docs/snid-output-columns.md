# SNID output columns: type, subtype and template

This page explains how the SNID-SAGE type columns in the bosscat results CSV are computed, and why they can disagree. Based on SNID-SAGE 1.2.5.

## How SNID-SAGE reaches a classification

1. The spectrum is correlated against every template. Each match has a redshift and a quality metric (HσLAP-CCC).
2. Matches are grouped by template type (Ia, Ib, Ic, II, Galaxy, ...). Within each type they are clustered in redshift.
3. Each cluster is scored by the mean of its five best HσLAP-CCC values. The top-scoring cluster is the **winning cluster**.

## Columns

| Column | Meaning | Source |
|---|---|---|
| `snid_consensus_type` | Broad type of the winning cluster (Ia, Ib, Ic, II, ...). `Unknown` if no match survives. | Winning cluster |
| `snid_best_subtype` | Most supported subtype in the winning cluster (for example `Ia-norm`, `IIn`). | Vote over all cluster members |
| `snid_template_type`, `snid_template_subtype`, `snid_template_name` | One template from the winning cluster. | First entry of the cluster's match list |
| `snid_redshift` | Redshift of that same template. | First entry of the cluster's match list |
| `snid_consensus_redshift` | Weighted redshift of the winning cluster. | Winning cluster |

## How `snid_best_subtype` is computed

For the members of the winning cluster, matches are grouped by subtype. Each subtype scores the mean of its five best HσLAP-CCC values. If it has fewer than five templates, the score is multiplied by n/5. The highest score wins.

This is a consensus over the whole cluster (over a thousand matches for a typical Ia).

## How the template columns are computed

The template columns describe a single template: the first entry of the winning cluster's match list. In SNID-SAGE 1.2.5 that list is ordered by redshift within the cluster, so the first entry is the **lowest-redshift member**, not the best-correlating match.

In three test spectra the template taken had a clearly lower metric than the best match in the cluster:

| Spectrum | `snid_best_subtype` | Template taken | Best-metric match |
|---|---|---|---|
| 100050 | Ia-norm | Ia-03fg (metric 5.6, z 0.067) | Ia-norm (metric 11.1, z 0.078) |
| 100407 | Ia-norm | Ia-02cx (metric 2.2, z 0.126) | Ia-91T (metric 6.3, z 0.156) |
| 100157 | Gal-SB | Gal-SB (metric 11.1) | Gal-SB (metric 17.4) |

This was observed on a small sample. It has not been confirmed as intended upstream behavior.

## Why the subtype columns can disagree

1. **Different statistics.** The subtype is a top-5-mean vote across the cluster. The template is one match.
2. **Ordering.** The template is chosen by redshift order, not by quality, so it can be a mediocre match.
3. **Real ambiguity.** When two subtypes score similarly in the vote, either can win. SNID-SAGE computes a second-best subtype and a margin, but bosscat does not currently report them.

## Practical guidance

- Treat `snid_consensus_type` as the headline classification and `snid_best_subtype` as the specific one.
- A disagreement between `snid_best_subtype` and `snid_template_subtype` is not by itself evidence that the transient is ambiguous.
- Prefer `snid_consensus_redshift` over `snid_redshift` for science, since `snid_redshift` comes from the lowest-redshift cluster member.
- Agreement of `snid_consensus_type` with other classifiers (NGSF, TNS) is evaluated at the broad-class level. Subtype labels differ between tools.
