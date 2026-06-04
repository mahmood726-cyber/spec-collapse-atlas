"""Built-in meta-analytic datasets (yi, vi pairs).

Lifted verbatim from MultiverseMA (C:\\Projects\\MultiverseMA\\multiverse-ma.html,
BUILT_IN_DATASETS) so the spec-collapse prototype operates on exactly the same
numbers the existing engine ships. yi = effect size (log-RR or log-OR), vi =
within-study variance.
"""

DATASETS = {
    "bcg": {
        # Canonical metafor/metadat dat.bcg, derived from the published Colditz
        # (1994) 2x2 counts. Engine reproduces metafor's published RE result on
        # these exactly: REML tau^2 = 0.3132, estimate = -0.7145 (see tests).
        # NB: MultiverseMA shipped a corrupted version of this dataset (9/13 rows
        # wrong); these authoritative values replace it.
        "name": "BCG Vaccine & Tuberculosis",
        "source": "Colditz et al. (1994); metadat dat.bcg",
        "measure": "log-RR",
        "yi": [-0.889311, -1.585389, -1.348073, -1.441551, -0.217547, -0.786116,
               -1.620898, 0.011952, -0.469418, -1.371345, -0.339359, 0.445913,
               -0.017314],
        "vi": [0.325585, 0.194581, 0.415368, 0.02001, 0.05121, 0.006906,
               0.223017, 0.003962, 0.056434, 0.073025, 0.012412, 0.532506,
               0.071405],
    },
    "aspirin": {
        "name": "Aspirin & Stroke Prevention",
        # Illustrative compilation (trial labels are real but the per-study
        # yi/vi are not verified against a single published dataset). Used only
        # as an interactive demo; the corpus result does not depend on it.
        "source": "Illustrative compilation (not source-verified)",
        "measure": "log-OR",
        "yi": [-0.29, -0.18, -0.34, 0.03, -0.02, -0.19, -0.08, -0.23],
        "vi": [0.0324, 0.0196, 0.0256, 0.0004, 0.0009, 0.0289, 0.0144, 0.0225],
    },
    "omega3": {
        "name": "Omega-3 & Cardiovascular Mortality",
        # Illustrative compilation (real trial labels, per-study yi/vi not
        # source-verified). Interactive demo only.
        "source": "Illustrative compilation (not source-verified)",
        "measure": "log-RR",
        "yi": [-0.20, -0.19, -0.01, -0.02, -0.04, -0.20, 0.01],
        "vi": [0.0064, 0.0144, 0.0049, 0.0064, 0.0100, 0.0100, 0.0081],
    },
    "magnesium": {
        # Canonical metadat dat.egger2001 (16 IV-magnesium-in-AMI trials),
        # log-OR from the published 2x2 counts; Bertschat 1989 (zero cell) uses
        # the Haldane-Anscombe 0.5 correction. Replaces MV's non-canonical
        # 10-trial compilation. The LIMIT-2 / ISIS-4 contrast is the classic
        # fragility controversy.
        "name": "IV Magnesium & Mortality in AMI",
        "source": "Egger et al. (2001); metadat dat.egger2001",
        "measure": "log-OR",
        "yi": [-0.8303, -1.0561, -1.2783, -0.0435, 0.2231, -2.4075, -1.2809,
               -1.1917, -0.6957, -2.2083, -2.0382, -0.8502, -0.7932, -0.2993,
               -1.5708, 0.0576],
        "vi": [1.5551, 0.1715, 0.6531, 2.0435, 0.2393, 1.1496, 1.425, 2.7599,
               0.2875, 1.2313, 0.6095, 0.3825, 0.3917, 0.0215, 0.3295, 0.001],
    },
    "corticosteroids": {
        "name": "Corticosteroids for Preterm Birth",
        # Illustrative compilation in the spirit of the antenatal-corticosteroids
        # review (real trial labels); per-study yi/vi not verified against the
        # published dataset, so not cited as verbatim. Interactive demo only.
        "source": "Illustrative compilation (not source-verified)",
        "measure": "log-OR",
        "yi": [-0.80, -0.54, -0.59, -0.13, -0.94, -0.35, -0.48, -0.27, -0.23,
               -0.31, -0.33, -0.22],
        "vi": [0.1681, 0.5329, 0.3721, 0.1764, 0.2116, 0.3136, 0.5776, 0.0576,
               0.0841, 0.1024, 0.0441, 0.0289],
    },
}
