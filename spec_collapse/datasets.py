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
    "warfarin": {
        # Canonical metadat dat.hart1999: warfarin vs control for stroke
        # prevention in atrial fibrillation, 6 trials, log incidence-rate ratio
        # from published strokes / person-years.
        "name": "Warfarin & Stroke Prevention (AF)",
        "source": "Hart et al. (1999); metadat dat.hart1999",
        "measure": "log-IRR",
        "yi": [-0.7842, -0.9359, -1.5793, -0.3887, -1.2019, -1.1409],
        "vi": [0.1637, 0.1776, 0.4103, 0.2778, 0.1863, 0.07],
    },
    "statins": {
        # Canonical metadat dat.cannon2006: intensive vs moderate statin therapy,
        # 4 trials, log-OR for coronary death or non-fatal MI.
        "name": "Intensive vs Moderate Statins (CHD events)",
        "source": "Cannon et al. (2006); metadat dat.cannon2006",
        "measure": "log-OR",
        "yi": [-0.1888, -0.1676, -0.2401, -0.1296],
        "vi": [0.0137, 0.0101, 0.0058, 0.0051],
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
    "sdd": {
        # Canonical metadat dat.damico2009: selective digestive decontamination
        # (topical + systemic antibiotics) vs control for respiratory tract
        # infections, 16 trials, log-OR; Jacobs 1992 (zero cell) uses the
        # Haldane-Anscombe 0.5 correction.
        "name": "Antibiotics & Respiratory Infections (SDD)",
        "source": "D'Amico et al. (2009); metadat dat.damico2009",
        "measure": "log-OR",
        "yi": [-2.4313, -3.2291, -1.2738, -0.3765, -1.2182, -0.499, -2.2654,
               -2.8362, -1.4271, -1.5945, -0.8417, -0.8073, -1.7494, -0.7628,
               -1.8944, -0.9248],
        "vi": [0.2855, 1.1038, 0.1239, 0.2525, 0.3633, 0.5142, 2.2677, 0.3175,
               0.2067, 0.2423, 0.0705, 0.0433, 0.2344, 0.0832, 0.4169, 0.0323],
    },
}
