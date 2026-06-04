"""spec-collapse-atlas: correct aggregation for multiverse / many-analyst meta-analysis."""
from .datasets import DATASETS
from .engine import enumerate_specs
from .aggregators import naive_concordance, naive_ivre_pool, weighted_likelihood
from .coverage import run_coverage

__all__ = [
    "DATASETS", "enumerate_specs",
    "naive_concordance", "naive_ivre_pool", "weighted_likelihood",
    "run_coverage",
]
