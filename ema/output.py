import os
import json
from datetime import datetime


class OutputManager:
    """Manages hierarchical output directory structure for PeakATail."""

    def __init__(self, base_dir="emaout"):
        self.base_dir = base_dir
        self.dirs = {
            "peak_calling": os.path.join(base_dir, "01_peak_calling"),
            "cb_filter": os.path.join(base_dir, "02_cb_filter"),
            "gtf_annotation": os.path.join(base_dir, "03_gtf_annotation"),
            "pas_gene": os.path.join(base_dir, "04_pas_gene_assignment"),
            "annotated": os.path.join(base_dir, "05_annotated_matrix"),
            "preprocessing": os.path.join(base_dir, "06_preprocessing"),
            "clustering": os.path.join(base_dir, "07_clustering"),
            "differential": os.path.join(base_dir, "08_differential"),
            "gtf_cache": os.path.join(base_dir, "gtf_cache"),
        }

    def setup(self):
        """Create all output directories."""
        for d in self.dirs.values():
            os.makedirs(d, exist_ok=True)

    def path(self, stage, filename):
        """Get path for a file in a specific stage directory."""
        return os.path.join(self.dirs[stage], filename)

    def save_stats(self, stage, stats):
        """Save statistics JSON for a pipeline stage."""
        path = self.path(stage, f"{stage}_stats.json")
        with open(path, 'w') as f:
            json.dump(stats, f, indent=2)

    def save_run_config(self, args_dict):
        """Save full run configuration."""
        config = {
            "timestamp": datetime.now().isoformat(),
            **args_dict
        }
        path = os.path.join(self.base_dir, "run_config.json")
        with open(path, 'w') as f:
            json.dump(config, f, indent=2)
