from ema.annotate.gtftobed import gtf_bed as gb
from ema.config import directory_config as dc
gb(gtfdir=dc.gtf_dir, endbeddir=dc.endbed, featuresdir=dc.raw_features)

