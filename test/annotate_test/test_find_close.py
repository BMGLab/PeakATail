from ema.annotate.find_close import find_close
from ema.config import directory_config


close = find_close(posbed_dir=directory_config.posbed, 
            negbe_dir=directory_config.negbed,
            genomebed_dir=directory_config.endbed,
            mergebed=directory_config.pasbed,
            features=directory_config.pas_geneid,
            annotatedbed_dir=directory_config.annotatedbed)


print(close)