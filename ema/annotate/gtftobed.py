from collections import defaultdict
from ema.config import directory_config as dc


def parse_gtf_attributes(attr_string):
    """Parse GTF attribute column (col 8) into a key-value dict.

    GTF attributes are semicolon-delimited key-value pairs like:
        gene_id "ENSG00000223972"; gene_type "transcribed_unprocessed_pseudogene";
    """
    attrs = {}
    for item in attr_string.split(';'):
        item = item.strip()
        if not item:
            continue
        if ' ' in item:
            key, val = item.split(' ', 1)
            attrs[key] = val.strip('"')
    return attrs


def gtf_bed(endbeddir=dc.endbed,
            gtfdir=dc.gtf_dir,
            featuresdir=dc.raw_features,
            utr_lengths_dir=None,
            biotypes=("protein_coding", "lncRNA", "snRNA", "antisense",
                      "miRNA", "processed_transcript", "lincRNA"),
            source_types=("gene",),
            id_prefixes=("ENSG",),
            score=0
            ) -> dict:
    """Convert a GTF file to BED format, extract features, and compute UTR lengths.

    Parses gene features for BED output and three_prime_utr features to compute
    per-gene UTR lengths for adaptive distance thresholding in find_close.

    Args:
        endbeddir: Path to write the gene BED file.
        gtfdir: Path to the input GTF file.
        featuresdir: Path to write the extracted gene features (id, symbol).
        utr_lengths_dir: Path to write per-gene UTR lengths TSV. If None, uses
            output_dir/utr_lengths.tsv.
        biotypes: Tuple of gene biotypes to include.
        source_types: Tuple of GTF feature types for gene entries.
        id_prefixes: Tuple of gene ID prefixes to include (e.g., "ENSG").
        score: BED score field value.

    Returns:
        dict: Mapping of gene_id -> max 3'UTR length (bp) for that gene.
    """
    if utr_lengths_dir is None:
        import os
        utr_lengths_dir = os.path.join(os.path.dirname(endbeddir), "utr_lengths.tsv")

    # Collect 3'UTR lengths per gene (a gene can have multiple UTR annotations)
    utr_lengths = defaultdict(int)

    with open(gtfdir, "r") as gtf, open(endbeddir, "w") as endbed, \
         open(featuresdir, "w") as features:

        for line in gtf:
            if line.startswith("#"):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            chro = fields[0]
            source = fields[2]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]

            # Parse attributes from column 8 (key-value, not positional)
            attrs = parse_gtf_attributes(fields[8])
            gene_id = attrs.get("gene_id", "")
            gene_biotype = attrs.get("gene_type", attrs.get("gene_biotype", ""))
            gene_name = attrs.get("gene_name", attrs.get("gene_symbol", gene_id))

            # Extract 3'UTR features for adaptive distance thresholding
            if source == "three_prime_utr" or source == "three_prime_UTR":
                utr_gene_id = gene_id.split('.')[0]  # strip version
                utr_length = end - start
                # Keep the longest UTR per gene
                if utr_length > utr_lengths[utr_gene_id]:
                    utr_lengths[utr_gene_id] = utr_length

            # Process gene features for BED output
            if source in source_types:
                # Strip version suffix for matching (ENSG00000223972.5 -> ENSG00000223972)
                gene_id_base = gene_id.split('.')[0]

                # Check prefix match
                prefix_match = any(gene_id_base.startswith(p) for p in id_prefixes)
                if not prefix_match or gene_biotype not in biotypes:
                    continue

                # Extend gene boundaries for nearby PAS detection
                match strand:
                    case "+":
                        end += 5000

                    case "-":
                        if start > 5000:
                            start -= 5000
                        else:
                            start = 1

                bed_info = f"{chro}\t{start}\t{end}\t{gene_id_base}\t{gene_name}\t{strand}\n"
                endbed.write(bed_info)

                feature = f"{gene_id_base}\t{gene_name}\n"
                features.write(feature)

    # Write UTR lengths to file
    with open(utr_lengths_dir, "w") as utr_file:
        utr_file.write("gene_id\tutr_length\n")
        for gene_id, length in sorted(utr_lengths.items()):
            utr_file.write(f"{gene_id}\t{length}\n")

    return dict(utr_lengths)


if __name__ == "__main__":
    gtf_bed()
