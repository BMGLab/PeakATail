def gtf_3utr_bed(bedfile, 
                 gtffile,
                 biotypes=("protein_coding", "lncRNA", "snRNA", "antisense", "miRNA", "processed_transcript", "lincRNA"), 
                 feature_type="three_prime_utr",
                 id_head=("ENST",),
                 chro_index=0,
                 feature_index=2,
                 start_index=3,
                 end_index=4,
                 strand_index=6,
                 biotype_key="gene_biotype",
                 transcript_id_key="transcript_id",
                 gene_name_key="gene_name",
                 score=0,
                 extension=5000  # <<< ADD extension amount as a parameter for flexibility
                 ) -> None:
    """
    Extract 3'-UTR regions from a GTF file and write to BED format.
    """

    def parse_gtf_attributes(attribute_text):
        """Helper to parse GTF attributes into a dictionary."""
        attributes = {}
        fields = attribute_text.strip().split(';')
        for field in fields:
            if field.strip():
                key, value = field.strip().split(' ', 1)
                attributes[key] = value.replace('"', '').replace(';', '')
        return attributes

    with open(gtffile, "r") as gtf, open(bedfile, "w") as bedout:
        for line in gtf:
            if line.startswith("#"):
                continue

            info = line.strip().split('\t')

            if len(info) < 9:
                continue  # skip incomplete lines

            chro = info[chro_index]
            feature = info[feature_index]
            start = int(info[start_index]) - 1  # BED is 0-based
            end = int(info[end_index])          # BED end is exclusive
            strand = info[strand_index]
            attributes = parse_gtf_attributes(info[8])  # Always column 9

            # Get necessary attributes
            transcript_id = attributes.get(transcript_id_key, None)
            gene_biotype = attributes.get(biotype_key, None)
            gene_name = attributes.get(gene_name_key, "NA")

            if transcript_id and transcript_id.startswith(id_head) and gene_biotype in biotypes and feature == feature_type:
                # Expand region 5kb downstream based on strand
                if strand == "+":
                    end += extension
                elif strand == "-":
                    start = max(0, start - extension)

                bed_line = f"{chro}\t{start}\t{end}\t{transcript_id}_{gene_name}\t{score}\t{strand}\n"
                bedout.write(bed_line)


if __name__ == "__main__":
    gtf_3utr_bed()

# Example usage:
#gtf_3utr_bed(
#    bedfile="data/regions_3utr.bed", 
#    gtffile="data/Homo_sapiens.GRCh38.99.gtf"
#)
