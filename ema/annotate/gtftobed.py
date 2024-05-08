from ema.config import directory_config as dc
def gtf_bed(endbeddir= dc.endbed,
             gtfdir=dc.gtf_dir,
             featuresdir=dc.raw_features,
             biotypes=("protein_coding", "lncRNA", "snRNA", "antisense",  "miRNA", "processed_transcript", "lincRNA"), 
             source_tuple=("three_prime_utr"),
             id_head=("ENSG"),
             transcript_id_head = ("ENST"),
             chro_index=0,
             source_index=2,
             start_index=3,
             end_index=4,
             strand_index=6,
             id_index=9,
             transcript_id = 13,
             biotype_index=21,
             symbol_index=17,
             score=0
             ) -> None:
    """
    Convert a GTF file to BED format and extract features.

    Args:
        endbeddir (str): Directory path to write the converted BED file.
        gtfdir (str): Directory path of the input GTF file.
        featuresdir (str): Directory path to write the extracted features.
        biotypes (tuple): Tuple of biotypes to include in the conversion.
        source_tuple (tuple): Tuple of source types to include in the conversion.
        id_head (tuple): Tuple of ID headers to include in the conversion.
        chro_index (int): Index of the chromosome column in the GTF file.
        source_index (int): Index of the source column in the GTF file.
        start_index (int): Index of the start position column in the GTF file.
        end_index (int): Index of the end position column in the GTF file.
        strand_index (int): Index of the strand column in the GTF file.
        id_index (int): Index of the ID column in the GTF file.
        biotype_index (int): Index of the biotype column in the GTF file.
        symbol_index (int): Index of the symbol column in the GTF file.
        score (int): Score value for the BED file.

    Returns:
        None: This function does not return anything.
    """
    with open(gtfdir, "r") as gtf, open(endbeddir, "w") as endbed, open(featuresdir, "w") as features:
        #ignore tuple
        ignore = "#"
        
        for line in gtf:
            if line.startswith(ignore):
                continue

            info = line.split()
            # get necessary data from gtf line
            # ex : suorce = "gene", id = "ENSG", gene_biotype = "protein_coding" 
            source, id,t_script_id, gene_biotype = info[source_index], info[id_index].replace('"','').replace(";",''),info[id_index].replace('"','').replace(";",''), info[biotype_index].replace('"','').replace(";",'')

            if source in source_tuple and t_script_id[:4] in transcript_id_head and gene_biotype in biotypes:
                chro, start, end, strand, symbol = info[chro_index], int(info[start_index]), int(info[end_index]), info[strand_index], info[symbol_index].replace('"','').replace(";",'')

                match strand:
                    case "+":
                        end += 5000
                    
                    case "-":
                        if start > 5000:
                            start -= 5000
                        else:
                            start = 1

                bed_info = f"{chro}\t{start}\t{end}\t{t_script_id}\t{symbol}\t{strand}\n"
                endbed.write(bed_info)
                
                feature = f"{t_script_id}\t{symbol}\n"
                features.write(feature)

if __name__ == "__main__":
    gtf_bed()