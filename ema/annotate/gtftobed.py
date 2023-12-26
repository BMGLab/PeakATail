from ema.config import directory_config as dc
def gtf_bed(endbeddir= dc.endbed,
             gtfdir=dc.gtf_dir,
             featuresdir=dc.raw_features,
             biotypes=("protein_coding", "lincRNA", "snRNA", "antisense",  "miRNA", "processed_transcript"), 
             source_tuple=("gene"),
             id_head=("ENSG"),
             chro_index=0,
             source_index=2,
             start_index=3,
             end_index=4,
             strand_index=6,
             id_index=9,
             biotype_index=17,
             symbol_index=13,
             score=0
             ) -> None:

    with open(gtfdir, "r") as gtf, open(endbeddir, "w") as endbed, open(featuresdir, "w") as features:
        #ignore tuple
        ignore = "#"
        

        for line in gtf:

            if line.startswith(ignore):
                continue

            info = line.split()
            # get necessary data from gtf line
            # ex : suorce = "gene", id = "ENSG", gene_biotype = "protein_coding" 
            source, id, gene_biotype = info[source_index], info[id_index].replace('"','').replace(";",''), info[biotype_index].replace('"','').replace(";",'')

            if source in source_tuple and id[:4] in id_head and gene_biotype in biotypes:
                chro, start, end, strand, symbol = info[chro_index], int(info[start_index]), int(info[end_index]), info[strand_index], info[symbol_index].replace('"','').replace(";",'')

                match strand:
                    case "+":
                        start = end - 1
                    
                    case "-":
                        end = start + 1

                bed_info = f"{chro}\t{start}\t{end}\t{id}\t{score}\t{strand}\n"
                endbed.write(bed_info)\
                
                feature = f"{id}\t{symbol}\n"
                features.write(feature)

if __name__ == "__main__":
    gtf_bed()