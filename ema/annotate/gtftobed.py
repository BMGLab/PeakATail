def gtf_bed(endbeddir:str,
             biotypes:tuple, 
             source_tuple:tuple,
             id_head:tuple,
             gtfdir:str,
             featuresdir:str,
             chro_index=1,
             source_index=2,
             start_index=3,
             end_index=4,
             strand_index=6,
             id_index=9,
             biotype_index=17,
             symbol_index=19,
             score=0
             ) -> None:

    with open(gtfdir, "r") as gtf, open(endbeddir, "w") as endbed, open(featuresdir, "w") as features:
        #ignore tuple
        ignore = "#"
        #index of id and biotype in gtf file format

        for line in gtf:

            if line.startswith(ignore):
                continue

            info = line.split()
            # get necessary data from gtf line
            # ex : suorce = "gene", id = "ENSG", gene_biotype = "protein_coding" 
            source, id, gene_biotype = info[source_index], info[id_index].replace('"','').replace(";",'')[:4], info[biotype_index].replace('"','').replace(";",'')

            if source in source_tuple and id in id_head and gene_biotype in biotypes:
                chro, start, end, strand, symbol = info[chro_index], int(info[start_index]), int(info[end_index]), info[strand_index], info[symbol_index]

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