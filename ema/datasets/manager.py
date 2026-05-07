from pathlib import Path
import pysam
import tempfile


class DatasetManager:
    def __init__(self, output_dir: Path):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def prepare(self, datasets: list[dict]) -> list[tuple[str, Path]]:
        """
        Given a list of dataset dicts from YAML config, returns list of (dataset_id, bam_path).
        For merge_strategy=before: merges all BAMs into one tagged BAM, returns single entry.
        For merge_strategy=after or none: returns one entry per BAM file.
        """
        result = []
        for ds in datasets:
            dataset_id = ds['id']
            bams = [Path(b) for b in ds['bams']]
            strategy = ds.get('merge_strategy', 'none')

            if strategy == 'before':
                merged = self._merge_with_rg(dataset_id, bams)
                result.append((dataset_id, merged))
            else:
                for bam in bams:
                    result.append((dataset_id, bam))
        return result

    def _merge_with_rg(self, dataset_id: str, bams: list[Path]) -> Path:
        """
        Tags all reads in all BAMs with RG=dataset_id, then merges and sorts.
        Returns path to merged sorted BAM in output_dir.
        """
        tagged_bams = []
        for bam_path in bams:
            tagged = self._tag_bam_with_rg(bam_path, dataset_id)
            tagged_bams.append(str(tagged))

        merged_path = self.output_dir / f"{dataset_id}_merged.bam"
        sorted_path = self.output_dir / f"{dataset_id}_merged_sorted.bam"

        if len(tagged_bams) == 1:
            pysam.sort("-o", str(sorted_path), tagged_bams[0])
        else:
            pysam.merge("-f", str(merged_path), *tagged_bams)
            pysam.sort("-o", str(sorted_path), str(merged_path))
            merged_path.unlink(missing_ok=True)

        pysam.index(str(sorted_path))

        # clean up tagged temp bams
        for t in tagged_bams:
            Path(t).unlink(missing_ok=True)

        return sorted_path

    def _tag_bam_with_rg(self, bam_path: Path, dataset_id: str) -> Path:
        """
        Writes a new BAM with RG header entry and RG tag set on every read.
        """
        out_path = self.output_dir / f"tagged_{dataset_id}_{bam_path.stem}.bam"
        with pysam.AlignmentFile(str(bam_path), "rb") as bam_in:
            header = bam_in.header.to_dict()
            rg_entry = {"ID": dataset_id, "SM": dataset_id}
            header.setdefault("RG", [])
            if not any(rg["ID"] == dataset_id for rg in header["RG"]):
                header["RG"].append(rg_entry)

            with pysam.AlignmentFile(str(out_path), "wb", header=header) as bam_out:
                for read in bam_in:
                    read.set_tag("RG", dataset_id)
                    bam_out.write(read)

        return out_path
