from pathlib import Path
import subprocess
import pysam
import tempfile


def _check_samtools_version(min_version: tuple[int, int] = (1, 10)) -> None:
    """Raise RuntimeError if samtools is missing or < min_version."""
    try:
        result = subprocess.run(
            ["samtools", "--version"], capture_output=True, check=True
        )
    except (FileNotFoundError, subprocess.CalledProcessError) as e:
        raise RuntimeError("samtools not found in PATH") from e
    # Parse first line: "samtools 1.13" or "samtools 1.10.1"
    # Decode with errors='replace' since samtools' license blurb may include non-UTF8 chars
    first_line = result.stdout.decode("utf-8", errors="replace").split("\n")[0]
    parts = first_line.split()
    if len(parts) < 2:
        raise RuntimeError(f"Cannot parse samtools version from: {first_line}")
    version_str = parts[1]
    version_parts = tuple(int(x) for x in version_str.split(".")[:2])
    if version_parts < min_version:
        raise RuntimeError(
            f"samtools >= {'.'.join(map(str, min_version))} required, got {version_str}"
        )


class DatasetManager:
    def __init__(self, output_dir: Path, threads: int = 4):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.threads = threads
        _check_samtools_version()

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
        """Tag all reads with RG=dataset_id via `samtools addreplacerg`.

        Falls back to a clear error if samtools is missing or version < 1.10.

        Args:
            bam_path: Path to the input BAM file.
            dataset_id: Read-group ID (and SM) to stamp on every read.

        Returns:
            Path to the newly written, tagged BAM file.

        Raises:
            RuntimeError: If samtools is not on PATH or the command exits non-zero.
        """
        out_path = self.output_dir / f"tagged_{dataset_id}_{bam_path.stem}.bam"
        threads = getattr(self, 'threads', 4)
        cmd = [
            "samtools", "addreplacerg",
            "-r", f"ID:{dataset_id}\tSM:{dataset_id}",
            "-m", "overwrite_all",   # ensures any existing RG is replaced, not appended
            "-@", str(threads),
            "-o", str(out_path),
            str(bam_path),
        ]
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
        except FileNotFoundError as e:
            raise RuntimeError(
                "samtools not found in PATH. Install samtools >= 1.10."
            ) from e
        except subprocess.CalledProcessError as e:
            raise RuntimeError(
                f"samtools addreplacerg failed: {e.stderr}"
            ) from e
        return out_path
