import subprocess
from pathlib import Path

# Path objektumok (fájlrendszer-útvonalak)
projectRoot = Path.home()/"projects"/"cytodna"

trimmedDir = projectRoot/"results"/"trimmed_fastq"
alignDir = projectRoot/"results"/"alignment"
logDir = projectRoot/"logs"/"bwa_samtools"
# BWA által indexelt referencia genom (FASTA)
reference = projectRoot/"reference_data"/"grch38_v112"/"Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"

# trimmelt R1 fájlok keresése
r1_files = sorted(trimmedDir.glob("*_R1_trimmed.fastq.gz")) # rendezett lista _R1... végű fájlokkal

for r1 in r1_files:
    r2 = r1.with_name(r1.name.replace("_R1_trimmed.fastq.gz", "_R2_trimmed.fastq.gz"))

    sample_id = r1.name.replace("_R1_trimmed.fastq.gz", "")

    log_file = logDir/f"{sample_id}.log"
    sam_file = alignDir / f"{sample_id}.sam"
    bam_tmp = alignDir / f"{sample_id}.unsorted.bam"
    bam_file = alignDir / f"{sample_id}.bam"

    with open(log_file, "w") as log:
        
        # alignment -> sam fájl
        with open(sam_file, "w") as sam_out:
            subprocess.run(
                ["bwa", "mem", "-t", "4", 
                "-R", f"@RG\\tID:{sample_id}\\tSM:{sample_id}\\tPL:ILLUMINA",
                reference, r1, r2],
                stdout=sam_out,
                stderr=log,
                check=True
            )
        # sam -> bam konvertálás
        subprocess.run(
            ["samtools", "view", "-b", "-o", bam_tmp, sam_file],
            stderr=log,
            check=True
        )
        # bam genom koordináta szerinti rendezése
        subprocess.run(
            ["samtools", "sort", "-o", bam_file, bam_tmp],
            stderr=log,
            check=True
        )
        # indexelés -> .bam.bai fájl
        subprocess.run(
            ["samtools", "index", bam_file],
            stderr=log,
            check=True
        )

    # létrejött sam és köztes fájl törlése
    if sam_file.exists():
        sam_file.unlink()

    if bam_tmp.exists():
        bam_tmp.unlink()

    print(f"{sample_id} alignment kesz.")

print("Alignment kesz az osszes mintara.")
