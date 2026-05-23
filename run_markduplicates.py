import subprocess
from pathlib import Path

projectRoot = Path.home()/"projects"/"cytodna"

alignDir = projectRoot/"results"/"alignment"
dedupDir = projectRoot/"results"/"dedup"
logDir = projectRoot/"logs"/"picard"

bam_files = sorted(alignDir.glob("*.bam")) # .bam fájlok keresése 
print(f"{len(bam_files)} BAM fajl talalva.")

for bam in bam_files:

    sample_id = bam.stem # fájlnév kiterjesztés nélkül
    
    output_bam = dedupDir/f"{sample_id}_marked_dup.bam"
    metrics_file = dedupDir/f"{sample_id}_marked_dup_metrics.txt"
    log_file = logDir/f"{sample_id}.log"

    with open(log_file, "w") as log:
        # MarkDuplicates futtatás
        subprocess.run(
            [
                "picard", "MarkDuplicates",
                f"I={bam}",
                f"O={output_bam}",
                f"M={metrics_file}",
                "REMOVE_DUPLICATES=true", # a jelölt duplikátumok eltávolítása
                "CREATE_INDEX=true" # új .bam.bai indexelt fájl létrehozása
            ],
            stdout=log,
            stderr=log,
            check=True
        )

    print(f"{sample_id} MarkDuplicates kesz.")

print("MarkDuplicates lefutott az osszes mintara.")
