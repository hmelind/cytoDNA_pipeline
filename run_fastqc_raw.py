import subprocess
from pathlib import Path

projectRoot = Path.home()/"projects"/"cytodna"

rawDir = projectRoot/"raw_data"
qcDir = projectRoot/"results"/"qc"/"raw"
logDir = projectRoot/"logs"/"fastqc"

# R1 fájlok keresése
r1_files = sorted(rawDir.glob("*_R1_001.fastq.gz"))

for r1 in r1_files:
    r2 = r1.with_name(r1.name.replace("_R1_001.fastq.gz", "_R2_001.fastq.gz"))

    sample_id = r1.name.replace("_R1_001.fastq.gz", "")

    log_file = logDir/f"{sample_id}_raw_fastqc.log"

    command = [
        "fastqc",
        "-t", "4",
        "-o", qcDir,
        r1,
        r2,
    ]

    with open(log_file, "w") as log:
        subprocess.run(command, stdout=log, stderr=log, check=True)

    print(f"{sample_id} FastQC kesz.")

print("FastQC lefutott az osszes mintara.")

# MultiQC
multiqc_log = logDir/"multiqc_raw.log"

multiqc_command = [
    "multiqc",
    qcDir,
    "-o", qcDir
]

with open(multiqc_log, "w") as log:
    subprocess.run(multiqc_command, stdout=log, stderr=log, check=True)

print("MultiQC kesz.")