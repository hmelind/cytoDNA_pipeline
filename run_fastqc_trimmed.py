import subprocess
from pathlib import Path

# Path objektumok (fájlrendszer-útvonalak)
projectRoot = Path.home()/"projects"/"cytodna"

trimmedDir = projectRoot/"results"/"trimmed_fastq"
qcDir = projectRoot/"results"/"qc"/"trimmed"/"trimmed_fastqc"
logDir = projectRoot/"logs"/"fastqc"

# trimmelt R1 fájlok keresése
r1_files = sorted(trimmedDir.glob("*_R1_trimmed.fastq.gz")) # rendezett lista _R1... végű fájlokkal

for r1 in r1_files:
    r2 = r1.with_name(r1.name.replace("_R1_trimmed.fastq.gz", "_R2_trimmed.fastq.gz"))

    sample_id = r1.name.replace("_R1_trimmed.fastq.gz", "")

    log_file = logDir/f"{sample_id}_trimmed_fastqc.log"

    # fastqc parancs lista
    command = [
        "fastqc",
        "-t", "4", # 4 szálas CPU használat
        "-o", qcDir, # output
        r1, # input
        r2, # input
    ]
    with open(log_file, "w") as log:
        subprocess.run(command, stdout=log, stderr=log, check=True)

    print(f"{sample_id} FastQC kesz.")

print("FastQC lefutott az osszes trimmelt mintara.")

# MultiQC futtatása
multiqc_log = logDir/"multiqc_trimmed.log"

multiqc_command = [
    "multiqc",
    qcDir, # input a mappa tartalma
    "-o", qcDir # output
]

with open(multiqc_log, "w") as log:
    subprocess.run(multiqc_command, stdout=log, stderr=log, check=True)

print("MultiQC report kesz.")
