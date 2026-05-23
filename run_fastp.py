import subprocess # külső program futtatásához
import re # reguláris kifejezésekhez - névmódosítás
from pathlib import Path

# Path objektumok (fájlrendszer-útvonalak)
projectRoot = Path.home()/"projects"/"cytodna"

rawDir = projectRoot/"raw_data"
trimmedDir = projectRoot/"results"/"trimmed_fastq"
qcDir = projectRoot/"results"/"qc"/"trimmed"
logDir = projectRoot/"logs"/"fastp"

# R1 fájlok keresése
r1_files = sorted(rawDir.glob("*_R1_001.fastq.gz")) # rendezett lista _R1... végű fájlokkal
"""print(f"{len(r1_files)} db R1 fajl van")"""

for r1 in r1_files:
    r2 = r1.with_name(r1.name.replace("_R1_001.fastq.gz","_R2_001.fastq.gz")) # R1 elérési útvonala marad, hozzárendeljük a megfelelő R2-t
    
    # minta alap neve
    sample_id = r1.name.replace("_R1_001.fastq.gz", "") # fastq suffix eltávolítás
    sample_id = re.sub(r"_S\d+_L001$", "", sample_id) # Illumina run meta eltávolítás
    sample_id = re.sub(r"^RPE-1-", "", sample_id) # sejtvonal prefix eltávolítás
    """print(f"Aktualis minta: {sample_id}")"""

    # output fájlok helye, elnevezése
    outp_r1 = trimmedDir/f"{sample_id}_R1_trimmed.fastq.gz"
    outp_r2 = trimmedDir/f"{sample_id}_R2_trimmed.fastq.gz"

    html_report = qcDir/f"{sample_id}_fastp.html"
    json_report = qcDir/f"{sample_id}_fastp.json"

    log_file = logDir/f"{sample_id}.log"

    # fastp parancs lista
    command = [
        "fastp",
        "-i", r1,
        "-I", r2,
        "-o", outp_r1,
        "-O", outp_r2,
        "--thread", "4", # 4 szálas CPU használat
        "--qualified_quality_phred", "28", # egy bázis akkor „qualified”, ha a Phred értéke eléri a 28-at
        "--trim_front1", "3", # R1 5' végen 3 bázis levágása
        "--trim_front2", "3", # R2 5' végen 3 bázis levágása
        "--trim_tail1", "3", # R1 3' végen 3 bázis levágása
        "--trim_tail2", "3", # R2 3' végen 3 bázis levágása
        "--cut_right", # readek rossz minőségű végének eltávolítása
        "--cut_window_size", "4", # a végek levágása 4-es mozgó ablakban történik
        "--cut_mean_quality", "15", # a végek mozgó átlagának legalacsonyabb elfogadható értéke 15
        "--length_required", "38", # minimum read hosszúság a trimming után 38bp
        "-x", # poly X farok vágása 
        "-g", # poly G farok vágása, elsőként ez történik, ha más vágás is aktív
        "--detect_adapter_for_pe", # automatikus adapter vágás
        "--html", html_report,
        "--json", json_report
    ]
    
    # fastp futtatása
    with open(log_file, "w") as log:
        subprocess.run(command, stdout=log, stderr=log, check=True)

    print(f"{sample_id} trimming kesz.")

print("A Fastp trimming lefuttott az osszes mintara.")
