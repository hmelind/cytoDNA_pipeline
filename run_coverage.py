import subprocess # külső programok futtatásához
from pathlib import Path # útvonalak kezelésére

# Path objektumok, mappák
projectRoot = Path.home() / "projects" / "cytodna"
referenceDir = projectRoot / "reference_data"
dedupDir = projectRoot / "results" / "dedup"
cleanDir = projectRoot / "results" / "blacklist_clean"
coverageDir = projectRoot / "results" / "coverage"

# konkrét fájlok
bins = referenceDir / "grch38_100kbp_bins.sorted.bed" # teljes genom 100k bp-os rácsa; makewindows output
blacklist = referenceDir / "hg38-blacklist.v2.nochr.bed" # kromoszóma nevezéktan módosított ENCODE blacklist

bam_files = sorted(dedupDir.glob("*_marked_dup.bam")) # bam fájlok megtalálása, rendezése
print(f"{len(bam_files)} BAM fajl van coverage szamitasra.")

for bam in bam_files[2:]:

    sample_id = bam.stem.replace("_marked_dup", "")

    clean_unsorted = cleanDir / f"{sample_id}_clean.unsorted.bam"
    clean_sorted = cleanDir / f"{sample_id}_clean.bam"
    coverage_output = coverageDir / f"{sample_id}_100kbp_coverage.txt"
    
    # blacklist szűrés read szinten
    with open(clean_unsorted, "wb") as out:
        subprocess.run(
            [
                "bedtools", "intersect",
                "-v",
                "-abam", str(bam),
                "-b", str(blacklist)
            ],
            stdout=out,
            check=True
        )
    print(f"{sample_id} blacklist torles kesz.")
    
    # tisztított bam fájl rendezése
    subprocess.run(
        [
            "samtools", "sort",
            "-o", str(clean_sorted),
            str(clean_unsorted)
        ],
        check=True
    )       
    # ideiglenes unsorted bam törlése
    clean_unsorted.unlink()

    # coverage
    with open(coverage_output, "w") as out:
        subprocess.run(
            [
                "bedtools", "coverage",
                "-a", str(bins), # bed fájl
                "-b", str(clean_sorted), # bam fájl
                "-sorted"
            ],
            stdout=out,
            check=True
        )

    print(f"{sample_id} coverage kesz.")

print("Blacklist torles és coverage szamitas kesz minden mintara.")
