# minimiro



## Simple Usage
```
minimap2 -x asm20 -s {SCORE} -p 0.01 -N 1000 --cs {input.ref} {input.query} > {output.paf}

minimiro.py --paf {input.paf} --bestn 1000 -o {output.ps} && ps2pdf {output.ps}
```
## Anotation pipeline usage 

First make a file called `minimiro.yaml` that looks like this:
```
scores:  # minimum score threshold to plot, can specify multiple and it will make multiple pdfs
    - 50000
    - 150000


DEF_hg38_vs_CHM13: # run name for this comparison, will name output accordingly
    ref: /net/eichler/vol26/projects/sda_assemblies/nobackups/assemblies/hg38/ucsc.hg38.no_alts.fasta
    regions:
        - chr8:6000000-13500000
    genes: /net/eichler/vol26/projects/sda_assemblies/nobackups/assemblies/hg38/ucsc.refseq.genes.bed # bed12 gene file with ref coordiantes
    query: /net/eichler/vol27/projects/AlphaSatelliteMapping/nobackups/FindingAlphaSat/t2t_chr8/t2t_rel3_glchr8/v8_hybrid/asm/HiCanu/chm13_hicanu_hifi_20k_glchr8v8.fa
    queryregions:
        - glchr8v8:6000000-13500000
    rc: False # if set to True the query will be reverse complemened before displaying
```

