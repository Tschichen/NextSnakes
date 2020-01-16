DEDUPBIN, DEDUPENV = env_bin_from_config2(SAMPLES,config,'DEDUP')

#wildcard_constraints:
#    feat="!os.sep()"

rule all:
    input:  "DEDUP/DONE",
            expand("DEDUP/{file}.bam", file=samplecond(SAMPLES,config))

rule remove_dups:
    input:  "SORTED_MAPPED/{file}_mapped_{cond}_sorted.bam"
    output: report("SORTED_MAPPED/{file}_mapped_{cond}_sorted_dedup.bam", category="DEDUP"),
            "LOGS/{file}_{cond}_dedup.txt",
            "SORTED_MAPPED/{file}_mapped_{cond}_sorted_dedup.bam.bai"
    log:    "LOGS/{file}/remove_dups_{cond}.log"
    conda: "../envs/picardtools.yaml"
    threads: 1
    shell: "picard MarkDuplicates I={input[0]} O={output[0]} M={output[1]} ASO=coordinate REMOVE_SEQUENCING_DUPLICATES=true REMOVE_DUPLICATES=true && samtools index {output[0]} 2> {log}"
#for newer picard verions       "picard MarkDuplicates -I {input[0]} -O {output[0]} -M {output[1]} -ASO coordinate -REMOVE_SEQUENCING_DUPLICATES true -REMOVE_DUPLICATES true && samtools index {output[0]}"

rule remove_dups_unique:
    input:  "UNIQUE_MAPPED/{file}_mapped_{cond}_sorted_unique.bam",
            "QC/{file}_mapped_{cond}_sorted_unique_fastqc.zip"
    output: report("UNIQUE_MAPPED/{file}_mapped_{cond}_sorted_unique_dedup.bam", category="UNIQUE"),
            "LOGS/{file}_{cond}_unique_dedup.txt",
            "UNIQUE_MAPPED/{file}_mapped_{cond}_sorted_unique_dedup.bam.bai"
    log:    "LOGS/{file}/remove_dups_unique_{cond}.log"
    conda: "../envs/picardtools.yaml"
    threads: 1
    shell:  "picard MarkDuplicates I={input[0]} O={output[0]} M={output[1]} ASO=coordinate REMOVE_SEQUENCING_DUPLICATES=true REMOVE_DUPLICATES=true && samtools index {output[0]} 2> {log}"

rule featurecount_unique:
    input:  u = "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam",
    output: c = "COUNTS/Featurecounter_{feat}s/{file}_mapped_sorted_unique.counts",
            t = temp("COUNTS/Featurecounter_{feat}s/{file}_unique.anno")
    log:    "LOGS/{file}/featurecount_{feat}s_unique.log"
    conda:  "snakes/envs/"+COUNTENV+".yaml"
    threads: MAXTHREAD
    params: count = COUNTBIN,
            anno = lambda wildcards: str.join(os.sep,[config["REFERENCE"],os.path.dirname(genomepath(wildcards.file, config)),tool_params(wildcards.file, None, config, 'COUNTING')['ANNOTATION']]),
            cpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "COUNTING")['OPTIONS'][0].items())+' -t '+wildcards.feat+' -g '+config['COUNTING']['FEATURES'][wildcards.feat],
            paired = lambda x: '-p' if paired == 'paired' else '',
            stranded = lambda x: '-s 1' if stranded == 'fr' else '-s 2' if stranded == 'rf' else ''
    shell:  "zcat {params.anno} > {output.t} && {params.count} -T {threads} {params.cpara} {params.paired} {params.stranded} -a {output.t} -o {output.c} {input.u} 2> {log}"

rule count_dedup_mappers:
    input:  u = "DEDUP/{file}_mapped_sorted_dedup.bam"
    output: u = "COUNTS/{file}_mapped_unique.count"
    log:    "LOGS/{file}/count_unique_mappers.log"
    conda:  "snakes/envs/samtools.yaml"
    threads: MAXTHREAD
    shell:  "export LC_ALL=C; arr=({input.u}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do samtools view -F 260 ${{arr[$i]}} | cut -d$'\t' -f1|sort --parallel={threads} -S 25% -T SORTTMP -u |wc -l > {output.u} ;done 2>> {log}"
