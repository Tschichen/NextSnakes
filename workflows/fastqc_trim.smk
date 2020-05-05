if paired == 'paired':
    log.info('Running paired mode QC')
    rule qc_raw:
        input:  r1 = expand("FASTQ/{rawfile}_{read}.fastq.gz", rawfile=list(SAMPLES), read=['R1','R2'])
        output: o1 = report(expand("QC/{rawfile}_{read}_fastqc.zip", rawfile=list(SAMPLES), read=['R1','R2']), category="QC")
        log:    expand("LOGS/{rawfile}/fastqc_{read}_raw.log", rawfile=list(SAMPLES), read=['R1','R2'])
        conda:  "snakes/envs/qc.yaml"
        threads: MAXTHREAD
        params:  qpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(SAMPLES[0], None ,config, 'QC')['OPTIONS'][0].items())
        shell: "arr=({input.r1}); orr=({output.o1});alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do OUT=$(dirname ${{orr[$i]}}); fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f fastq  ${{arr[$i]}} 2> {log};done"

    rule qc_trimmed:
        input:  r1 = expand("TRIMMED_FASTQ/{file}_{read}_trimmed.fastq.gz", file=samplecond(SAMPLES,config), read=["R1", "R2"])
        output: o1 = report(expand("QC/{file}_{read}_trimmed_fastqc.zip", file=samplecond(SAMPLES,config), read=['R1','R2']), category="QC")
        log:    expand("LOGS/{file}/fastqc_{read}_trimmed.log", file=samplecond(SAMPLES,config), read=['R1','R2'])
        conda:  "snakes/envs/qc.yaml"
        threads: MAXTHREAD
        params:  qpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(SAMPLES[0], None ,config, 'QC')['OPTIONS'][0].items())
        shell: "arr=({input.r1}); orr=({output.o1});alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do OUT=$(dirname ${{orr[$i]}}); fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f fastq  ${{arr[$i]}} 2> {log};done"

    rule multiqc:
        input:  expand(rules.qc_raw.output.o1, rawfile=list(SAMPLES), read=['R1','R2']),
                expand(rules.qc_trimmed.output.o1, file=samplecond(SAMPLES,config), read=['R1','R2'])
        output: html = report("QC/Multi/TRIMMED_RAW/{condition}/multiqc_report.html", category="QC"),
                tmp = temp("QC/Multi/TRIMMED_RAW/{condition}/tmp"),
                lst = "QC/Multi/TRIMMED_RAW/{condition}/qclist.txt"
        log:    "LOGS/{condition}/multiqc_trimmed_raw.log"
        conda:  "snakes/envs/qc.yaml"
        threads: 1
        shell:  "OUT=$(dirname {output.html}); for i in {input};do echo $(dirname \"${{i}}\") >> {output.tmp};done; cat {output.tmp} |sort -u > {output.lst};export LC_ALL=en_US.utf8; export LC_ALL=C.UTF-8; multiqc -f --exclude picard --exclude gatk -k json -z -o $OUT -l {output.lst} 2> {log}"

else:
    rule qc_raw:
        input:  r1 = expand("FASTQ/{rawfile}.fastq.gz", rawfile=list(SAMPLES))
        output: o1 = report(expand("QC/{rawfile}_fastqc.zip", rawfile=list(SAMPLES)), category="QC")
        log:    expand("LOGS/{rawfile}/fastqc_raw.log", rawfile=list(SAMPLES))
        conda:  "snakes/envs/qc.yaml"
        threads: MAXTHREAD
        params:  qpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(SAMPLES[0], None ,config, 'QC')['OPTIONS'][0].items())
        shell: "arr=({input.r1}); orr=({output.o1});alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do OUT=$(dirname ${{orr[$i]}}); fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f fastq  ${{arr[$i]}} 2>> {log};done"

    rule qc_trimmed:
        input:  r1 = expand("TRIMMED_FASTQ/{file}_trimmed.fastq.gz", file=samplecond(SAMPLES,config))
        output: o1 = report(expand("QC/{file}_trimmed_fastqc.zip", file=samplecond(SAMPLES,config)), category="QC")
        log:    expand("LOGS/{file}/fastqc_trimmed.log", file=samplecond(SAMPLES,config))
        conda:  "snakes/envs/qc.yaml"
        threads: MAXTHREAD
        params:  qpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(SAMPLES[0], None ,config, 'QC')['OPTIONS'][0].items())
        shell: "arr=({input.r1}); orr=({output.o1});alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do OUT=$(dirname ${{orr[$i]}}); fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f fastq  ${{arr[$i]}} 2>> {log};done"

    rule multiqc:
        input: expand(rules.qc_raw.output.o1, rawfile=list(SAMPLES)),
               expand(rules.qc_trimmed.output.o1, file=samplecond(SAMPLES,config))
        output: html = report("QC/Multi/TRIMMED_RAW/{condition}/multiqc_report.html", category="QC"),
                tmp = temp("QC/Multi/TRIMMED_RAW/{condition}/tmp"),
                lst = "QC/Multi/TRIMMED_RAW/{condition}/qclist.txt"
        log:    "LOGS/{condition}/multiqc_trimmed_raw.log"
        conda:  "snakes/envs/qc.yaml"
        threads: 1
        shell:  "OUT=$(dirname {output.html}); for i in {input};do echo $(dirname \"${{i}}\") >> {output.tmp};done; cat {output.tmp} |sort -u > {output.lst};export LC_ALL=en_US.utf8; export LC_ALL=C.UTF-8; multiqc -f --exclude picard --exclude gatk -k json -z -o $OUT -l {output.lst} 2> {log}"
