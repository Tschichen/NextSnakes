rule trimgalore_trim:
    input:  r1 = lambda wildcards: "FASTQ/{rawfile}.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0])
    output: o1 = "TRIMMED_FASTQ/{file}_trimmed.fq.gz"
    log:    "LOGS/{file}_trim.log"
    conda: "../envs/"+TRIMENV+".yaml"
    threads: int(MAXTHREAD/2)
    params: odir=lambda wildcards,output: os.path.dirname(output.o1),
            tpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "TRIMMING")[0].items()),
            trim=TRIMBIN,
            cores = min(int(MAXTHREAD/4),4)
    shell:  "{params.trim} --cores {params.cores} --paired --no_report_file --gzip {params.tpara} -o {params.odir} {input.r1} &> {log}"

rule trimgalore_rename:
    input:  o1 = rules.trimgalore_trim.output.o1
    output: r1 = "TRIMMED_FASTQ/{file}_trimmed.fastq.gz"
    conda: "../envs/"+TRIMENV+".yaml"
    threads: 1
    shell:  "mv {input.o1} {output.r1}"
