# -------- hifiasm rules -------- #
rule run_hifiasm:
    output:
        hp1_gfa="results_hs/hs-{k}/{sample_id}/hifiasm/{sample_id}.bp.hap1.p_ctg.gfa",
        hp2_gfa="results_hs/hs-{k}/{sample_id}/hifiasm/{sample_id}.bp.hap2.p_ctg.gfa",
        ec_fq="results_hs/hs-{k}/{sample_id}/hifiasm/{sample_id}.ec.fq"
    input:
        ont_r10_sequence=config["raw_reads_dir"] + "/{sample_id}/{sample_id}.fastq"
    benchmark: "benchmarks/{sample_id}/hs-{k}/run_hifiasm.benchmark.txt"
    log: "logs/{sample_id}/hs-{k}/run_hifiasm.log"
    threads: 128
    container: "docker://quay.io/shnegi/hifiasm:0.25.0"
    shell:
        """
        hifiasm --ont -t {threads} -o results_hs/hs-{wildcards.k}/{wildcards.sample_id}/hifiasm/{wildcards.sample_id} --write-ec {input.ont_r10_sequence} > {log} 2>&1
        """

rule convert_ec_fq_to_fasta:
	output:
		fasta="results_hs/hs-{k}/reads/{sample_id}.ec.fasta"
	input:
		ec_fq="results_hs/hs-{k}/{sample_id}/hifiasm/{sample_id}.ec.fq"
	shell:
		"""
		# Convert FASTQ to FASTA
		cat {input.ec_fq} | awk '{{if(NR%4==1) {{printf(">%s\\n",substr($0,2));}} else if(NR%4==2) print;}}' > {output.fasta}
		"""

rule convert_gfa_to_fasta:
    output:
        hp1_fasta="results_hs/hs-{k}/{sample_id}/hifiasm/{sample_id}.bp.hap1.p_ctg.fasta",
        hp2_fasta="results_hs/hs-{k}/{sample_id}/hifiasm/{sample_id}.bp.hap2.p_ctg.fasta"
    input:
        hp1_gfa="results_hs/hs-{k}/{sample_id}/hifiasm/{sample_id}.bp.hap1.p_ctg.gfa",
        hp2_gfa="results_hs/hs-{k}/{sample_id}/hifiasm/{sample_id}.bp.hap2.p_ctg.gfa"
    benchmark: "benchmarks/{sample_id}/hs-{k}/convert_gfa_to_fasta.benchmark.txt"
    log: "logs/{sample_id}/hs-{k}/convert_gfa_to_fasta.log"
    threads: 128
    container: "docker://quay.io/shnegi/hifiasm:0.25.0"
    shell:
        """
        gfatools gfa2fa {input.hp1_gfa} > {output.hp1_fasta}
        gfatools gfa2fa {input.hp2_gfa} > {output.hp2_fasta}
        """

rule combine_hifiasm_hap1_and_hap2_fastas:
    output:
        hifiasm_assembly="results_hs/hs-{k}/{sample_id}/hifiasm_alignment/{sample_id}.hifiasm.fasta"
    input:
        hap1_fasta="results_hs/hs-{k}/{sample_id}/hifiasm/{sample_id}.bp.hap1.p_ctg.fasta",
        hap2_fasta="results_hs/hs-{k}/{sample_id}/hifiasm/{sample_id}.bp.hap2.p_ctg.fasta"
    shell:
        """
        cat {input.hap1_fasta} {input.hap2_fasta} > {output.hifiasm_assembly}
        """

rule align_hifiasm_assembly_to_reference:
    output:
        bam="results_hs/hs-{k}/{sample_id}/hifiasm_alignment/{sample_id}_hifiasm_to_hg002_minimap_{asm_preset}.bam",
        bai="results_hs/hs-{k}/{sample_id}/hifiasm_alignment/{sample_id}_hifiasm_to_hg002_minimap_{asm_preset}.bam.bai"
    input:
        hg002_reference=config["HG002v101_ref"],
        hifiasm_assembly="results_hs/hs-{k}/{sample_id}/hifiasm_alignment/{sample_id}.hifiasm.fasta"
    params:
        asm_preset=config["MINIMAP"]["asmPreset"]
    benchmark: "benchmarks/{sample_id}/hs-{k}/align_hifiasm_assembly_to_reference_{asm_preset}.benchmark.txt"
    log: "logs/{sample_id}/hs-{k}/align_hifiasm_assembly_to_reference_{asm_preset}.log"
    threads: 64
    container: "docker://mkolmogo/card_minimap2:2.23"
    shell:
        """
        minimap2 -t {threads} -I 20G -ax {params.asm_preset} -K 1M --eqx --cs {input.hg002_reference} {input.hifiasm_assembly} | \
            samtools sort -@ {threads} -o {output.bam} -
        samtools index {output.bam}
        """

#-------- hifiasm r_utg assembly rules --------#
rule convert_r_utg_gfa_to_fasta:
    output:
        r_utg_fasta="results_hs/hs-{k}/{sample_id}/hifiasm/{sample_id}.bp.r_utg.fasta"
    benchmark: "benchmarks/{sample_id}/hs-{k}/convert_r_utg_gfa_to_fasta.benchmark.txt"
    log: "logs/{sample_id}/hs-{k}/convert_r_utg_gfa_to_fasta.log"
    container: "docker://quay.io/shnegi/hifiasm:0.25.0"
    shell:
        """
        R_UTG_GFA=results_hs/hs-{wildcards.k}/{wildcards.sample_id}/hifiasm/{wildcards.sample_id}.bp.r_utg.gfa
        gfatools gfa2fa $R_UTG_GFA > {output.r_utg_fasta}
        """

rule align_r_utg_hifiasm_assembly_to_reference:
    output:
        bam="results_hs/hs-{k}/{sample_id}/hifiasm_alignment/{sample_id}_r_utg_to_hg002_minimap_{asm_preset}.bam",
        bai="results_hs/hs-{k}/{sample_id}/hifiasm_alignment/{sample_id}_r_utg_to_hg002_minimap_{asm_preset}.bam.bai"
    input:
        hg002_reference=config["HG002v101_ref"],
        r_utg_fasta="results_hs/hs-{k}/{sample_id}/hifiasm/{sample_id}.bp.r_utg.fasta"
    params:
        asm_preset=config["MINIMAP"]["asmPreset"]
    benchmark: "benchmarks/{sample_id}/hs-{k}/align_r_utg_hifiasm_assembly_to_reference_{asm_preset}.benchmark.txt"
    log: "logs/{sample_id}/hs-{k}/align_r_utg_hifiasm_assembly_to_reference_{asm_preset}.log"
    threads: 64
    container: "docker://mkolmogo/card_minimap2:2.23"
    shell:
        """
        minimap2 -t {threads} -I 20G -ax {params.asm_preset} -K 1M --eqx --cs {input.hg002_reference} {input.r_utg_fasta} | \
            samtools sort -@ {threads} -o {output.bam} -
        samtools index {output.bam}
        """
