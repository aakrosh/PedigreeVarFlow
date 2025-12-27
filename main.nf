// Description: Run a workflow that aligns a CRAM file using BWA-MEM.
nextflow.enable.dsl = 2
include { bgzip_index_variants as bgzip_index_raw;
          bgzip_index_variants as bgzip_index_standard;
          bgzip_index_variants as bgzip_index_vep } from "./modules/index"

process validate_inputs {
    tag "Input Validation"
    executor 'local'

    input:
    path samplesheet
    path reference
    path reference_in
    path pedigree
    path phenotype
    path target_regions
    path exclusion_bed

    output:
    val true

    script:
    """
    #!/usr/bin/env python3
    import sys
    import csv
    from pathlib import Path

    def error_exit(message):
        print(f"ERROR: {message}", file=sys.stderr)
        sys.exit(1)

    def warning(message):
        print(f"WARNING: {message}", file=sys.stderr)

    print("=" * 70)
    print("PedigreeVarFlow Input Validation")
    print("=" * 70)

    # Check reference genome
    print("\\n[1/7] Validating reference genome...")
    ref = Path("${reference}")
    if not ref.exists():
        error_exit(f"Reference genome not found: {ref}")

    # Check for reference index (.fai)
    fai = Path(f"{ref}.fai")
    if not fai.exists():
        error_exit(f"Reference index (.fai) not found: {fai}\\n" +
                   f"       Create with: samtools faidx {ref}")

    # Check for BWA index
    bwa_index = Path(f"{ref}.bwt")
    if not bwa_index.exists():
        error_exit(f"BWA index not found: {bwa_index}\\n" +
                   f"       Create with: bwa index {ref}")

    print(f"   [OK] Reference: {ref.name}")
    print(f"   [OK] Reference index found")
    print(f"   [OK] BWA index found")

    # Check input reference (for CRAM)
    print("\\n[2/7] Validating input reference (for CRAM)...")
    ref_in = Path("${reference_in}")
    if not ref_in.exists():
        error_exit(f"Input reference not found: {ref_in}")
    ref_in_fai = Path(f"{ref_in}.fai")
    if not ref_in_fai.exists():
        error_exit(f"Input reference index (.fai) not found: {ref_in_fai}")
    print(f"   [OK] Input reference: {ref_in.name}")

    # Check samplesheet
    print("\\n[3/7] Validating samplesheet...")
    samplesheet = Path("${samplesheet}")
    if not samplesheet.exists():
        error_exit(f"Samplesheet not found: {samplesheet}")

    samples = []
    with open(samplesheet) as f:
        reader = csv.DictReader(f)

        # Check required columns
        if 'sample' not in reader.fieldnames:
            error_exit("Samplesheet missing required column: 'sample'")
        if 'cram_location' not in reader.fieldnames:
            error_exit("Samplesheet missing required column: 'cram_location'")

        # Validate each row
        for i, row in enumerate(reader, start=2):
            sample = row.get('sample', '').strip()
            cram_path = row.get('cram_location', '').strip()

            if not sample:
                error_exit(f"Empty sample name at line {i}")

            if sample in samples:
                error_exit(f"Duplicate sample ID '{sample}' at line {i}")
            samples.append(sample)

            if not cram_path:
                error_exit(f"Empty CRAM path for sample '{sample}' at line {i}")

            # Check if CRAM file exists (only if not a URL)
            if not cram_path.startswith(('http://', 'https://', 's3://', 'gs://')):
                cram = Path(cram_path)
                if not cram.exists():
                    error_exit(f"CRAM file not found for sample '{sample}': {cram_path}")

                # Check for CRAM index
                crai = Path(f"{cram_path}.crai")
                if not crai.exists():
                    warning(f"CRAM index not found for '{sample}': {crai}")

    if not samples:
        error_exit("Samplesheet is empty (no samples found)")

    print(f"   [OK] Samplesheet valid")
    print(f"   [OK] Found {len(samples)} samples: {', '.join(samples[:5])}" +
          (f" ..." if len(samples) > 5 else ""))

    # Check pedigree file
    print("\\n[4/7] Validating pedigree file...")
    pedigree = Path("${pedigree}")
    if not pedigree.exists():
        error_exit(f"Pedigree file not found: {pedigree}")

    pedigree_samples = set()
    with open(pedigree) as f:
        for i, line in enumerate(f, start=1):
            if line.startswith('#'):
                continue
            fields = line.strip().split()
            if len(fields) < 6:
                error_exit(f"Pedigree line {i} has fewer than 6 columns")

            # Extract individual ID (column 2)
            individual = fields[1]
            pedigree_samples.add(individual)

            # Validate sex (column 5)
            sex = fields[4]
            if sex not in ('0', '1', '2'):
                warning(f"Non-standard sex value '{sex}' for {individual} (expected 0/1/2)")

            # Validate phenotype (column 6)
            phenotype_val = fields[5]
            if phenotype_val not in ('0', '1', '2', '-9'):
                warning(f"Non-standard phenotype '{phenotype_val}' for {individual}")

    print(f"   [OK] Pedigree file valid")
    print(f"   [OK] Found {len(pedigree_samples)} individuals in pedigree")

    # Check proband is in pedigree
    proband = "${params.proband}"
    if proband and proband not in pedigree_samples:
        error_exit(f"Proband '{proband}' not found in pedigree file")
    if proband:
        print(f"   [OK] Proband '{proband}' found in pedigree")

    # Check samplesheet vs pedigree consistency
    print("\\n[5/7] Checking samplesheet vs pedigree consistency...")
    sample_set = set(samples)

    # Samples in samplesheet but not in pedigree
    missing_in_pedigree = sample_set - pedigree_samples
    if missing_in_pedigree:
        warning(f"Samples in samplesheet but not in pedigree: {missing_in_pedigree}")

    # Samples in pedigree but not in samplesheet
    missing_in_samplesheet = pedigree_samples - sample_set
    if missing_in_samplesheet:
        warning(f"Samples in pedigree but not in samplesheet: {missing_in_samplesheet}")

    if not missing_in_pedigree and not missing_in_samplesheet:
        print(f"   [OK] All samples match between samplesheet and pedigree")

    # Check phenotype file
    print("\\n[6/7] Validating phenotype file...")
    phenotype = Path("${phenotype}")
    if not phenotype.exists():
        error_exit(f"Phenotype file not found: {phenotype}")

    hpo_terms = []
    with open(phenotype) as f:
        for i, line in enumerate(f, start=1):
            term = line.strip()
            if not term:
                continue
            if not term.startswith('HP:'):
                warning(f"Line {i}: '{term}' doesn't look like HPO term (expected HP:XXXXXXX)")
            hpo_terms.append(term)

    if not hpo_terms:
        error_exit("Phenotype file is empty (no HPO terms found)")

    print(f"   [OK] Phenotype file valid")
    print(f"   [OK] Found {len(hpo_terms)} HPO terms")

    # Check BED files
    print("\\n[7/7] Validating BED files...")
    target = Path("${target_regions}")
    if not target.exists():
        error_exit(f"Target regions BED file not found: {target}")

    exclusion = Path("${exclusion_bed}")
    if not exclusion.exists():
        error_exit(f"Exclusion BED file not found: {exclusion}")

    print(f"   [OK] Target regions BED: {target.name}")
    print(f"   [OK] Exclusion BED: {exclusion.name}")

    # Success
    print("\\n" + "=" * 70)
    print("SUCCESS: All validation checks passed!")
    print("=" * 70)
    """
}

process extract_fastq {
    tag "$sample"
    container 'aakrosh/bwa-suite:alpine'

    input:
    tuple val(sample), path(cram)

    output:
    tuple val(sample),path("${sample}_1.fastq.gz"),path("${sample}_2.fastq.gz") 
    
    script:
    """
    samtools sort -n --reference ${params.reference_in} -O BAM \
        -o ${sample}.sorted.bam -@ ${task.cpus} ${cram} 
    samtools fastq -@ ${task.cpus} -1 ${sample}_1.fastq.gz \
        -2 ${sample}_2.fastq.gz ${sample}.sorted.bam
    """
}

process align_fastq {
    tag "$sample"
    container "aakrosh/bwa-suite:alpine"
    publishDir 'results', mode: 'copy'

    input:
    tuple val(sample), path(fq1), path(fq2)

    output:
    tuple val(sample), path("${sample}.bam"), path("${sample}.bam.bai")

    script:
    // Split resources: more CPUs to BWA (bottleneck), rest to sort
    def bwa_threads = Math.max(1, task.cpus - 4)
    def sort_threads = Math.min(4, Math.max(1, task.cpus - bwa_threads))
    def sort_mem = "${Math.max(1, Math.floor(task.memory.toGiga() / sort_threads - 2))}G"

    """
    mkdir -p tmp_sort

    bwa mem -t ${bwa_threads} -Y \
      -R \"@RG\\tID:${sample}\\tSM:${sample}\\tPL:illumina\" \
      -K 100000000 ${params.reference} ${fq1} ${fq2} \
    | samblaster --addMateTags \
    | samtools view -u -h - \
    | samtools sort -@ ${sort_threads} -m ${sort_mem} \
        -T tmp_sort/${sample} -O BAM -o ${sample}.bam -

    samtools index -@ ${task.cpus} ${sample}.bam

    rm -rf tmp_sort
    """
}

process get_stats {
    tag "${sample}"
    container "aakrosh/bwa-suite:alpine" 
    publishDir 'results', mode: 'copy'

    input:
    tuple val(sample), path(bam), path(bai) 

    output:
    tuple val(sample), path("${sample}.stats.txt") 

    script:
    """
    alignstats -i ${bam} -o ${sample}.stats.txt -p -W -P ${task.cpus} \
        -t ${params.target_regions}
    """
}

process call_small_variants {
    tag "joint calling"
    container 'community.wave.seqera.io/library/freebayes:1.3.9--8b4bf7c06d07bf77'

    input:
    path(bams_and_bais)

    output:
    path("variants.vcf")

    script:
    def bam_files = bams_and_bais.findAll { it.name.endsWith('.bam') }
    def bam_files_str = bam_files.join(' ')

    """
    freebayes-parallel <(awk '{{print \$1":"\$2"-"\$3}}' ${params.target_regions}) \
        ${task.cpus} --genotype-qualities \
        --fasta-reference ${params.reference} ${bam_files_str}  \
    > variants.vcf
    """
}

process filter_variants {
    tag "filter"
    container "community.wave.seqera.io/library/bcftools:1.21--4335bec1d7b44d11"

    input:
    path(vcf)
    path(exclude_bed)

    output:
    path("${vcf.simpleName}.flt.vcf")

    script:
    """
    bcftools filter -i 'QUAL>1 && QUAL/INFO/AO>10 && INFO/SAF>0 && INFO/SAR>0 && INFO/RPR>1 && INFO/RPL>1' ${vcf} \
    | bcftools +setGT -- -t q -n . -i 'FMT/GQ<20' \
    | bcftools annotate -a ${exclude_bed} -c CHROM,FROM,TO \
        -h <(echo '##FILTER=<ID=BED_OVERLAP,Description="Variant overlaps BED regions">') \
        -m BED_OVERLAP -Ov -o ${vcf.simpleName}.flt.vcf
    """
}

process calculate_variant_stats {
    tag "variant stats"
    container "openjdk:21-jdk"
    publishDir 'results', mode: 'copy'

    input:
    tuple path(vcf), path(tbi)
    path(pedigree_file)

    output:
    path("variant_qc.html")

    script:
    """
    java -jar ${params.variantqc_jar} VariantQC -R ${params.reference} \
        -ped ${pedigree_file} -V ${vcf} -O variant_qc.html
    """
}

process leftalign_split {
    tag "standardize"
    container "community.wave.seqera.io/library/bcftools:1.21--4335bec1d7b44d11"
    publishDir 'results', mode: 'copy'    

    input:
    tuple path(vcf), path(tbi)

    output:
    path("${vcf.baseName}.norm.vcf")

    script:
    """
    bcftools norm -a -d all -f ${params.reference} -m -both -Ov ${vcf} \
    | bcftools view --trim-alt-alleles -Ov -o ${vcf.baseName}.norm.vcf
    """
}

process prioritize_variants {
    tag "exomiser_prioritisation"
    container "openjdk:21-jdk"
    publishDir 'results', mode: 'copy'
    
    input:
    tuple path("variants.norm.vcf.gz"), path("variants.norm.vcf.gz.tbi") 
    path(phenotype_file)  // A file containing HPO terms
    path(pedigree_file)
    
    output:
    path("exomiser_results/prioritized.html")
    path("exomiser_results/prioritized.json")
    path("exomiser_results/prioritized.genes.tsv")
    path("exomiser_results/prioritized.variants.tsv")
    path("exomiser_results/prioritized.vcf.gz")
    path("exomiser_results/prioritized.vcf.gz.tbi")
    
    script:
    // Dynamically allocate memory based on SLURM allocation
    def max_mem = task.memory ? "-Xmx${task.memory.toMega() - 2048}m" : "-Xmx4g"

    """
    # Create an Exomiser analysis YAML configuration file
    cat > analysis.yml << EOF
    analysis:
        genomeAssembly: ${params.genome_assembly}
        vcf: variants.norm.vcf.gz
        ped: ${pedigree_file}
        proband: ${params.proband}
        hpoIds: [ \$(cat ${phenotype_file} | paste -sd, -) ]
        analysisMode: PASS_ONLY
        frequencySources: [
            THOUSAND_GENOMES, 
            TOPMED, 
            UK10K, 
            GNOMAD_E_AFR, 
            GNOMAD_E_AMR, 
            GNOMAD_E_ASJ, 
            GNOMAD_E_EAS, 
            GNOMAD_E_FIN, 
            GNOMAD_E_NFE, 
            GNOMAD_E_OTH, 
            GNOMAD_E_SAS, 
            GNOMAD_G_AFR, 
            GNOMAD_G_AMR, 
            GNOMAD_G_ASJ, 
            GNOMAD_G_EAS, 
            GNOMAD_G_FIN, 
            GNOMAD_G_NFE, 
            GNOMAD_G_OTH, 
            GNOMAD_G_SAS]
        pathogenicitySources: [REVEL, MVP, SPLICE_AI, ALPHA_MISSENSE]
        inheritanceModes: {
            AUTOSOMAL_DOMINANT: 0.1,
            AUTOSOMAL_RECESSIVE_HOM_ALT: 0.1,
            AUTOSOMAL_RECESSIVE_COMP_HET: 2.0,
            X_DOMINANT: 0.1,
            X_RECESSIVE_HOM_ALT: 0.1,
            X_RECESSIVE_COMP_HET: 2.0,
            MITOCHONDRIAL: 0.2
        }
        steps: [
            variantEffectFilter: {
                remove: [
                    FIVE_PRIME_UTR_EXON_VARIANT,
                    FIVE_PRIME_UTR_INTRON_VARIANT,
                    THREE_PRIME_UTR_EXON_VARIANT,
                    THREE_PRIME_UTR_INTRON_VARIANT,
                    NON_CODING_TRANSCRIPT_EXON_VARIANT,
                    UPSTREAM_GENE_VARIANT,
                    INTERGENIC_VARIANT,
                    REGULATORY_REGION_VARIANT,
                    CODING_TRANSCRIPT_INTRON_VARIANT,
                    NON_CODING_TRANSCRIPT_INTRON_VARIANT,
                    DOWNSTREAM_GENE_VARIANT
                ]
            },
            frequencyFilter: {maxFrequency: 2.0},
            pathogenicityFilter: {keepNonPathogenic: true},
            inheritanceFilter: {},
            omimPrioritiser: {},
            priorityScoreFilter: {priorityType: HIPHIVE_PRIORITY, minPriorityScore: 0.501},
            hiPhivePrioritiser: {runParams: 'human'},
        ]
    outputOptions:
        outputContributingVariantsOnly: true
        numGenes: 30
        outputPrefix: exomiser_results/prioritized
        outputFormats: [HTML, JSON, TSV_GENE, TSV_VARIANT, VCF]
    EOF
    
    # Create the output directory
    mkdir -p exomiser_results

    # Run Exomiser CLI with dynamic memory allocation
    java -Xms2g ${max_mem} -jar ${params.exomiser_cli_jar} \
        --analysis analysis.yml \
        --spring.config.location=${params.exomiser_config}
    """
}

process partition_variants {
    tag 'partition'
    container 'community.wave.seqera.io/library/bcftools:1.21--4335bec1d7b44d11'

    input:
    tuple path(vcf), path(tbi)

    output:
    path('*.vcf')

    script:
    """
    for sample in \$(bcftools query -l ${vcf}); do
        bcftools view -c1 -s "\${sample}" -Ou --threads ${task.cpus} ${vcf} \
        | bcftools annotate -x ^FORMAT/GT,FORMAT/AD,FORMAT/DP,FORMAT/GQ -Ov \
          -o "\${sample}.tmp.vcf" --threads ${task.cpus} 

        awk 'BEGIN {OFS="\\t"}
        /^#/ {print; next}
        {
         split(\$9, fmt, ":");
         split(\$10, dat, ":");
         for (i in fmt) map[fmt[i]] = dat[i];
         new_fmt = "GT:AD:DP:GQ";
         \$9 = new_fmt; 
         \$10 = map["GT"]":"map["AD"]":"map["DP"]":"map["GQ"];
         print
        }' "\${sample}.tmp.vcf" > "\${sample}.vcf"

        rm "\${sample}.tmp.vcf"
    done
    """
}

process annotate_vep {
    tag "VEP"
    container 'ensemblorg/ensembl-vep:release_114.0'

    input:
    path(vcf)
   
    output:
    tuple val("${vcf.baseName}"), val("VEP"), path("${vcf.simpleName}.vep.vcf")

    script:
    """
    vep --offline --cache --fork ${task.cpus} \
    --dir_cache ${params.vep_cachedir} --fasta ${params.reference} \
    --use_given_ref --species homo_sapiens --assembly ${params.genome_assembly} \
    --xref_refseq --hgvs --hgvsg --canonical --symbol --distance 0 \
    --exclude_predicted --flag_pick --lookup_ref --force --input_file ${vcf} \
    --output_file ${vcf.simpleName}.vep.vcf --format vcf --vcf --no_stats \
    --numbers
    """
}

process annotate_annovar {
    tag "Annovar"

    input:
    tuple val(sample), val(type), path(vcf)

    output:
    tuple val(sample), val("ANNOVAR"), path("${sample}.*_multianno.txt")

    script:
    def suffix = params.genome_assembly == 'GRCh38' ? 'hg38' : 'hg19'
    """
    perl ${params.annovar_dir}/table_annovar.pl ${vcf} \
        ${params.annovar_database} --buildver ${suffix} --out ${sample} \
        --remove --protocol gnomad211_exome,gnomad211_genome --operation f,f \
        --vcfinput --thread ${task.cpus}
    """
}

process annotate_intervar {
    tag "Intervar"
    container 'aakrosh/intervar:latest'

    input:
    tuple val(sample), val(type), path(vcf)

    output:
    tuple val(sample), val("INTERVAR"), path("${vcf.simpleName}.*_multianno.txt.intervar")

    script:
    def suffix = params.genome_assembly == 'GRCh38' ? 'hg38' : 'hg19'
    """
    /usr/bin/python3 /Intervar.py -b ${suffix} -i ${vcf} --input_type=VCF \
        -o ${vcf.simpleName} -t ${params.intervar_database} \
        --table_annovar ${params.annovar_dir}/table_annovar.pl \
        --database_locat ${params.annovar_database} \
        --convert2annovar ${params.annovar_dir}/convert2annovar.pl \
        --annotate_variation ${params.annovar_dir}/annotate_variation.pl
    """
}

process annotate_autopvs1 {
    tag "autopvs1"
    container 'aakrosh/autopvs1:latest'    

    input:
    tuple val(sample), val(type), path(vcf)
    path(autopvs1_data)
    path(autopvs1_config)

    output:
    tuple val(sample), val("AUTOPVS1"), path("${vcf.simpleName}.autopvs1.txt")

    script:
    def suffix = params.genome_assembly == 'GRCh38' ? 'hg38' : 'hg19'
    """
    /usr/bin/python3 /autopvs1/autoPVS1_from_VEP_vcf.py \
        --genome_version ${suffix} \
        --vep_vcf ${vcf} > ${vcf.simpleName}.autopvs1.txt
    """
}

process select_clinvar {
    tag "clinvar"
    container 'pgc-images.sbgenomics.com/diskin-lab/autogvp:v1.0.5'    
    containerOptions = { "--bind ${params.autogvp_dir}:/home/rstudio/AutoGVP" }

    output:
    path("ClinVar-selected-submissions.tsv")

    script:
    """
    Rscript /home/rstudio/AutoGVP/scripts/select-clinVar-submissions.R \
        --variant_summary /home/rstudio/AutoGVP/data/variant_summary.txt.gz \
        --submission_summary /home/rstudio/AutoGVP/data/submission_summary.txt.gz \
        --outdir . \
        --conceptID_list /home/rstudio/AutoGVP/data/clinvar_cpg_concept_ids.txt \
        --conflict_res "latest"
    """
}

process run_autogvp {
    tag "autogvp"
    container 'pgc-images.sbgenomics.com/diskin-lab/autogvp:v1.0.5'    
    containerOptions = { "--bind ${params.autogvp_dir}:/home/rstudio/AutoGVP" }
    publishDir 'results', mode: 'copy'
   
    input:
    tuple val(sample), path(vep_file), path(annovar_file), 
            path(intervar_file), path(autopvs1_file)
    path(clinvar_file)

    output:
    path("${vep_file.baseName}*")

    script:
    """
    bash /home/rstudio/AutoGVP/run_autogvp.sh --workflow="custom" \
    --vcf=${vep_file} \
    --clinvar=/home/rstudio/AutoGVP/data/clinvar.vcf.gz \
    --intervar=${intervar_file} \
    --multianno=${annovar_file} \
    --autopvs1=${autopvs1_file} \
    --outdir=. \
    --out="${vep_file.baseName}" \
    --selected_clinvar_submissions=${clinvar_file} \
    --variant_summary=/home/rstudio/AutoGVP/data/variant_summary.txt.gz \
    --submission_summary=/home/rstudio/AutoGVP/data/submission_summary.txt.gz \
    --conceptIDs=/home/rstudio/AutoGVP/data/clinvar_cpg_concept_ids.txt \
    --conflict_res="latest"   
    """
}

process multiqc {
    tag "QC Report"
    container 'quay.io/biocontainers/multiqc:1.21--pyhdfd78af_0'
    publishDir 'results/multiqc', mode: 'copy'

    input:
    path('*')

    output:
    path "multiqc_report.html"
    path "multiqc_data"

    script:
    """
    multiqc . \
        --title "PedigreeVarFlow QC Report" \
        --filename multiqc_report.html \
        --force
    """
}

workflow {
    // Validate all inputs before starting the pipeline
    validate_inputs(
        file(params.samplesheet),
        file(params.reference),
        file(params.reference_in),
        file(params.pedigree),
        file(params.phenotype_file),
        file(params.target_regions),
        file(params.exclusion_bed)
    )

    Channel.fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample, file(row.cram_location)) }
        .set { input_channel }

    // extract the fastq files from the CRAM
    extract_fastq(input_channel)
        .set { fq_channel }

    // align the reads into a BAM file
    align_fastq(fq_channel)
        .set { bam_channel }

    // calculate the alignment stats
    get_stats(bam_channel)
        .set { stats_channel }

    // call variants using FreeBayes
    bam_channel
        .map { sample, bam, bai -> [bam, bai] }
        .flatten()
        .collect()
        .set { all_bams_for_calling }
    call_small_variants(all_bams_for_calling)
        .set { variants_channel }

    // filter the variants to remove low-quality variants and genotypes
    filter_variants(variants_channel, file(params.exclusion_bed))
        .set { filtered_variants_channel }

    // bgzip and index the variants
    bgzip_index_raw(filtered_variants_channel)
        .set { zipped_variants_channel }

    // collect some stats
    calculate_variant_stats(zipped_variants_channel, file(params.pedigree))

    // left-align and normalize indels, split multiallelic sites
    leftalign_split(zipped_variants_channel)
        .set { standardize_channel }

    // bzgip and index the variants
    bgzip_index_standard(standardize_channel)
        .set { zipped_standard_channel }

    // prioritize the variants using exomiser
    prioritize_variants(zipped_standard_channel,
                        file(params.phenotype_file),
                        file(params.pedigree))

    // here on, we want to calculate variant effect per sample in the VCF file
    // applying ACMG guidelines. First partition the multi-vcf into per-sample
    // VCF. Also, AUTOGVP only expects GT,AD,DP,GQ for sample
    partition_variants(zipped_standard_channel)
        .set { per_sample_vcfs_ch }

    // annotate the file using VEP
    per_sample_vcfs_ch
        .flatten()
        .set { single_vcf_ch }

    annotate_vep(single_vcf_ch)
        .set { vep_annotated_channel }

    // annotate the variants using ANNOVAR
    annotate_annovar(vep_annotated_channel)    
        .set { annovar_output_channel }

    // annotate the variants using INTERVAR
    annotate_intervar(vep_annotated_channel)
        .set { intervar_output_channel }

    // annotate the variants using AUTOPVS1
    annotate_autopvs1(vep_annotated_channel, 
                      file(params.autopvs1_data), 
                      file(params.autopvs1_config))
        .set { autopvs1_output_channel }

    // create the select clinvar file
    select_clinvar()
        .set { clinvar_channel }

    // run the AutoGVP workflow
    all_outputs = vep_annotated_channel
        .mix(annovar_output_channel)
        .mix(intervar_output_channel)
        .mix(autopvs1_output_channel)
        .map { sample, tool, file -> tuple(sample, tuple(tool, file)) }

    grouped = all_outputs
        .groupTuple()
        .map { sample, records ->
            def file_map = records.collectEntries { [(it[0]): it[1]] }
            def vep      = file_map.get('VEP')
            def annovar  = file_map.get('ANNOVAR')
            def intervar = file_map.get('INTERVAR')
            def autopvs1 = file_map.get('AUTOPVS1')

            if ([vep, annovar, intervar, autopvs1].any { it == null }) return null
            return tuple(sample, vep, annovar, intervar, autopvs1)
        }
        .filter { it != null }

    run_autogvp(grouped, clinvar_channel)

    // Collect QC outputs for MultiQC report
    stats_channel
        .map { sample, stats -> stats }
        .collect()
        .set { qc_files }

    multiqc(qc_files)

}
