#!/usr/bin/env nextflow

params.sample_id = "no_id_provided"
params.sample_gvcf = null
params.interval_vcf = null

process genotype_gvcf {
    container 'intelliseqngs/gatk-4.3.0.0-grch38-no-alt:1.0.0'
    
    input:
    path sample_gvcf
    path interval_vcf

    output:
    path "${params.sample_id}_genotyped-by-vcf.vcf.gz", emit: genotyped_vcf
    path "${params.sample_id}_genotyped-by-vcf.vcf.gz.tbi", emit: genotyped_vcf_index

    script:
    """
    set -e
    mkdir sample-files intervals

    # Find reference FASTA in the Docker image
    ref_fasta=\$(ls /resources/reference-genomes/*/*.fa)

    # Index input files
    tabix -p vcf ${sample_gvcf}
    tabix -p vcf ${interval_vcf}
    
    # Split intervals and sample gvcf by chromosome
    chroms=\$(tabix -l ${interval_vcf})
    for chrom in \$chroms; do
        tabix -h ${interval_vcf} \$chrom | bgzip > intervals/\${chrom}.intervals.vcf.gz
        tabix -h ${sample_gvcf} \$chrom | bgzip > sample-files/\${chrom}.sample.gvcf.gz
        tabix -p vcf intervals/\${chrom}.intervals.vcf.gz
        tabix -p vcf sample-files/\${chrom}.sample.gvcf.gz
    done

    # Genotype GVCFs
    for chrom in \$chroms; do
        gatk --java-options "-Xmx9g -Xms4g" GenotypeGVCFs \
            --allow-old-rms-mapping-quality-annotation-data \
            --lenient \
            --include-non-variant-sites \
            --intervals intervals/\${chrom}.intervals.vcf.gz \
            --variant sample-files/\${chrom}.sample.gvcf.gz \
            --output sample-files/\${chrom}.genotyped-sample.vcf.gz \
            --reference \$ref_fasta
    done

    # Get sample name
    sample_name=\$(bcftools query -l ${sample_gvcf})

    # Merge, normalize, and process VCFs
    for chrom in \$chroms; do
        bcftools merge -m all intervals/\${chrom}.intervals.vcf.gz sample-files/\${chrom}.genotyped-sample.vcf.gz \
            -R intervals/\${chrom}.intervals.vcf.gz -o sample-files/\${chrom}.merged.vcf.gz -O z

        bcftools norm -m -any -f \$ref_fasta sample-files/\${chrom}.merged.vcf.gz \
            | bcftools view -s \$sample_name -o sample-files/\${chrom}.genotyped-by-vcf.vcf.gz -O z

        tabix -p vcf sample-files/\${chrom}.genotyped-by-vcf.vcf.gz
    done

    # Concatenate chromosome-wise genotyped VCF files
    gatk --java-options "-Xmx9g -Xms4g" MergeVcfs \
        \$(ls sample-files/*.genotyped-by-vcf.vcf.gz | sed 's/^/-I /') \
        -O ${params.sample_id}_genotyped-by-vcf_concat.vcf.gz

    # Set missing genotypes to ./.
    bcftools filter -e 'GT="ref" & FORMAT/DP=0' --set-GTs . \
        ${params.sample_id}_genotyped-by-vcf_concat.vcf.gz \
        -o ${params.sample_id}_genotyped-by-vcf.vcf.gz -Oz
    
    tabix -p vcf ${params.sample_id}_genotyped-by-vcf.vcf.gz
    """
}

workflow {
    sample_gvcf = channel.fromPath(params.sample_gvcf)
    interval_vcf = channel.fromPath(params.interval_vcf)

    genotype_gvcf(sample_gvcf, interval_vcf)
}
