#!/usr/bin/env nextflow

params.sample_id = "no_id_provided"
params.g_vcf_gz = null
params.interval_vcf_gz = null
params.interval_vcf_gz_tbi = null

process genotype_gvcf {
    container 'intelliseqngs/gatk-4.3.0.0-grch38-no-alt:1.0.0'
    
    input:
    path g_vcf_gz
    path interval_vcf_gz
    path interval_vcf_gz_tbi

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
    tabix -p vcf ${g_vcf_gz}

    # Get sample name
    sample_name=\$(bcftools query -l ${g_vcf_gz})

    # Genotype GVCFs
    gatk --java-options "-Xmx12g -Xms12g" GenotypeGVCFs \
        --allow-old-rms-mapping-quality-annotation-data \
        --lenient \
        --include-non-variant-sites \
        --intervals ${interval_vcf_gz} \
            --variant ${g_vcf_gz} \
            --output sample-files/\${sample_name}.genotyped-sample.vcf.gz \
            --reference \${ref_fasta}

    tabix -p vcf sample-files/\${sample_name}.genotyped-sample.vcf.gz
    bcftools norm -m -any -f \${ref_fasta} sample-files/\${sample_name}.genotyped-sample.vcf.gz \
        | bcftools view -s \${sample_name} -o sample-files/\${sample_name}.normalized-sample.vcf.gz -O z

    tabix -p vcf sample-files/\${sample_name}.normalized-sample.vcf.gz

    # Set missing genotypes to ./.
    bcftools filter -e 'GT="ref" & FORMAT/DP=0' --set-GTs . \
        sample-files/\${sample_name}.normalized-sample.vcf.gz \
        -o ${params.sample_id}_genotyped-by-vcf.vcf.gz -Oz
    
    tabix -p vcf ${params.sample_id}_genotyped-by-vcf.vcf.gz
    """
}

workflow {
    g_vcf_gz = channel.fromPath(params.g_vcf_gz)
    interval_vcf_gz = channel.fromPath(params.interval_vcf_gz)
    interval_vcf_gz_tbi = channel.fromPath(params.interval_vcf_gz_tbi)

    genotype_gvcf(g_vcf_gz, interval_vcf_gz, interval_vcf_gz_tbi)
}
