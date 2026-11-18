process ALIGN {

	tag "mapping on $name using $task.cpus CPUs and $task.memory memory"
        publishDir "${params.outDirectory}/${sample.run}/mapped/", mode:'copy'

        label "l_cpu"
        label "l_mem"

        shell:
       '/bin/bash'
   
        input:
        tuple val(name), val(sample), path(fwd), path(rev)

        output:
        tuple val(name), val(sample), path("${name}.Aligned.sortedByCoord.out.bam"), path("${name}.Aligned.sortedByCoord.out.bam.bai"), path("${name}.Chimeric.out.junction")

        script:
        """
       source activate star
       STAR --runMode alignReads --genomeDir ${params.STAR_INDEX} \
            --readFilesIn $fwd $rev --readFilesCommand zcat \
            --sjdbOverhang 100  \
            --outFileNamePrefix ${name}. \
            --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
            --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 1.0 --outFilterMismatchNoverLmax 0.3 \
            --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 \
            --outFilterMatchNmin 0.66 --outFilterScoreMinOverLread 0.66 --outFilterMatchNminOverLread 0.66 \
            --outSAMattrRGline ID:${name} PL:Illumina PU:${name} SM:${name} \
            --outSAMunmapped Within --outFilterType BySJout --outSAMattributes All \
            --outWigStrand Stranded --quantMode GeneCounts TranscriptomeSAM --sjdbScore 1 --twopassMode Basic \
            --outMultimapperOrder Random \
            --outSAMtype BAM SortedByCoordinate \
            --chimSegmentMin 10 \
            --chimOutType Junctions WithinBAM HardClip \
            --chimJunctionOverhangMin 10 \
            --chimScoreDropMax 30 \
            --chimScoreJunctionNonGTAG 0 \
            --chimScoreSeparation 1 \
            --chimSegmentReadGapMax 3 \
            --chimMultimapNmax 50

     samtools index ${name}.Aligned.sortedByCoord.out.bam

         """
}

process ARRIBAfusion {

        tag "mapping on $name using $task.cpus CPUs and $task.memory memory"
        publishDir "${params.outDirectory}/${sample.run}/arriba/", mode:'copy'

        label "m_cpu"
        label "m_mem"

        shell:
       '/bin/bash'

        input:
        tuple val(name), val(sample), path("${name}.Aligned.sortedByCoord.out.bam"), path("${name}.Aligned.sortedByCoord.out.bam.bai"), path("${name}.Chimeric.out.junction")

        output:
        tuple val(name), val(sample), path("${name}.fusion.tsv"), path("${name}.fusions.discarded.tsv")

        script:
        """
        source activate arriba
        arriba -x ${name}.Aligned.sortedByCoord.out.bam -a ${params.ref}.fa -g ${params.GTF} -b ${params.blacklist} -k ${params.known_fusions} -p ${params.protein_domains} -o ${name}.fusion.tsv -O ${name}.fusions.discarded.tsv

        """
}


process STARfusion {

        tag "mapping on $name using $task.cpus CPUs and $task.memory memory"
        publishDir "${params.outDirectory}/${sample.run}/starfusion/", mode:'copy'

        label "l_cpu"
        label "l_mem"

        shell:
       '/bin/bash'

        input:
        tuple val(name), val(sample), path(fwd), path(rev)

        output:
        tuple val(name), val(sample), path("star-fusion-${name}")

        script:
        """
        echo STARFUSIOn $name

        export TMPDIR=\${PWD}/tmp
        mkdir -p \$TMPDIR

        source activate star-fusion
        STAR-Fusion --genome_lib_dir ${params.GRCh38_gencode_v37_CTAT_lib} --left_fq $fwd --right_fq $rev --output_dir \${PWD}/star-fusion-${name}
        """
} 


process Featurecounts {

        tag "mapping on $name using $task.cpus CPUs and $task.memory memory"
        publishDir "${params.outDirectory}/${sample.run}/Featurecounts/", mode:'copy'

        label "m_cpu"
        label "m_mem"

        shell:
       '/bin/bash'

        input:
        tuple val(name), val(sample), path("${name}.Aligned.sortedByCoord.out.bam"), path("${name}.Aligned.sortedByCoord.out.bam.bai"), path("${name}.Chimeric.out.junction")

        output:
        tuple val(name), val(sample), path("featurecounts-${name}")

        script:
        """
        export TMPDIR=\${PWD}/tmp
        mkdir -p \$TMPDIR
 
        source activate rna
        featureCounts -M -t exon -g gene_id -O -s 0 -T 12 -Q 0 -d 1 -D 1000 -a ${params.GTF} -o \${PWD}/featurecounts-${name} ${name}.Aligned.sortedByCoord.out.bam # -f parametr pak da vysledek na jednotlive exony, pro jednotlive  transkripty pak toto -g transcript_id, pripadne -s 0 je unstranded, ignoruje se smer a pocitaji se vsechna cteni

        """
}

process MergeFeatureCounts {
        
        publishDir "${params.outDirectory}/${run}/Featurecounts/", mode:'copy'
  
        input:
        tuple val(run), path(featurecounts_files)

        output:
        tuple val(run), path("featurecounts-ALL.txt"), path("featurecounts-ALL-genes.txt")

        script:
        """
        source activate arriba
        first=\$(echo ${featurecounts_files} | awk '{print \$1}')

        tail -n +2 "\$first" | cut -f1 > genes.txt

        cp genes.txt featurecounts-ALL.txt

        for f in ${featurecounts_files}; do
        tail -n +2 "\$f" | cut -f7 > counts.tmp
        paste featurecounts-ALL.txt counts.tmp > merged.tmp
        mv merged.tmp featurecounts-ALL.txt
        done 

        cp ${params.GTF_engs_gene} Homo_sapiens.GRCh38.95-gene-engs.txt

        Rscript ${projectDir}/scripts/engs-to-gene.R

        """
}

process COVERAGE1 {
        tag "COVERAGE1 on $name"
        publishDir "${params.outDirectory}/${sample.run}/coverage/", mode:'copy'
        container "staphb/samtools:1.20"

        label "l_cpu"
        label "l_mem"


        input:
         tuple val(name), val(sample), path("${name}.Aligned.sortedByCoord.out.bam"), path("${name}.Aligned.sortedByCoord.out.bam.bai"), path("${name}.Chimeric.out.junction")

        output:
        tuple val(name), val(sample), path("${name}.coveragefin.txt")

        script:
        """
        echo COVERAGE1 $name
        samtools bedcov ${params.varbed1} ${name}.Aligned.sortedByCoord.out.bam -d 20 > ${name}.COV
        awk '{print \$5/(\$3-\$2)}'  ${name}.COV >  ${name}.COV-mean
        awk '{print (\$6/(\$3-\$2))*100"%"}' ${name}.COV > ${name}-procento-nad-20
        paste ${name}.COV-mean ${name}-procento-nad-20 > vysledek
        echo "chr" "start" "stop" "name" ${name}.COV-mean ${name}-procento-nad-20 > hlavicka
        sed -i 's/ /\t/'g hlavicka
        paste ${params.varbed1} vysledek > coverage
        cat hlavicka coverage > ${name}.coveragefin.txt
        sed -i -e "s/\r//g" ${name}.coveragefin.txt
        """
}

process COMBINECOVERAGEMEAN {
    tag "COMBINECOVERAGEMEAN"

    input:
    tuple val(run_name), path(coverage_files)

    publishDir { "${params.outDirectory}/${run_name}/coverage/" }, mode: 'copy'

    output:
    tuple val(run_name), path("coveragemeanALL")


    script:
    """
    i=0
    tmpfiles=""
    for file in ${coverage_files}; do
        awk '{print \$5}' "\$file" > column_\$i.txt
        tmpfiles="\$tmpfiles column_\$i.txt"
        i=\$((i+1))
    done
    echo "chr" "start" "stop" "name" > hlavicka
    sed -i 's/ /\t/'g hlavicka
    cat hlavicka ${params.varbed1} > bedshlavickou
    paste  bedshlavickou \$tmpfiles > coveragemeanALL
    """
}


process GATK {
         container 'broadinstitute/gatk:4.6.0.0'
         tag "variantcalling on $name using $task.cpus CPUs and $task.memory memory"
         publishDir "${params.outDirectory}/${sample.run}/varianty/", mode:'copy'


         input:
         tuple val(name), val(sample), path("${name}.Aligned.sortedByCoord.out.bam"), path("${name}.Aligned.sortedByCoord.out.bam.bai"), path("${name}.Chimeric.out.junction")

         output:
         tuple val(name), val(sample), path("${name}.vcf")

    script:
    """
    gatk --java-options "-Xmx4g" SplitNCigarReads -R ${params.ref}.fa -I ${name}.Aligned.sortedByCoord.out.bam -O ${name}.split.bam
    gatk --java-options "-Xmx4g" HaplotypeCaller -R ${params.ref}.fa -I ${name}.split.bam -L ${params.varbed1} --dont-use-soft-clipped-bases true -A StrandBiasBySample -minimum-mapping-quality 0 --mapping-quality-threshold-for-genotyping 0 --enable-dynamic-read-disqualification-for-genotyping true --flow-filter-alleles-qual-threshold 0 -O ${name}.vcf

    """
}

process VAFaNORMALIZACE {
        tag "VAFaNORMALIZACE on $name"
        publishDir "${params.outDirectory}/${sample.run}/varianty/", mode:'copy'
        container "xreblova/bcftoolsbgziptabix:2"

        input:
        tuple val(name), val(sample), path(gatk)

        output:
        tuple val(name), val(sample), path("${name}.norm.vcf.gz"), path("${name}.norm.vcf.gz.tbi")

        script:
        """
        echo VAFaNORMALIZACE $name

        bcftools +fill-tags $gatk -Ob -o ${name}.pom2.bcf -- -t FORMAT/VAF
        bcftools convert -O v -o ${name}.vaf.vcf ${name}.pom2.bcf

        bcftools norm -f ${params.ref}.fa -m -both ${name}.vaf.vcf -o ${name}.norm.vcf
        bgzip ${name}.norm.vcf
        tabix ${name}.norm.vcf.gz

        """
}


process LOH {
        tag "LOH"
        publishDir { "${params.outDirectory}/${sample.run}/LOH/fily" }, mode: 'copy'
        container "xreblova/bcftoolsbgziptabix:2"

        input:
        tuple val(name), val(sample), path("${name}.norm.vcf.gz"), path("${name}.norm.vcf.gz.tbi")

        output:
        tuple val(name), val(sample), path("${name}.chr1.smooth"), path("${name}.chr2.smooth"), path("${name}.chr3.smooth"),path("${name}.chr4.smooth"),path("${name}.chr5.smooth"),path("${name}.chr6.smooth"),path("${name}.chr7.smooth"),path("${name}.chr8.smooth"),path("${name}.chr9.smooth"),path("${name}.chr10.smooth"),path("${name}.chr11.smooth"),path("${name}.chr12.smooth"),path("${name}.chr13.smooth"),path("${name}.chr14.smooth"),path("${name}.chr15.smooth"),path("${name}.chr16.smooth"),path("${name}.chr17.smooth"),path("${name}.chr18.smooth"),path("${name}.chr19.smooth"),path("${name}.chr20.smooth"),path("${name}.chr21.smooth"),path("${name}.chr22.smooth"), path("${name}.chrX.smooth")

        script:
        """
        bcftools query   -f '%CHROM\t%POS\t[%GT]\t[%VAF]\n' ${name}.norm.vcf.gz > ${name}.baf.tsv

        awk '{print "chr"\$0}' ${name}.baf.tsv > ${name}.baf-chr.tsv

       grep -w chr1 ${name}.baf-chr.tsv | awk '{print \$1"\t"\$2"\t"\$3"\t"(\$4>0.5?\$4-0.5:0.5-\$4)}' | awk '{bin=int(\$2/500000); sum[bin]+=\$4; n[bin]++} END{for(i in sum)print i*500000, sum[i]/n[i]}'  | sort -n > ${name}.chr1.smooth
       grep -w chr2 ${name}.baf-chr.tsv | awk '{print \$1"\t"\$2"\t"\$3"\t"(\$4>0.5?\$4-0.5:0.5-\$4)}' | awk '{bin=int(\$2/500000); sum[bin]+=\$4; n[bin]++} END{for(i in sum)print i*500000, sum[i]/n[i]}'  | sort -n > ${name}.chr2.smooth
       grep -w chr3 ${name}.baf-chr.tsv | awk '{print \$1"\t"\$2"\t"\$3"\t"(\$4>0.5?\$4-0.5:0.5-\$4)}' | awk '{bin=int(\$2/500000); sum[bin]+=\$4; n[bin]++} END{for(i in sum)print i*500000, sum[i]/n[i]}'  | sort -n > ${name}.chr3.smooth
       grep -w chr4 ${name}.baf-chr.tsv | awk '{print \$1"\t"\$2"\t"\$3"\t"(\$4>0.5?\$4-0.5:0.5-\$4)}' | awk '{bin=int(\$2/500000); sum[bin]+=\$4; n[bin]++} END{for(i in sum)print i*500000, sum[i]/n[i]}'  | sort -n > ${name}.chr4.smooth
       grep -w chr5 ${name}.baf-chr.tsv | awk '{print \$1"\t"\$2"\t"\$3"\t"(\$4>0.5?\$4-0.5:0.5-\$4)}' | awk '{bin=int(\$2/500000); sum[bin]+=\$4; n[bin]++} END{for(i in sum)print i*500000, sum[i]/n[i]}'  | sort -n > ${name}.chr5.smooth
       grep -w chr6 ${name}.baf-chr.tsv | awk '{print \$1"\t"\$2"\t"\$3"\t"(\$4>0.5?\$4-0.5:0.5-\$4)}' | awk '{bin=int(\$2/500000); sum[bin]+=\$4; n[bin]++} END{for(i in sum)print i*500000, sum[i]/n[i]}'  | sort -n > ${name}.chr6.smooth
       grep -w chr7 ${name}.baf-chr.tsv | awk '{print \$1"\t"\$2"\t"\$3"\t"(\$4>0.5?\$4-0.5:0.5-\$4)}' | awk '{bin=int(\$2/500000); sum[bin]+=\$4; n[bin]++} END{for(i in sum)print i*500000, sum[i]/n[i]}'  | sort -n > ${name}.chr7.smooth
       grep -w chr8 ${name}.baf-chr.tsv | awk '{print \$1"\t"\$2"\t"\$3"\t"(\$4>0.5?\$4-0.5:0.5-\$4)}' | awk '{bin=int(\$2/500000); sum[bin]+=\$4; n[bin]++} END{for(i in sum)print i*500000, sum[i]/n[i]}'  | sort -n > ${name}.chr8.smooth
       grep -w chr9 ${name}.baf-chr.tsv | awk '{print \$1"\t"\$2"\t"\$3"\t"(\$4>0.5?\$4-0.5:0.5-\$4)}' | awk '{bin=int(\$2/500000); sum[bin]+=\$4; n[bin]++} END{for(i in sum)print i*500000, sum[i]/n[i]}'  | sort -n > ${name}.chr9.smooth
       grep -w chr10 ${name}.baf-chr.tsv | awk '{print \$1"\t"\$2"\t"\$3"\t"(\$4>0.5?\$4-0.5:0.5-\$4)}' | awk '{bin=int(\$2/500000); sum[bin]+=\$4; n[bin]++} END{for(i in sum)print i*500000, sum[i]/n[i]}'  | sort -n > ${name}.chr10.smooth
       grep -w chr11 ${name}.baf-chr.tsv | awk '{print \$1"\t"\$2"\t"\$3"\t"(\$4>0.5?\$4-0.5:0.5-\$4)}' | awk '{bin=int(\$2/500000); sum[bin]+=\$4; n[bin]++} END{for(i in sum)print i*500000, sum[i]/n[i]}'  | sort -n > ${name}.chr11.smooth
       grep -w chr12 ${name}.baf-chr.tsv | awk '{print \$1"\t"\$2"\t"\$3"\t"(\$4>0.5?\$4-0.5:0.5-\$4)}' | awk '{bin=int(\$2/500000); sum[bin]+=\$4; n[bin]++} END{for(i in sum)print i*500000, sum[i]/n[i]}'  | sort -n > ${name}.chr12.smooth
       grep -w chr13 ${name}.baf-chr.tsv | awk '{print \$1"\t"\$2"\t"\$3"\t"(\$4>0.5?\$4-0.5:0.5-\$4)}' | awk '{bin=int(\$2/500000); sum[bin]+=\$4; n[bin]++} END{for(i in sum)print i*500000, sum[i]/n[i]}'  | sort -n > ${name}.chr13.smooth
       grep -w chr14 ${name}.baf-chr.tsv | awk '{print \$1"\t"\$2"\t"\$3"\t"(\$4>0.5?\$4-0.5:0.5-\$4)}' | awk '{bin=int(\$2/500000); sum[bin]+=\$4; n[bin]++} END{for(i in sum)print i*500000, sum[i]/n[i]}'  | sort -n > ${name}.chr14.smooth
       grep -w chr15 ${name}.baf-chr.tsv | awk '{print \$1"\t"\$2"\t"\$3"\t"(\$4>0.5?\$4-0.5:0.5-\$4)}' | awk '{bin=int(\$2/500000); sum[bin]+=\$4; n[bin]++} END{for(i in sum)print i*500000, sum[i]/n[i]}'  | sort -n > ${name}.chr15.smooth
       grep -w chr16 ${name}.baf-chr.tsv | awk '{print \$1"\t"\$2"\t"\$3"\t"(\$4>0.5?\$4-0.5:0.5-\$4)}' | awk '{bin=int(\$2/500000); sum[bin]+=\$4; n[bin]++} END{for(i in sum)print i*500000, sum[i]/n[i]}'  | sort -n > ${name}.chr16.smooth
       grep -w chr17 ${name}.baf-chr.tsv | awk '{print \$1"\t"\$2"\t"\$3"\t"(\$4>0.5?\$4-0.5:0.5-\$4)}' | awk '{bin=int(\$2/500000); sum[bin]+=\$4; n[bin]++} END{for(i in sum)print i*500000, sum[i]/n[i]}'  | sort -n > ${name}.chr17.smooth
       grep -w chr18 ${name}.baf-chr.tsv | awk '{print \$1"\t"\$2"\t"\$3"\t"(\$4>0.5?\$4-0.5:0.5-\$4)}' | awk '{bin=int(\$2/500000); sum[bin]+=\$4; n[bin]++} END{for(i in sum)print i*500000, sum[i]/n[i]}'  | sort -n > ${name}.chr18.smooth
       grep -w chr19 ${name}.baf-chr.tsv | awk '{print \$1"\t"\$2"\t"\$3"\t"(\$4>0.5?\$4-0.5:0.5-\$4)}' | awk '{bin=int(\$2/500000); sum[bin]+=\$4; n[bin]++} END{for(i in sum)print i*500000, sum[i]/n[i]}'  | sort -n > ${name}.chr19.smooth
       grep -w chr20 ${name}.baf-chr.tsv  |awk '{print \$1"\t"\$2"\t"\$3"\t"(\$4>0.5?\$4-0.5:0.5-\$4)}' | awk '{bin=int(\$2/500000); sum[bin]+=\$4; n[bin]++} END{for(i in sum)print i*500000, sum[i]/n[i]}'  | sort -n > ${name}.chr20.smooth
       grep -w chr21 ${name}.baf-chr.tsv | awk '{print \$1"\t"\$2"\t"\$3"\t"(\$4>0.5?\$4-0.5:0.5-\$4)}' | awk '{bin=int(\$2/500000); sum[bin]+=\$4; n[bin]++} END{for(i in sum)print i*500000, sum[i]/n[i]}'  | sort -n > ${name}.chr21.smooth
       grep -w chr22 ${name}.baf-chr.tsv | awk '{print \$1"\t"\$2"\t"\$3"\t"(\$4>0.5?\$4-0.5:0.5-\$4)}' | awk '{bin=int(\$2/500000); sum[bin]+=\$4; n[bin]++} END{for(i in sum)print i*500000, sum[i]/n[i]}'  | sort -n > ${name}.chr22.smooth
       grep -w chrX ${name}.baf-chr.tsv | awk '{print \$1"\t"\$2"\t"\$3"\t"(\$4>0.5?\$4-0.5:0.5-\$4)}' | awk '{bin=int(\$2/500000); sum[bin]+=\$4; n[bin]++} END{for(i in sum)print i*500000, sum[i]/n[i]}'  | sort -n > ${name}.chrX.smooth



       """
}

process GNUPLOT {
        tag "GNUPLOT"
        publishDir { "${params.outDirectory}/${sample.run}/LOH/" }, mode: 'copy'
        container "xreblova/gnuplot:1"

        input:
        tuple val(name), val(sample), path("${name}.chr1.smooth"), path("${name}.chr2.smooth"), path("${name}.chr3.smooth"),path("${name}.chr4.smooth"),path("${name}.chr5.smooth"),path("${name}.chr6.smooth"),path("${name}.chr7.smooth"),path("${name}.chr8.smooth"),path("${name}.chr9.smooth"),path("${name}.chr10.smooth"),path("${name}.chr11.smooth"),path("${name}.chr12.smooth"),path("${name}.chr13.smooth"),path("${name}.chr14.smooth"),path("${name}.chr15.smooth"),path("${name}.chr16.smooth"),path("${name}.chr17.smooth"),path("${name}.chr18.smooth"),path("${name}.chr19.smooth"),path("${name}.chr20.smooth"),path("${name}.chr21.smooth"),path("${name}.chr22.smooth"), path("${name}.chrX.smooth")

        output:
        tuple val(name), val(sample), path("${name}.loh.png")

        script:
        """
        echo 'set term png size 2000,900' > ${name}.gnuplot.script
        echo 'set output "${name}.loh.png"' >> ${name}.gnuplot.script
        echo 'set xtics font "Arial,9"' >> ${name}.gnuplot.script
        echo 'set multiplot layout 6,4' >> ${name}.gnuplot.script

        echo 'plot "${name}.chr1.smooth" using 1:2 w l title "CHR1"' >> ${name}.gnuplot.script
        echo 'plot "${name}.chr2.smooth" using 1:2 w l title "CHR2"' >> ${name}.gnuplot.script
        echo 'plot "${name}.chr3.smooth" using 1:2 w l title "CHR3"' >> ${name}.gnuplot.script
        echo 'plot "${name}.chr4.smooth" using 1:2 w l title "CHR4"' >> ${name}.gnuplot.script
        echo 'plot "${name}.chr5.smooth" using 1:2 w l title "CHR5"' >> ${name}.gnuplot.script
        echo 'plot "${name}.chr6.smooth" using 1:2 w l title "CHR6"' >> ${name}.gnuplot.script
        echo 'plot "${name}.chr7.smooth" using 1:2 w l title "CHR7"' >> ${name}.gnuplot.script
        echo 'plot "${name}.chr8.smooth" using 1:2 w l title "CHR8"' >> ${name}.gnuplot.script
        echo 'plot "${name}.chr9.smooth" using 1:2 w l title "CHR9"' >> ${name}.gnuplot.script
        echo 'plot "${name}.chr10.smooth" using 1:2 w l title "CHR10"' >> ${name}.gnuplot.script
        echo 'plot "${name}.chr11.smooth" using 1:2 w l title "CHR11"' >> ${name}.gnuplot.script
        echo 'plot "${name}.chr12.smooth" using 1:2 w l title "CHR12"' >> ${name}.gnuplot.script
        echo 'plot "${name}.chr13.smooth" using 1:2 w l title "CHR13"' >> ${name}.gnuplot.script
        echo 'plot "${name}.chr14.smooth" using 1:2 w l title "CHR14"' >> ${name}.gnuplot.script
        echo 'plot "${name}.chr15.smooth" using 1:2 w l title "CHR15"' >> ${name}.gnuplot.script
        echo 'plot "${name}.chr16.smooth" using 1:2 w l title "CHR16"' >> ${name}.gnuplot.script
        echo 'plot "${name}.chr17.smooth" using 1:2 w l title "CHR17"' >> ${name}.gnuplot.script
        echo 'plot "${name}.chr18.smooth" using 1:2 w l title "CHR18"' >> ${name}.gnuplot.script
        echo 'plot "${name}.chr19.smooth" using 1:2 w l title "CHR19"' >> ${name}.gnuplot.script
        echo 'plot "${name}.chr20.smooth" using 1:2 w l title "CHR20"' >> ${name}.gnuplot.script
        echo 'plot "${name}.chr21.smooth" using 1:2 w l title "CHR21"' >> ${name}.gnuplot.script
        echo 'plot "${name}.chr22.smooth" using 1:2 w l title "CHR22"' >> ${name}.gnuplot.script
        echo 'plot "${name}.chrX.smooth" using 1:2 w l title "CHRX"' >> ${name}.gnuplot.script

        echo 'unset multiplot' >> ${name}.gnuplot.script

       gnuplot ${name}.gnuplot.script
 
      """
}


workflow {
        rawfastq = Channel.fromPath("${params.homeDir}/samplesheet.csv")
    .splitCsv(header: true)
    .map { row ->
        def baseDir = new File("${params.baseDir}")
        def runDir = baseDir.listFiles(new FilenameFilter() {
            public boolean accept(File dir, String name) {
                return name.endsWith(row.run)
            }
        })[0] //get the real folderName that has prepended date

        def fileR1 = file("${runDir}/processed_fastq/${row.name}_R1.fastq.gz", checkIfExists: true)
        def fileR2 = file("${runDir}/processed_fastq/${row.name}_R2.fastq.gz", checkIfExists: true)

                def meta = [name: row.name, run: row.run]
        [
            meta.name,
            meta,
            fileR1,
            fileR2,
                ]
    }
     . view()

aligned = ALIGN(rawfastq)
starfusion = STARfusion(rawfastq)
arriba = ARRIBAfusion(aligned)

exprese = Featurecounts(aligned)
exprese_collected = exprese
    .map { name, sample, f -> tuple(sample.run, file(f)) }
    .groupTuple() // groups by sample.run automatically!
spojeni = MergeFeatureCounts(exprese_collected)

coverage_results = COVERAGE1(aligned)

coverage_files_collected = coverage_results
    .map { name, sample, f -> tuple(sample.run, file(f)) }
    .groupTuple() // groups by sample.run automatically!
finalcoverage = COMBINECOVERAGEMEAN(coverage_files_collected)

variantcalling = GATK(aligned)
normalizovany = VAFaNORMALIZACE(variantcalling)


loh = LOH(normalizovany)
graf = GNUPLOT(loh)



}
