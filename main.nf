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
        publishDir "${params.outDirectory}/${sample.run}/starfusion/", mode:'copy'

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
}
