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
            --sjdbOverhang 100 --sjdbGTFfile ${params.GTF} \
            --outFileNamePrefix ${name}. \
            --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
            --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 1.0 --outFilterMismatchNoverLmax 0.05 \
            --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 \
            --outFilterMatchNmin 0 --outFilterScoreMinOverLread 0.66 --outFilterMatchNminOverLread 0.66 \
            --chimSegmentMin 12 --chimJunctionOverhangMin 12 \
            --chimOutType Junctions WithinBAM SoftClip --chimScoreMin 1 \
            --chimScoreDropMax 30 --chimScoreSeparation 10 --chimSegmentReadGapMax 3 \
            --outSAMattrRGline ID:${name} PL:Illumina PU:${name} SM:${name} \
            --outSAMunmapped Within --outFilterType BySJout --outSAMattributes All \
            --outWigStrand Stranded --quantMode GeneCounts TranscriptomeSAM --sjdbScore 1 --twopassMode Basic \
            --outMultimapperOrder Random \
            --outSAMtype BAM SortedByCoordinate

     samtools index ${name}.Aligned.sortedByCoord.out.bam

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
}
