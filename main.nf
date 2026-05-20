process ALIGN {

	tag "mapping on $name using $task.cpus CPUs and $task.memory memory"
        publishDir "${params.outDirectory}/${sample.run}/mapped/", mode:'copy'

        label "l_cpu"
        label "l_mem"


        input:
        tuple val(name), val(sample), path(fwd), path(rev)

        output:
        tuple val(name), val(sample), path("${name}.bam"), path("${name}.bai")

        script:
        """

       STAR --runMode alignReads --genomeDir ${params.STAR} \
            --readFilesIn <(gunzip -c $fwd) <(gunzip -c $rev) \
            --sjdbOverhang 100 --sjdbGTFfile ${params.GTF} \
            --outFileNamePrefix ${name}.bam
            --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
            --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 1.0 --outFilterMismatchNoverLmax 0.05 \
            --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 \
            --outFilterMatchNmin 0 --outFilterScoreMinOverLread 0.66 --outFilterMatchNminOverLread 0.66 \
            --outSAMheaderHD @HD VN:1.4 SO:coordinate --chimSegmentMin 30 --chimOutType SeparateSAMold \
            --outSAMattrRGline ID:${name} PL:Illumina PU:${name} SM:${name} \
            --outSAMunmapped Within --outFilterType BySJout --outSAMattributes All \
            --outWigStrand Stranded --quantMode GeneCounts TranscriptomeSAM --sjdbScore 1 --twopassMode None \
            --outMultimapperOrder Random --outSAMtype BAM SortedByCoordinate

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
