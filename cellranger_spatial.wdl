version 1.0

# Task for running Cell Ranger for spatial transcriptomics
task cellranger_spatial {
    input {
        # Required inputs
        File fastq_tar       # Tarball of gzipped FASTQ files
        File ref_tar         # Tarball of gzipped reference genome
        String sample_id     # Sample ID for Cell Ranger
        String spacename     # Name of the spatial tissue section
        Int memory_count     # Memory to allocate for this task
        Int cpu_count        # Number of CPUs to allocate for this task
        Int storage_count    # Disk space required for this task

        # Optional inputs
        Boolean? nosecondary # Flag to skip secondary analysis
        Int? localcores      # Number of local cores to use
        String? localmem     # Local memory to allocate
        String? localvmem    # Local virtual memory to allocate
    }

    parameter_meta {
        fastq_tar: "Tarball of gzipped FASTQ files"
        ref_tar: "Tarball of gzipped reference genome"
        sample_id: "Sample ID for Cell Ranger"
        spacename: "Name of the spatial tissue section"
        nosecondary: "Flag to skip secondary analysis"
        localcores: "Number of local cores"
        localmem: "Local memory"
        localvmem: "Local virtual memory"
    }

    command {
        set -exo pipefail
        
        # Unpack FASTQ files
        mkdir fastq_bundle
        echo ~{fastq_tar}
        tar -xzf ~{fastq_tar} -C fastq_bundle --no-same-owner
        fastq_path=$(realpath fastq_bundle)
        echo $fastq_path
        mv fastq_bundle/*/* fastq_bundle
        fastq_file_path="$(dirname $fastq_path)"
        echo $fastq_file_path

        # Unpack reference genome
        mkdir ref_bundle
        echo ~{ref_tar}
        tar -xzf ~{ref_tar} -C ref_bundle --no-same-owner
        ref_path=$(realpath ref_bundle)
        echo $ref_path
        mv ref_bundle/*/* ref_bundle
        ref_file_path="$(dirname $ref_path)"
        echo $ref_file_path

        # Construct Cell Ranger command
        cmd="cellranger count --id=${sample_id} --fastqs=$fastq_path --transcriptome=$ref_path --tissue=$spacename"
        
        # Add optional parameters to the command if specified
        if ${nosecondary}; then cmd+=" --nosecondary"; fi
        if [ -n "${localcores}" ]; then cmd+=" --localcores ${localcores}"; fi
        if [ -n "${localmem}" ]; then cmd+=" --localmem ${localmem}"; fi
        if [ -n "${localvmem}" ]; then cmd+=" --localvmem ${localvmem}"; fi
        echo $cmd

        # Run the Cell Ranger command
        /software/reboot-utils/cellranger/bin/cellranger $cmd

        # Package the results
        tar -czvf ${sample_id}_result.tar.gz "$ref_file_path/${sample_id}"
    }

    output {
        File out = '${sample_id}_result.tar.gz'
    }

    runtime {
        memory: memory_count + "G"
        cpu: cpu_count
        disk: storage_count + "GB"
        docker: "docker.io/man4ish/cellranger:latest"
    }
}

# Workflow for the entire spatial transcriptomics pipeline
workflow spatial_transcriptomics {
    input {
        File fastq_tar         # Tarball of gzipped FASTQ files
        File ref_tar           # Tarball of gzipped reference genome
        String sample_id       # Sample ID for Cell Ranger
        String spacename       # Name of the spatial tissue section
        Int memory_count       # Memory to allocate for the task
        Int cpu_count          # Number of CPUs to allocate for the task
        Int storage_count      # Disk space required for the task

        # Optional inputs
        Boolean? nosecondary   # Flag to skip secondary analysis
        Int? localcores        # Number of local cores to use
        String? localmem       # Local memory to allocate
        String? localvmem      # Local virtual memory to allocate
    }

    call cellranger_spatial {
        input: 
            fastq_tar=fastq_tar,
            ref_tar=ref_tar,
            sample_id=sample_id,
            spacename=spacename,
            memory_count=memory_count,
            cpu_count=cpu_count,
            storage_count=storage_count,
            nosecondary=nosecondary,
            localcores=localcores,
            localmem=localmem,
            localvmem=localvmem
    }
}

