version 1.0

# Define the task to run spaceranger count
task spaceranger_count {
  input {
    String spaceranger_path
    String reference
    String fastqs
    String output_dir
    String sample_id
    String tissue_image
    String white_list
    String additional_flags
  }

  command {
    ${spaceranger_path}/bin/spaceranger count \
      --id=${sample_id} \
      --transcriptome=${reference} \
      --fastqs=${fastqs} \
      --output-dir=${output_dir} \
      --image=${tissue_image} \
      --whitelist=${white_list} \
      ${additional_flags}
  }

  output {
    String output_dir = "${output_dir}"
  }

  runtime {
    docker: "spaceranger:3.0.1"  # Docker image that includes spaceranger
  }

  # Optionally define resources, such as CPU and memory
  # runtime {
  #   memory: "8 GB"
  #   cpu: 4
  # }
}

# Define the workflow that runs the spaceranger count task
workflow spaceranger_workflow {
  input {
    String spaceranger_path
    String reference
    String fastqs
    String output_dir
    String sample_id
    String tissue_image
    String white_list
    String additional_flags
  }

  call spaceranger_count {
    input:
      spaceranger_path = spaceranger_path,
      reference = reference,
      fastqs = fastqs,
      output_dir = output_dir,
      sample_id = sample_id,
      tissue_image = tissue_image,
      white_list = white_list,
      additional_flags = additional_flags
  }

  output {
    String result_dir = spaceranger_count.output_dir
  }
}

