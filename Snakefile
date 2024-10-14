import os
import subprocess
import pathlib
from time import localtime
from snakemake.utils import min_version

PROJECT = "hybrid-sentry-434815-v3"
MACHINE_TYPE = "c2d-highcpu-56"
IMAGE = "ubuntu-2204-lts"
BOOT_DISK_SIZE = "10GB"

min_version("8.17.0")

rule all:
    input: "instance_deleted.txt"

rule create_instance:
    params:
        INSTANCE_NAME=config.get("instance-name", f"vaxpress-instance-{localtime().tm_hour}-{localtime().tm_min}-{localtime().tm_sec}"),
        ZONE=config.get("zone", "asia-south1-a")
    output: "instance_created.txt"
    shell:
        """
        gcloud compute --project={PROJECT} instances create {params.INSTANCE_NAME} \
            --zone={params.ZONE} --machine-type={MACHINE_TYPE} \
            --image-family={IMAGE} --image-project=ubuntu-os-cloud \
            --boot-disk-size={BOOT_DISK_SIZE} \
            && touch {output}
        """

rule add_ssh_key_to_instance:
    input: "instance_created.txt"
    params:
        INSTANCE_NAME=config.get("instance-name", f"vaxpress-instance-{localtime().tm_hour}-{localtime().tm_min}-{localtime().tm_sec}"),
        ZONE=config.get("zone", "asia-south1-a")
    output: "ssh_key_added.txt"
    shell:
        """  
        gcloud compute instances add-metadata {params.INSTANCE_NAME} \
            --zone={params.ZONE} --metadata "ssh-keys=$(whoami):$(cat ~/.ssh/id_rsa.pub)"

        maxtries=12
        count=0

        until gcloud compute ssh {params.INSTANCE_NAME} --zone={params.ZONE} --command "echo 'SSH key applied'" >/dev/null 2>&1 || [ $count -ge $maxtries ]; do
            sleep 3
            count=$((count+1))
        done

        if [ $count -ge $maxtries ]; then
            echo "SSH key propagation failed after $maxtries attempts." && exit 1;
        fi

        touch {output}
        """

rule broadcast_data:
    input: "ssh_key_added.txt"
    params: 
        seq=config['input'],
        config=config['vaxpress-config'],
        SCRIPT = workflow.source_path("./gcloud-setup.sh"),
        INSTANCE_NAME=config.get("instance-name", f"vaxpress-instance-{localtime().tm_hour}-{localtime().tm_min}-{localtime().tm_sec}"),
        ZONE=config.get("zone", "asia-south1-a")
    output: "copied_data.txt"
    shell:
        """
        gcloud compute scp --strict-host-key-checking=no {params.seq} {params.INSTANCE_NAME}:sequence.fasta\
        --zone={params.ZONE} \
        && gcloud compute scp --strict-host-key-checking=no {params.config} {params.INSTANCE_NAME}:config.json \
            --zone={params.ZONE} \
        && gcloud compute scp --strict-host-key-checking=no {params.SCRIPT} {params.INSTANCE_NAME}:~ \
            --zone={params.ZONE} \
        && touch {output}
        """


rule run_shell_script:
    input: "copied_data.txt"
    params:
        INSTANCE_NAME=config.get("instance-name", f"vaxpress-instance-{localtime().tm_hour}-{localtime().tm_min}-{localtime().tm_sec}"),
        ZONE=config.get("zone", "asia-south1-a")
    output: "ran_script.txt"
    shell:
        """
        gcloud compute ssh --strict-host-key-checking=no {params.INSTANCE_NAME} --zone={params.ZONE} --command \
            'bash gcloud-setup.sh' \
            && touch {output}
        """

rule run_vaxpress:
    input: "ran_script.txt"
    params:
        INSTANCE_NAME=config.get("instance-name", f"vaxpress-instance-{localtime().tm_hour}-{localtime().tm_min}-{localtime().tm_sec}"),
        ZONE=config.get("zone", "asia-south1-a")
    output: "simulation_complete.txt"
    shell:
        """
        gcloud compute ssh --strict-host-key-checking=no {params.INSTANCE_NAME} --zone={params.ZONE} --command \
            'source VaxPress/.venv/bin/activate && vaxpress -p $(nproc) -i sequence.fasta -o results --preset config.json' \
            && touch {output}
        """

rule download_output:
    input: "simulation_complete.txt"
    params:
        output=config['output-dir'],
        INSTANCE_NAME=config.get("instance-name", f"vaxpress-instance-{localtime().tm_hour}-{localtime().tm_min}-{localtime().tm_sec}"),
        ZONE=config.get("zone", "asia-south1-a")
    output: "output_downloaded.txt"
    shell:
        """
        gcloud compute scp --strict-host-key-checking=no --recurse {params.INSTANCE_NAME}:results {params.output} --zone={params.ZONE}\
            && touch {output}
        """

rule delete_instance:
    input: "output_downloaded.txt"
    params:
        INSTANCE_NAME=config.get("instance-name", f"vaxpress-instance-{localtime().tm_hour}-{localtime().tm_min}-{localtime().tm_sec}"),
        ZONE=config.get("zone", "asia-south1-a")
    output: "instance_deleted.txt"
    shell:
        """
        gcloud compute instances delete {params.INSTANCE_NAME} \
            --zone={params.ZONE} --quiet \
            && touch {output}
        """


