docker build -t us.gcr.io/alignments-65005/sentieon .

docker run -v /mnt/disks/:/mnt/data --env-file env.list -it --rm us.gcr.io/alignments-65005/sentieon /bin/bash -c "source mito-workflow.sh"

