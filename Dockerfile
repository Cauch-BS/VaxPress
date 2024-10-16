# Use the official Python image as the base
FROM python:3.11-slim

# Set the working directory
WORKDIR /app

# Update and install dependencies
RUN apt-get update && apt-get upgrade -y && \
    apt-get install -y python3-pip python3-venv build-essential && \
    rm -rf /var/lib/apt/lists/*

COPY . /app/VaxPress

WORKDIR /app/VaxPress

RUN python3 -m venv .venv && \
    .venv/bin/pip install . && \
    .venv/bin/pip install linearpartition-unofficial linearfold-unofficial

ENTRYPOINT ["/app/VaxPress/.venv/bin/vaxpress"]

#docker build -t gcr.io/your-project-id/vaxpress:0.9.1 .
#gcloud auth configure-docker
#docker push gcr.io/your-project-id/vaxpress:0.9.1
#gcloud compute ssh vaxpress-vm --zone=us-central1-a --command "sudo docker run --rm gcr.io/your-project-id/vaxpress:latest -h"