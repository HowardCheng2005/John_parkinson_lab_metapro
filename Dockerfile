# Dockerfile
FROM python:3.11-slim

WORKDIR /pipeline

# Copy only the pipeline script into the image
COPY iss_pipeline.py .

# Install insilicoseq inside the image
RUN pip install --no-cache-dir insilicoseq

# Default command (we will override with args in Singularity job)
CMD ["python", "/pipeline/iss_pipeline.py"]