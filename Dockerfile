# Dockerfile
FROM python:3.11-slim

WORKDIR /pipeline

# Copy the pipeline script into the image
COPY iss_pipeline.py .
COPY header_editor.py .

# Install inside the image
RUN pip install --no-cache-dir insilicoseq
RUN pip install --no-cache-dir biopython

# Default command (we will override with args in Singularity job)
CMD ["python", "/pipeline/iss_pipeline.py"]