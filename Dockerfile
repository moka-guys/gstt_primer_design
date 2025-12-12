FROM python:3.11-slim

# Set work directory
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    libpq-dev \
    bcftools \
    bowtie2 \
    && rm -rf /var/lib/apt/lists/*

# Install Python packages
COPY requirements.txt .
RUN pip install --upgrade pip setuptools wheel \
    && pip install --no-cache-dir --force-reinstall -r requirements.txt

# define env

ENV DOWNLOAD_FOLDER /app/output
ENV SESSION_FILE_DIR /var/local/gstt_primer_design/flask_session

# Copy the rest of the app except .dockerignore
COPY . .


# Expose flask port
EXPOSE 5000

# Run the Flask app with gunicorn
CMD ["gunicorn", "--bind", "0.0.0.0:5000", "--workers", "2", "--timeout", "120", "--keep-alive", "5", "wsgi:app"]

