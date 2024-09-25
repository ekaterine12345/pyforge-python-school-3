FROM continuumio/miniconda3

# Install necessary system dependencies for psycopg2 and other packages
RUN apt-get update && apt-get install -y libpq-dev gcc

# Create a conda environment file to install RDKit and other packages
COPY environment.yml /tmp/environment.yml

# Update conda to the latest version and configure solver options for faster environment solving
RUN conda update -n base -c defaults conda && \
    conda config --set channel_priority strict && \
    conda install mamba -n base -c conda-forge && \
    mamba env create -f /tmp/environment.yml

# Set the working directory
WORKDIR /app

# Copy the requirements.txt file first to leverage Docker cache
COPY requirements.txt /app/requirements.txt

# Install dependencies from requirements.txt
RUN pip install --no-cache-dir -r requirements.txt

# Copy the current directory into the container at /app
COPY . /app

# Command to run the application
CMD ["uvicorn", "src.main:app", "--host", "0.0.0.0", "--port", "8000"]
