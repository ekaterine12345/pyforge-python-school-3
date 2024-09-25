FROM continuumio/miniconda3

# Install necessary system dependencies for psycopg2 and other packages
RUN apt-get update && apt-get install -y libpq-dev gcc

# Install mamba and update conda in separate steps to avoid memory overload
RUN conda update -n base -c defaults conda && \
    conda install mamba -n base -c conda-forge

# Create the environment using mamba in a separate step
COPY environment.yml /tmp/environment.yml
RUN mamba env create -f /tmp/environment.yml

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
