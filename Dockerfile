FROM continuumio/miniconda3

# Install necessary system dependencies for psycopg2 and other packages
RUN apt-get update && apt-get install -y libpq-dev gcc

RUN conda install conda-forge::rdkit

#RUN pip install fastapi uvicorn
#RUN pip install pytest httpx pytest-asyncio


# Set the working directory
WORKDIR /app

# Copy the requirements.txt file first to leverage Docker cache
COPY requirements.txt /app/requirements.txt

# Install dependencies from requirements.txt
RUN pip install --no-cache-dir -r requirements.txt

# Copy the current directory into the container at /app
#COPY ./src /app/src
COPY . /app

# Command to run the application
CMD ["uvicorn", "src.main:app", "--host", "0.0.0.0", "--port", "8000"]