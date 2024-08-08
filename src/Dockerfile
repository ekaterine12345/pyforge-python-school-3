FROM continuumio/miniconda3
RUN conda install conda-forge::rdkit

RUN pip install fastapi uvicorn
RUN pip install pytest httpx pytest-asyncio


# Set the working directory
WORKDIR /app

# Copy the current directory contents into the container at /app
COPY . /app

# Command to run the application
CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]