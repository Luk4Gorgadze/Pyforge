# Use the Miniconda3 image
FROM continuumio/miniconda3

# Set the working directory
WORKDIR /app

# Copy the requirements file
COPY requirements.txt .

# Install RDKit and other dependencies
RUN conda install -c conda-forge -y rdkit && \
    pip install --no-cache-dir -r requirements.txt

# Copy the application code
COPY . /app


ENV PYTHONPATH=/app/src

RUN mkdir -p /src/data

# Expose the port the app runs on
EXPOSE 8000

# Run the FastAPI application
CMD ["uvicorn", "main.main:app", "--host", "0.0.0.0", "--port", "8000", "--reload"]