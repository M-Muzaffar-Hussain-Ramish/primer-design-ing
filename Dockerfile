FROM python:3.12-slim

# Install system-level MSA tools and dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
        clustalo \
        && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

# Create default cache and results directories
RUN mkdir -p cache results

EXPOSE 8000

CMD ["uvicorn", "api.main:app", "--host", "0.0.0.0", "--port", "8000"]
