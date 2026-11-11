FROM python:3.13-slim-bookworm
RUN mkdir -p /app
RUN mkdir -p /app/submission
RUN mkdir -p /app/output
RUN pip install --no-cache-dir pytest==8.4.1
COPY . /tests
RUN chown -R root:root /tests && chmod -R go-rwx /tests && mkdir -p /logs/verifier && chmod 700 /logs/verifier
ENV OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
WORKDIR /app
