FROM python:3.13-slim-bookworm
RUN mkdir -p /app
RUN mkdir -p /app/workbench /app/submission /app/output /data
COPY data/ /data/
COPY workbench/ /app/workbench/
ENV OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
WORKDIR /app
