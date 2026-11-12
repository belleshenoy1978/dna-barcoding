# Container for the web interface. Works both locally and on Hugging Face
# Spaces (Docker SDK), which runs the image as uid 1000 and expects the app on
# the port named in the Space README's `app_port`.
#
#   docker build -t barcode-pipeline .
#   docker run -p 8000:8000 barcode-pipeline                  # public defaults
#   docker run -p 8000:8000 -e BARCODE_PUBLIC=0 barcode-pipeline   # unrestricted
FROM python:3.12-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    KMP_DUPLICATE_LIB_OK=TRUE

# Containers are assumed to be reachable by other people: identification is off,
# uploads are capped and jobs expire. Override at run time for a private server.
ENV BARCODE_PUBLIC=1

RUN useradd -m -u 1000 user
USER user
ENV HOME=/home/user \
    PATH=/home/user/.local/bin:$PATH
WORKDIR /home/user/app

COPY --chown=user:user pyproject.toml README.md ./
COPY --chown=user:user barcode_pipeline ./barcode_pipeline
COPY --chown=user:user webapp ./webapp
COPY --chown=user:user examples ./examples
RUN pip install --no-cache-dir --user ".[web]"

EXPOSE 8000
CMD ["uvicorn", "webapp.app:app", "--host", "0.0.0.0", "--port", "8000"]
