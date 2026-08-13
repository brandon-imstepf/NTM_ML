# Dockerfile, Image, Container

FROM python:3.12.7

WORKDIR /app

ADD nn-hyp-tuning.py /app/
ADD my_nn.py /app/
ADD nn_hyp_utils.py /app/
ADD training_data /app/training_data/
ADD run_script.sh /app/

COPY requirements.txt /app/requirements.txt

RUN pip install --no-cache-dir -r /app/requirements.txt
RUN chmod +x ./run_script.sh

ENTRYPOINT [ "/app/run_script.sh" ]
