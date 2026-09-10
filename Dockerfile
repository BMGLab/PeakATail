FROM python:3.10.12

WORKDIR /app 

COPY ./requirments.txt ./

RUN apt-get update && apt-get install -y bedtools


RUN  pip install --no-cache-dir -r requirments.txt

COPY . /app

ENV PYTHONPATH="/app"

RUN  pip install .

#ENTRYPOINT ["ema"]
CMD ["python" , "./ema/main.py"]
