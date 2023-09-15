FROM python:3.10.12

WORKDIR /app 

COPY ./requirements.txt ./
RUN  pip install --no--cach-dir -r requirements.txt

COPY ./ema /app/ema

ENTRYPOINT ["ema"]

CMD ["ema", "--bamDir", "BAM_DIR", "--sequenceLen", "SEQLEN", "--CellBarcodeLen", "CB_LEN", "--BarcodeTag", "BARCODE_TAG"]
