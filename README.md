# PeakATail
<img src="https://github.com/user-attachments/assets/9fcae81d-59f9-4bd4-aa22-1d09a03bf2b1" width="400" height="300">


### Base installation 
- you should have samtools and bedtools installed

    on your system
    ```bash 
    sudo apt-get install bedtools
    ```
    ```bash
    sudo apt-get install samtools
    ```

### install

- clone repo 
    ```bash
    git clone https://github.com/BMGLab/PeakATail.git
    ```
- go to project directory
    ```bash
    cd PeakATail
    ```
- open python env
    ```bash
    python3.10 -m venv .venv
    ```
[Venv](https://docs.python.org/3/library/venv.html)
- activate environment 
    ```bash 
    source ,venv/bin/activate
    ```
- install requierments
    ```bash
    pip install -r requirments.txt
    ```
- install ema 
- use -e flag to could able automatic updates
    ```bash
    pip install -e .
    ```

### Run 

- ema merge
- use to merge multiple bam files 
    ```bash
    ema_merge --bamFiles bamfiles.yaml --threads 4
    ```

- [Example yaml File](example.yaml)

- ema 
- this command peakcall and cluster base on PAS 
    ```bash
    ema --bamDir path/to/merged.bam \
        --sequenceLen 98 \
        --CellBarcodeLen 16 \
        --BarcodeTag CB \
        --gtfDir path/to/human.gtf

    ```
- or you could have all in a single line 
