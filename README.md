# ProjectEMA
![image](https://github.com/user-attachments/assets/da4284d7-5aa8-4301-bef1-72c8c446ec20)

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
    git clone https://github.com/BMGLab/ProjectEMA.git
    ```
- go to project directory
    ```bash
    cd ProjectEMA
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
