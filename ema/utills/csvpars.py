import pandas as pd

def csv_to_list(file_path:str) -> list:
    data = pd.read_csv(file_path, header=None)
    return data.values.tolist()