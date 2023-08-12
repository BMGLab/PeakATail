cb_total = {}
cb_index = 0

def indexing(cb:str) -> int:
    global cb_index
    try:
        col = cb_total[cb]
        return col
    except KeyError:
        cb_index += 1
        cb_total[cb], col = cb_index, cb_index
        return col 