cb_total = {}
cb_index = 0

def indexing(cb:str) -> int:
    global cb_index
    global cb_total
    idx = cb_total.get(cb)
    if idx is not None:
        return idx
    else:
        cb_index += 1
        cb_total[cb], col = cb_index, cb_index
        return col