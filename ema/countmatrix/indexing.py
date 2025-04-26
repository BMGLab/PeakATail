cb_total = {}
cb_index = 0

"""def indexing(cb:str) -> int:
    index_map = {cb: i+1 for i, cb in enumerate(dict.fromkeys(cb_list))}

    global cb_index
    global cb_total
    try:
        col = cb_total[cb]
        return col
    except KeyError:
        cb_index += 1
        cb_total[cb], col = cb_index, cb_index
        return col"""

def indexing(cb:list) -> int:
    cb_num = len(cb_total)
    index_map = {cb: cb_num+i+1 for i, cb in enumerate(dict.fromkeys(cb_list))}
    cb_total.update(index_map)