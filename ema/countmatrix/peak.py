class Peak():
    # use to merge peaks 
    last_peak_end  = None
    '''
    TODO
    '''

    def __init__(self, peak_list:list, peak_start:int, peak_strand:bool):
        self.peak_list = peak_list
        self.peal_start = peak_start
        self.peak_strand = peak_strand

    # each time data_array slicing ubdate peak_add 
    def peak_add(self, data_array:list, slice_loc:int):
        array_len = len(data_array)

        for i in range(slice_loc):
            self.peak_list.append([data_array[i], array_len-i])

    def pasfind():
        