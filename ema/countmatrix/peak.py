class Peak():
    # use to merge peaks 
    last_peak_end  = None
    '''
    TODO
    '''

    def __init__(self, peak_list:list, peak_start:int, peak_strand:bool, cb_dict:dict):
        '''
        :param peak_list: [[read_end1, height1], [read_end2, height]...., [read_endn, heightendn]]
        :param peak_start: this point is one read endpoint that detect as peak start 
        :param peak_strand: presents gene strand 
        :param cb_dict: collect all CellBarcodes are in peak
        '''
        self.peak_list = peak_list
        self.peal_start = peak_start
        self.peak_strand = peak_strand
        self.cb_dict = cb_dict


    # each time data_array slicing ubdate peak_add 
    def peak_add(self, data_array:list, slice_loc:int):
        array_len = len(data_array)

        for i in range(slice_loc):
            self.peak_list.append([data_array[i], array_len-i])


    def cb_counting(self, cb:str):
        try:
            self.cb_dict[cb] += 1
        except KeyError:
            self.cb_dict[cb] = 1


    
    def pasfind(self):
        '''
        pasfind method loop on peak_list and fidn max_height
        find 5% then find value larger than max_heght 5% and assigne it as pas 

        if method return False mean peak is not valid 

        TODO will develope an static model
        '''

        try:
            max_height = max(self.peak_list, key=lambda x:x[1])[1]
            pas_1, pas_2, pas_cov = 0, 0, 0
            # 5% must be atleast 1 so max height must be at least 20 
            
            if max_height <= 20:
                return False
            
            else:
                # diffirent two block for different strands
                #TODO code repitition 
                ther = max_height*0.05

                if self.peak_strand == False:
                    for item in reversed(self.peak_list[:-1]):
                        if item[1] >= ther and pas_1 == 0:
                            pas_1, pas_cov=  item[0], item[1]
                        elif item[0] <= pas_1 and item[1] > pas_cov:
                            pas_2 = item[0]
                            return pas_1, pas_2
                    
                    else:
                        return False
                    
                else:
                    for item in self.peak_list:
                        if item[1] >= ther and pas_1 == 0:
                            pas_1, pas_cov=  item[0], item[1]
                        elif item[0] >= pas_1 and item[1] > pas_cov:
                            pas_2 = item[0]
                            return pas_1, pas_2
                    else:
                        return False

        except:
            return False
        