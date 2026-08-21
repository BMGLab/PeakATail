from ema.countmatrix.polya import molecule_cb


class Peak():

    pasnumber = 0
    '''
    TODO
    '''

    @classmethod
    def reset_pasnumber(cls):
        """Reset the global PAS number counter to 0.

        Call between independent peak-calling runs (e.g. between datasets in
        multi-sample mode) so per-run pasnumbers start at 1. Safe to call;
        downstream code that keys on (dataset_id, pasnumber) tuples remains
        correct because dataset_id provides global uniqueness.
        """
        cls.pasnumber = 0

    def __init__(self, peak_list=None, peak_start=0, last_peak_end=0, peak_strand=True, cb_dict=None, cb_positions=None, polya_sites=None):
        '''
        :param peak_list: [[read_end1, height1], [read_end2, height]...., [read_endn, heightendn]]
        :param peak_start: this point is one read endpoint that detect as peak start but it not mean firat start of peak cause coulde merge multiple peaks
            peak_start assigne when signal turn True in peakcalling function
        :param peak_strand: presents gene strand
        :param cb_dict: collect all CellBarcodes are in peak
        :param polya_sites: poly(A) soft-clip evidence accumulated into this
            peak: {cleavage_site: [n_reads, umi_set, n_reads_without_umi,
            n_reads_f3844, umi_set_f3844, n_reads_without_umi_f3844]}.
            Populated via polya_counting() when --polya-evidence is on.
            Support is reported in DISTINCT MOLECULES; the f3844 slots hold
            the same counts restricted to reads passing samtools -F 3844
            (see ema.countmatrix.polya.clip_read_ok).
        '''
        self.peak_list = peak_list if peak_list is not None else []
        self.peak_start = peak_start
        self.last_peak_end = last_peak_end
        self.peak_strand = peak_strand
        self.cb_dict = cb_dict if cb_dict is not None else {}
        self.cb_positions = cb_positions if cb_positions is not None else {}
        self.polya_sites = polya_sites if polya_sites is not None else {}

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


    def cb_position_counting(self, end_pos:int, cb:str):
        '''Track which CB contributed a read at which position.'''
        if end_pos not in self.cb_positions:
            self.cb_positions[end_pos] = {}
        try:
            self.cb_positions[end_pos][cb] += 1
        except KeyError:
            self.cb_positions[end_pos][cb] = 1


    def polya_counting(self, site: int, cb: str, umi=None, primary: bool = True):
        '''Record one poly(A)-clipped read's inferred cleavage site.

        Mirrors cb_position_counting: called from the peak-calling loop for
        reads counted into this peak that carry a qualifying terminal
        poly(A) soft clip (see ema.countmatrix.polya.clip_site).

        :param site: 0-based cleavage coordinate from clip_site().
        :param cb: cell barcode (composite sample_cb string).
        :param umi: UMI (UB tag) or None. Distinct molecules per site are
            len(umi_set) + n_reads_without_umi.
        :param primary: clip_read_ok(read) -- False for secondary /
            supplementary / duplicate / qcfail alignments, which are
            tracked separately so both support units are available.
        '''
        rec = self.polya_sites.get(site)
        if rec is None:
            rec = [0, set(), 0, 0, set(), 0]
            self.polya_sites[site] = rec
        rec[0] += 1
        mol = None if umi is None else (molecule_cb(cb), umi)
        if mol is None:
            rec[2] += 1
        else:
            rec[1].add(mol)
        if primary:
            rec[3] += 1
            if mol is None:
                rec[5] += 1
            else:
                rec[4].add(mol)

    def polya_support(self, pas_1: int, pas_2: int, strand: bool, window: int):
        '''Clip support for one emitted PAS: reads and distinct molecules
        whose cleavage site lies within +/-window of the PAS's strand-aware
        3' base (forward: bed_end - 1; reverse: bed_start).

        Returns (n_clip_reads, n_distinct_molecules, n_clip_reads_f3844,
        n_distinct_molecules_f3844).  BED column 5 carries the MOLECULE
        count; the sidecar carries all four.
        '''
        if not self.polya_sites:
            return 0, 0, 0, 0
        bed_start = min(pas_1, pas_2)
        bed_end = max(pas_1, pas_2)
        three = bed_start if strand else bed_end - 1
        nreads = 0
        n_no_umi = 0
        umis = set()
        nreads_f = 0
        n_no_umi_f = 0
        umis_f = set()
        for site, rec in self.polya_sites.items():
            if abs(site - three) <= window:
                nreads += rec[0]
                umis |= rec[1]
                n_no_umi += rec[2]
                nreads_f += rec[3]
                umis_f |= rec[4]
                n_no_umi_f += rec[5]
        return (nreads, len(umis) + n_no_umi,
                nreads_f, len(umis_f) + n_no_umi_f)

    def pasfind(self) -> int:
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
                return 0, 0
            
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
                        return 0, 0
                    
                else:
                    for item in self.peak_list:
                        if item[1] >= ther and pas_1 == 0:
                            pas_1, pas_cov=  item[0], item[1]
                        elif item[0] >= pas_1 and item[1] > pas_cov:
                            pas_2 = item[0]
                            return pas_1, pas_2
                    else:
                        return 0, 0

        except Exception:
            return 0, 0
        