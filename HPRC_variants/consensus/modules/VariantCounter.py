import sys
# Define the Variant counter
class VariantCounter:
    # for each site(0, 1, ..., ref_len-1), count the number of variants
    # variant type: SNV, INS, DEL
    # and the read depth: DP
    def __init__(self, ref_len):
        self.ref_len = ref_len
        self.variant_count = [0] * ref_len
        self.variant_count_SNV = [0] * ref_len
        self.variant_count_INS = [0] * ref_len
        self.variant_count_DEL = [0] * ref_len
        self.DP = [0] * (ref_len + 1) # read depth 
        
        # variant frequency
        self.SNV_freq = [0] * ref_len
        self.INS_freq = [0] * ref_len
        self.DEL_freq = [0] * ref_len
    def add_variant(self, variant_type, start, end):
        # variant_type: SNV, INS, DEL
        # start, end: 0-based
        # if variant_type == "SNV":
        if variant_type == "*":
            for i in range(start, end+1):
                self.variant_count_SNV[i] += 1
        # elif variant_type == "INS":
        elif variant_type[0] == "+":
            # start should be equal to end
            for i in range(start, end+1):
                self.variant_count_INS[i] += 1
        # elif variant_type == "DEL":
        elif variant_type == "-":
            for i in range(start, end+1):
                self.variant_count_DEL[i] += 1
        elif variant_type == ":": # match
            pass
        else:
            print("Error: variant_type should be SNV, INS or DEL")
            sys.exit(1)
        for i in range(start, end+1):
            self.variant_count[i] += 1

    def add_DP(self, start, end):
        # start, end: 0-based
        for i in range(start, end+1):
        # for i in range(start, max(end+1, self.ref_len)):
            # ここでずるいことをしている。PAF では len = 13332 なのに
            # ref_start = 0 and ref_end = 13332 というマッピングがある。
            # そのリードの長さは 13333 になるはず。どういうこと？
            self.DP[i] += 1
    def print(self):
        print(self.variant_count)
        print(self.variant_count_SNV)
        print(self.variant_count_INS)
        print(self.variant_count_DEL)
        print(self.DP)

    def calc_freq(self):
        self.SNV_freq = [x / y for x, y in zip(self.variant_count_SNV, self.DP)]
        self.INS_freq = [x / y for x, y in zip(self.variant_count_INS, self.DP)]
        self.DEL_freq = [x / y for x, y in zip(self.variant_count_DEL, self.DP)]
        # print(self.SNV_freq)
    def to_csv(self, output_csv, samplename):
        # columns: POS(1-based), SNV_COUNT, INS_COUNT, DEL_COUNT, DP, SNV_FREQ, INS_FREQ, DEL_FREQ
        self.calc_freq()
        with open(output_csv, "w") as f:
            f.write(f"POS(1-based),SNV_COUNT,INS_COUNT,DEL_COUNT,DP,SNV_FREQ,INS_FREQ,DEL_FREQ\n")
            for i in range(self.ref_len):
                # print(f"{i+1},{self.variant_count_SNV[i]},{self.variant_count_INS[i]},{self.variant_count_DEL[i]},{self.DP[i]},{self.SNV_freq[i]},{self.INS_freq[i]},{self.DEL_freq[i]}")
                f.write(f"{i},{self.variant_count_SNV[i]},{self.variant_count_INS[i]},{self.variant_count_DEL[i]},{self.DP[i]},{self.SNV_freq[i]},{self.INS_freq[i]},{self.DEL_freq[i]}\n")

def load_cstag(parsed_cstag, ref_start, variant_counter):
    position = ref_start
    for cstag_unit in parsed_cstag:
        cstag_type = cstag_unit[0]
        # print(f"position: {position}, cstag_unit: {cstag_unit}")
        
        if cstag_type == ":":
            # match
            match_len = int(cstag_unit[1:])
            position += match_len
        elif cstag_type == "-":
            # deletion
            del_seq = cstag_unit[1:]
            del_len = len(del_seq)
            variant_counter.add_variant("-", position, position + del_len - 1)
            position += del_len
        elif cstag_type == "+":
            # insertion
            ins_seq = cstag_unit[1:]
            ins_len = len(ins_seq)
            variant_counter.add_variant("+", position, position)
            # position += 1
        elif cstag_type == "*":
            # SNV
            snv_unit = cstag_unit[1:]
            # snv_before = snv_unit[0]
            # snv_after = snv_unit[1]
            variant_counter.add_variant("*", position, position)
            position += 1
        else:
            print("Error: cstag_type should be :+-*")
            sys.exit(1)

    return variant_counter