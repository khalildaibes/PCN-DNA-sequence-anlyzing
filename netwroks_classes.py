class NucleotidesNetNodes:
    def __init__(self, sub_seq=None, sequence=None, subject=None, seq_id=None, position=None):
        self.sub_seq = sub_seq
        self.sequence = sequence
        self.subject = subject
        self.seq_id = seq_id
        self.position = position


# A CLASS TO HOLD THE AMINO ACID NODES DATA AND METADATA
class AminoAcidNetNodes:
    def __init__(self, sequence=None, subject=None):
        self.sequence = sequence
        self.subject = subject