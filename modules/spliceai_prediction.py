import numpy as np
from modules.spliceai.utils import Annotator
from modules.spliceai.utils import get_delta_scores
from dataclasses import dataclass

def getSpliceAI(chr, pos, ref, alts, buildver, ref_fasta_file):
    @dataclass
    class Record:
        chrom: str
        pos: int
        ref: str
        alts: list

    # Create an instance
    record = Record(chr, pos, ref, alts)
    dist_var = 50  #  Maximum distance between the variant and gained/lost splice site (SpliceAI Default)
    mask = np.array([True] * len(record.alts))  

    # Create an instance of the Annotator class
    ref_fasta = ref_fasta_file
    annotations = 'grch37' if buildver == 'hg19' else 'grch38'
    ann = Annotator(ref_fasta, annotations)
    DS_AG = -1
    DS_AL = -1
    DS_DG = -1
    DS_DL = -1
    delta_scores = get_delta_scores(record, ann, dist_var, mask)
    if delta_scores:
        DS_AG = float(delta_scores[0].split('|')[2])
        DS_AL = float(delta_scores[0].split('|')[3])
        DS_DG = float(delta_scores[0].split('|')[4])
        DS_DL = float(delta_scores[0].split('|')[5])

    return DS_AG, DS_AL, DS_DG, DS_DL
