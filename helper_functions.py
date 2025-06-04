import math

def jukes_cantor(reference_sequence: str, distant_sequence: str) -> float:
    """The Jukes-Cantor correction for estimating genetic distances
    calculated with Hamming distance

    Parameters
    ----------
    referene_sequence: str
        A string of nucleotides in a sequence used as a reference
        in an alignment with other (i.e. ATGGC-TAG)
    distant_sequence: str
        A string of nucleotides in a sequence after the alignment
        with a reference (i.e. ATG-CTTAG)

    Returns
    -------
    float
        The Jukes-Cantor corrected genetic distance using Hamming distance.

    """
    refSeqLen = len(reference_sequence)
    distSeqLen = len(distant_sequence)
    mutationCounter = 0
    extra = 0
    for i in range(min(refSeqLen,distSeqLen)):
        refEl = reference_sequence[i]
        distEl = distant_sequence[i]
        if (refEl == "-" or distEl == "-"):
            extra += 1
        else:
            if refEl != distEl:
                mutationCounter += 1
    p = mutationCounter / (min(refSeqLen,distSeqLen) - extra)
    correctedDistance = (-3.0/4.0) * math.log(1 - ((4.0/3.0) * p))
    return correctedDistance


def kimura_two_parameter(reference_sequence: str, distant_sequence: str) -> float:
    """The Kimura Two Parameter correction for estimating genetic distances
    calculated with Hamming distance

    Parameters
    ----------
    referene_sequence: str
        A string of nucleotides in a sequence used as a reference
        in an alignment with other (i.e. ATGGC-TAG)
    distant_sequence: str
        A string of nucleotides in a sequence after the alignment
        with a reference (i.e. ATG-CTTAG)

    Returns
    -------
    float
        The Kimura corrected genetic distance using Hamming distance.

    """
    refSeqLen = len(reference_sequence)
    distSeqLen = len(distant_sequence)
    mutationCounterTransitions = 0
    mutationCounterTransversions = 0
    extra = 0
    transitions = [("A","G"), ("G","A"),("T","C"),("C","T")]
    for i in range(min(refSeqLen,distSeqLen)):
        refEl = reference_sequence[i]
        distEl = distant_sequence[i]
        if (refEl == "-" or distEl == "-"):
            extra += 1
        else:
            if refEl != distEl:
                if (refEl,distEl) in transitions:
                    mutationCounterTransitions += 1
                else:
                    mutationCounterTransversions += 1

    p = mutationCounterTransitions / (min(refSeqLen,distSeqLen) - extra)
    q = mutationCounterTransversions / (min(refSeqLen,distSeqLen) - extra)

    correctedDistance = (-1/2) * math.log((1 - 2 * p - q) * math.sqrt(1 - 2 * q))
    return correctedDistance