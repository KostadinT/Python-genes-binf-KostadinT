import pickle
from os import path
from typing import Tuple, Generator, List

from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def load(organism_id: str) -> SeqRecord:
    """Load the NCBI record, use cached files if possible."""
    if not path.exists(path.join("data", f"{organism_id}.pkl.gz")):
        with Entrez.efetch(db="nucleotide", rettype="gb", id=organism_id) as handle:
            record = SeqIO.read(handle, "gb")
            with open(path.join("data", f"{organism_id}.pkl.gz"), "wb") as f:
                pickle.dump(record, f)
    else:
        with open(path.join("data", f"{organism_id}.pkl.gz"), "rb") as f:
            record = pickle.load(f)

    return record


def codons(seq: str) -> Generator[str, None, None]:
    """Walk along the string, three nucleotides at a time. Cut off excess."""
    for i in range(0, len(seq) - 2, 3):
        yield seq[i:i + 3]


def extract_gt_orfs(record, start_codons, stop_codons, validate_cds=True, verbose=False):
    """Extract the ground truth ORFs as indicated by the NCBI annotator in the
    gene coding regions (CDS regins) of the genome.

    Parameters
    ----------
    record: SeqRecord
    start_codons: List[str]
    stop_codons: List[str]
    validate_cds: bool
        Filter out NCBI provided ORFs that do not fit our ORF criteria.
    verbose: bool

    Returns
    -------
    List[Tuple[int, int, int]]
        tuples of form (strand, start_loc, stop_loc). Strand should be either 1
        for reference strand and -1 for reverse complement.

    """
    cds_regions = [f for f in record.features if f.type == "CDS"]

    orfs = []
    for region in cds_regions:
        loc = region.location
        seq = record.seq[loc.start.position:loc.end.position]
        if region.strand == -1:
            seq = seq.reverse_complement()
            
        if not validate_cds:
            orfs.append((region.strand, loc.start.position, loc.end.position))
            continue

        try:
            assert seq[:3] in start_codons, "Start codon not found!"
            assert seq[-3:] in stop_codons, "Stop codon not found!"
            # Make sure there are no stop codons in the middle of the sequence
            for codon in codons(seq[3:-3]):
                assert (
                    codon not in stop_codons
                ), f"Stop codon {codon} found in the middle of the sequence!"

            # The CDS looks fine, add it to the ORFs
            orfs.append((region.strand, loc.start.position, loc.end.position))

        except AssertionError as ex:
            if verbose:
                print(
                    "Skipped CDS at region [%d - %d] on strand %d"
                    % (loc.start.position, loc.end.position, region.strand)
                )
                print("\t", str(ex))

    return orfs


def find_orfs(sequence, start_codons, stop_codons):
    """Find possible ORF candidates in a single reading frame.

    Parameters
    ----------
    sequence: Seq
    start_codons: List[str]
    stop_codons: List[str]

    Returns
    -------
    List[Tuple[int, int]]
        tuples of form (start_loc, stop_loc)

    """
    # TODO
    i=0
    sr=-1;
    sp=-1;
    start=[];
    stop=[];
    total=[];
    #["ATG"], ["TAA", "TAG", "TGA"]
    while i <len(sequence)-2:
        k=sequence[i:i+3]
        #print("k= ",k," i= ",i," sr= ",sr," sp= ",sp," total= ",total)
        if(k in start_codons and sp==-1 and sr==-1): 
            sr=i
        elif(sr!=-1 and k in stop_codons): 
            sp=i
        if(sr!=-1 and sp!=-1):
            total.append((sr,sp+3))
            sr=-1
            sp=-1
        i+=3
    return total


def find_all_orfs(sequence, start_codons, stop_codons):
    """Find ALL the possible ORF candidates in the sequence using all six
    reading frames.

    Parameters
    ----------
    sequence: Seq
    start_codons: List[str]
    stop_codons: List[str]

    Returns
    -------
    List[Tuple[int, int, int]]
        tuples of form (strand, start_loc, stop_loc). Strand should be either 1
        for reference strand and -1 for reverse complement.

    """
    # TODO
    total=[];
    sequence1=str(sequence)
    srev=sequence1[::-1]
    srev=srev.replace("A","1")
    srev=srev.replace("T","A")
    srev=srev.replace("1","T")
    srev=srev.replace("C","1")
    srev=srev.replace("G","C")
    srev=srev.replace("1","G")
    totall=[]
    for i in range(6):
        if(i<3):
            kk=sequence[i:len(sequence)]
        else:
            kk=srev[(i%3):len(srev)]
        total=find_orfs(kk, start_codons, stop_codons)
        for j in range(len(total)):
            total[j]=list(total[j])
            for jj in range(len(total[j])):
                total[j][jj]+=(i%3)
        for j in total:
            tempL=[1]
            if(i>2):
                tempL=[-1]
            tempL.extend(j)
            totall.append(tuple(tempL))
    return totall
    

def translate_to_protein(seq):
    """Translate a nucleotide sequence into a protein sequence.

    Parameters
    ----------
    seq: str

    Returns
    -------
    str
        The translated protein sequence.

    """
    # TODO
    codonToProtein = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TGC':'C', 'TGT':'C', 'TGG':'W'
    }
    #  'TAA':'_', 'TAG':'_', 'TGA':'_',
    resProt=""
    i=0
    while i<len(seq)-2:
        k=seq[i:i+3]
        if k in codonToProtein.keys():
            resProt+=codonToProtein[k]
        i+=3
    return resProt
