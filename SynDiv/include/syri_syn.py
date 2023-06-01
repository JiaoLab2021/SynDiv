#!/home/ltan/anaconda3/envs/syri_env/bin/python3
# Created on 2023/4/13
# @author: Zezhen Du
# email: dzz0539@gmail.com or dzz0539@163.com


import os
import sys
import logging
import numpy as np
from collections import deque, defaultdict
import pandas as pd
import psutil
from gc import collect


class alignmentBlock:
    def __init__(self, id, children, data):
        self.id = id
        self.children = list(children)
        self.parents = []
        self.score = (data.aLen + data.bLen) * data.iden
        self.bestParentID = -1
    
    def addParent(self, parentID):
        self.parents.append(parentID)
        
    def bestParent(self,parentID,maxScore):
        self.bestParentID = parentID
        self.score = self.score + maxScore


def revcomp(seq):
    assert type(seq) == str
    old = 'ACGTRYKMBDHVacgtrykmbdhv'
    rev = 'TGCAYRMKVHDBtgcayrmkvhdb'
    tab = str.maketrans(old, rev)
    return seq.translate(tab)[::-1]


def readfasta(f):
    from gzip import open as gzopen
    from gzip import BadGzipFile
    from collections import deque
    import sys
    out = {}
    chrid = ''
    chrseq = deque()
    # Test if the file is Gzipped or not
    with gzopen(f, 'rb') as fin:
        try:
            fin.read(1)
            isgzip = True
        except BadGzipFile:
            isgzip = False
    try:
        if isgzip:
            with gzopen(f, 'rb') as fin:
                for line in fin:
                    if b'>' in line:
                        if chrid != '':
                            out[chrid] = ''.join(chrseq)
                            chrid = line.strip().split(b'>')[1].split(b' ')[0].decode()
                            chrseq = deque()
                        else:
                            chrid = line.strip().split(b'>')[1].split(b' ')[0].decode()
                        if chrid in out.keys():
                            # sys.exit(" Duplicate chromosome IDs are not accepted. Chromosome ID {} is duplicated. Provided chromosome with unique IDs".format(chrid))
                            return out
                    else:
                        chrseq.append(line.strip().decode())
        else:
            with open(f, 'r') as fin:
                for line in fin:
                    if '>' in line:
                        if chrid != '':
                            out[chrid] = ''.join(chrseq)
                            chrid = line.strip().split('>')[1].split(' ')[0]
                            chrseq = deque()
                        else:
                            chrid = line.strip().split('>')[1].split(' ')[0]
                        if chrid in out.keys():
                            # sys.exit(" Duplicate chromosome IDs are not accepted. Chromosome ID {} is duplicated. Provided chromosome with unique IDs".format(chrid))
                            return out
                    else:
                        chrseq.append(line.strip())
    except Exception as e:
        raise Exception(e)
    if chrid != '':
        out[chrid] = ''.join(chrseq)
    # TODO: add check for the validation of input fasta files
    return out
# END

def cgtpl(cg):
    """
    Takes a cigar string as input and returns a cigar tuple
    """
    for i in "MIDNSHPX=":
        cg = cg.replace(i, ';'+i+',')
    return [i.split(';') for i in cg.split(',')[:-1]]
#end


def samtocoords(f):
    from pandas import DataFrame
    from collections import deque
    logger = logging.getLogger('SAM reader')
    rc = {}        # Referece chromosomes
    rcs = {}        # Selected chromosomes
    al = deque()    # Individual alignment
    try:
        with open(f, 'r') as fin:
            for l in fin:
                if l[:3] == '@SQ':
                    c, s = 0, 0
                    for h in l.strip().split()[1:]:
                        h = h.split(':')
                        if h[0] == 'SN': c = h[1]
                        if h[0] == 'LN': s = int(h[1])
                    rcs[c] = s
                    continue
                elif l[0] == '@': continue

                l = l.split('\t')[:6]
                # if l[1] == '2064': break
                if l[2] == '*':
                    # logger.warning(l[0]+ ' do not align with any reference sequence and cannot be analysed. Remove all unplaced scaffolds and contigs from the assemblies.')  # Skip rows corresponding to non-mapping sequences (contigs/scaffolds)
                    continue

                if 'M' in l[5]:
                    logger.error('Incorrect CIGAR string found. CIGAR string can only have I/D/H/S/X/=. CIGAR STRING: ' + l[5])
                    # sys.exit()
                    return pd.DataFrame()
                cgt = [[int(j[0]), j[1]] for j in [i.split(';') for i in l[5].replace('S', ';S,').replace('H', ';H,').replace('=', ';=,').replace('X', ';X,').replace('I', ';I,').replace('D', ';D,').split(',')[:-1]]]
                if len(cgt) > 2:
                    if True in [True if i[1] in ['S', 'H'] else False for i in cgt[1:-1]]:
                        logger.error("Incorrect CIGAR string found. Clipped bases inside alignment. H/S can only be in the terminal. CIGAR STRING: " + aln.cigarstring)
                        # sys.exit()
                        return pd.DataFrame()

                bf = '{:012b}'.format(int(l[1]))

                rs = int(l[3])
                re = rs - 1 + sum([i[0] for i in cgt if i[1] in ['X', '=', 'D']])

                if bf[7] == '0':    # forward alignment
                    if cgt[0][1] == '=':
                        qs = 1
                    elif cgt[0][1] in ['S', 'H']:
                        qs = cgt[0][0] + 1
                    else:
                        print('ERROR: CIGAR string starting with non-matching base')
                    qe = qs - 1 + sum([i[0] for i in cgt if i[1] in ['X', '=', 'I']])
                elif bf[7] == '1':  # inverted alignment
                    if cgt[-1][1] == '=':
                        qs = 1
                    elif cgt[-1][1] in ['S', 'H']:
                        qs = cgt[-1][0] + 1
                    else:
                        print('ERROR: CIGAR string starting with non-matching base')
                    qe = qs - 1 + sum([i[0] for i in cgt if i[1] in ['X', '=', 'I']])
                    qs, qe = qe, qs

                al.append([
                    rs,
                    re,
                    qs,
                    qe,
                    abs(re-rs) + 1,
                    abs(qs-qe) + 1,
                    format((sum([i[0] for i in cgt if i[1] == '=']) / sum(
                        [i[0] for i in cgt if i[1] in ['=', 'X', 'I', 'D']])) * 100, '.2f'),
                    1,
                    1 if bf[7] == '0' else -1,
                    l[2],
                    l[0],
                    "".join([str(i[0])+i[1] for i in cgt if i[1] in ['=', 'X', 'I', 'D']])
                ])
                rcs[l[2]] = 1
            rcs = list(rcs.keys())
            for k in list(rc.keys()):
                if k not in rcs: logger.warning(l[0]+ ' do not align with any query sequence and cannot be analysed. Remove all unplaced scaffolds and contigs from the assemblies.')
    except Exception as e:
        logger.error('Error in reading SAM file: ' + str(e))
        # sys.exit()
        return pd.DataFrame()
    
    # 如果没有比对上，返回空值
    if len(list(al)) == 0:
        return pd.DataFrame()

    al = DataFrame(list(al))
    al[6] = al[6].astype('float')
    al.sort_values([9,0,1,2,3,10], inplace = True, ascending=True)
    al.index = range(len(al.index))
    return al
# END


def readSAMBAM(fin, type='B'):
    import pysam
    logger = logging.getLogger('Reading BAM/SAM file')
    try:
        if type == 'B':
            findata = pysam.AlignmentFile(fin,'rb')
        elif type == 'S':
            return samtocoords(fin)
        else:
            raise ValueError("Wrong parameter")
    except ValueError as e:
        logger.error("Error in opening BAM/SAM file. " + str(e))
        # sys.exit()
        return pd.DataFrame()
    except OSError as e:
        logger.error("Error in reading input file." + str(e))
        # sys.exit()
        return pd.DataFrame()
    except Exception as e:
        logger.error("Unexpected error in opening BAM/SAM file. " + str(e))
        # sys.exit()
        return pd.DataFrame()

    try:
        qry_prim = {}
        ref_prim = {}
        cgdict = {1:'I', 2:'D', 7:'=', 8:'X'}
        coords = {}
        index = 0
        for aln in findata:
            index += 1
            ## Check whether every sequence has at least one primary alignment
            if aln.reference_name is not None:
                if aln.reference_name not in ref_prim.keys():
                    ref_prim[aln.reference_name] = False
            if aln.query_name not in qry_prim.keys():
                qry_prim[aln.query_name] = False
            if aln.reference_name is not None:
                if not ref_prim[aln.reference_name]:
                    if aln.flag < 256:
                        ref_prim[aln.reference_name] = True
            if not qry_prim[aln.query_name]:
                if aln.flag < 256:
                    qry_prim[aln.query_name] = True

            ## Pass non-alinging chromosomes
            if aln.cigarstring is None:
                logger.warning(aln.query_name + ' do not align with any reference chromosome and cannot be analysed')
                continue

            ## Check CIGAR:
            if False in [False if i[0] not in [1,2,4,5,7,8] else True for i in aln.cigartuples]:
                logger.error("Incorrect CIGAR string found. CIGAR string can only have I/D/H/S/X/=. CIGAR STRING: " + str(aln.cigarstring))
                # sys.exit()
                return pd.DataFrame()
            if len(aln.cigartuples) > 2:
                if True in [True if i[0] in [4,5] else False for i in aln.cigartuples[1:-1]]:
                    logger.error("Incorrect CIGAR string found. Clipped bases inside alignment. H/S can only be in the terminal. CIGAR STRING: " + aln.cigarstring)
                    # sys.exit()
                    return pd.DataFrame()

            ## Parse information from the aln object
            astart = aln.reference_start+1
            aend = aln.reference_end
            is_inv = True if np.binary_repr(aln.flag,12)[7] == '1' else False
            if not is_inv:
                if aln.cigartuples[0][0] in [4,5]:
                    bstart = aln.cigartuples[0][1]+1
                else:
                    bstart = 1
                bend = bstart + aln.query_alignment_length - 1
            else:
                if aln.cigartuples[-1][0] in [4,5]:
                    bend = aln.cigartuples[-1][1]+1
                else:
                    bend = 1
                bstart = bend + aln.query_alignment_length - 1
            alen = abs(aend - astart) + 1
            blen = abs(bend - bstart) + 1
            iden = format((sum([i[1] for i in aln.cigartuples if i[0] == 7])/sum([i[1] for i in aln.cigartuples if i[0] in [1,2,7,8]]))*100, '.2f')
            adir = 1
            bdir = -1 if is_inv else 1
            achr = aln.reference_name
            bchr = aln.query_name
            cg = "".join([str(i[1]) + cgdict[i[0]] for i in aln.cigartuples if i[0] not in [4,5]])
            coords[index] = [astart, aend, bstart, bend, alen, blen, iden, adir, bdir, achr, bchr, cg]

        ## Give warning for chromosomes which do not have any primary alignment
        for k,v in ref_prim.items():
            if not v:
                logger.warning('No primary alignment found for reference sequence ' + k +'. This could mean that the entire chromosome '+ k +' is reapeated.')
        for k,v in qry_prim.items():
            if not v:
                logger.warning('No primary alignment found for query sequence ' + k +'. This could mean that the entire chromosome '+ k + ' is reapeated.')

        ## Return alignments
        coords = pd.DataFrame.from_dict(coords, orient= 'index')
        coords.sort_values([9,0,1,2,3,10], inplace = True, ascending=True)
        coords.index = range(len(coords.index))
        coords[6] = coords[6].astype('float')
        return coords
    except Exception as e:
        logger.error("Error in reading BAM/SAM file. " + str(e))
        # sys.exit()
        return pd.DataFrame()
# END


def readPAF(paf):
    coords = deque()
    logger = logging.getLogger('Reading BAM/SAM file')
    try:
        with open(paf, 'r') as fin:
            for line in fin:
                line = line.strip().split()
                astart = int(line[7]) + 1
                aend = int(line[8])
                adir = 1
                bdir = 1 if line[4] == '+' else -1
                bstart = int(line[2]) + 1 if bdir == 1 else int(line[3])
                bend = int(line[3]) if bdir == 1 else int(line[2]) + 1
                alen = abs(aend - astart) + 1
                blen = abs(bend - bstart) + 1 if bdir == 1 else bstart - bend + 1
                cg = [i.split(":")[-1] for i in line[12:] if i[:2] == 'cg']
                if len(cg) != 1:
                    logger.error("CIGAR string is not present in PAF at line {}. Exiting.".format("\t".join(line)))
                    # sys.exit()
                    return pd.DataFrame()
                cg = cg[0]
                ## Check CIGAR:
                if not all([True if i[1] in {'I', 'D', 'H', 'S', 'X', '='} else False for i in cgtpl(cg)]):
                    logger.error("Incorrect CIGAR string found. CIGAR string can only have I/D/H/S/X/=. CIGAR STRING: " + str(cg))
                    # sys.exit()
                    return pd.DataFrame()
                if len(cgtpl(cg)) > 2:
                    if any([True if i[1] in {'H', 'S'} else False for i in cgtpl(cg)]):
                        logger.error("Incorrect CIGAR string found. Clipped bases inside alignment. H/S can only be in the terminal. CIGAR STRING: " + str(cg))
                        # sys.exit()
                        return pd.DataFrame()

                iden = round((sum([int(i[0]) for i in cgtpl(cg) if i[1] == '='])/sum([int(i[0]) for i in cgtpl(cg) if i[1] in {'=', 'X', 'D', 'I'}]))*100, 2)
                achr = line[5]
                bchr = line[0]
                coords.append([astart, aend, bstart, bend, alen, blen, iden, adir, bdir, achr, bchr, cg])
        coords = pd.DataFrame(coords)
        coords.sort_values([9,0,1,2,3,10], inplace = True, ascending=True)
        coords.index = range(len(coords.index))
        coords[6] = coords[6].astype('float')
        return coords
    except FileNotFoundError:
        logger.error("Cannot open {} file. Exiting".format(paf))
        # sys.exit()
        return pd.DataFrame()
    except ValueError as e:
        logger.error("Error in reading PAF: {}. Exiting".format(e))
        # sys.exit()
        return pd.DataFrame()
# END

def readCoords(coordsfin, args, cigar = False):
    logger = logging.getLogger('Reading Coords')
    logger.debug(args.ftype)
    chrlink = {}
    if args.ftype == 'T':
        logger.info("Reading input from .tsv file")
        try:
            coords = pd.read_table(coordsfin, header = None)
        except pd.errors.ParserError:
            coords = pd.read_table(coordsfin, header = None, engine = "python")
        except Exception as e:
            logger.error("Error in reading the alignment file. " + e)
            # sys.exit()
            return pd.DataFrame(), {}
    elif args.ftype == 'S':
        logger.info("Reading input from SAM file")
        try:
            coords = readSAMBAM(coordsfin, type='S')
        except Exception as e:
            logger.error("Error in reading the alignment file. " + e)
            # sys.exit()
            return pd.DataFrame(), {}
    elif args.ftype == 'B':
        logger.info("Reading input from BAM file")
        try:
            coords = readSAMBAM(coordsfin, type='B')
        except Exception as e:
            logger.error("Error in reading the alignment file" + e)
            # sys.exit()
            return pd.DataFrame(), {}
    elif args.ftype == 'P':
        logger.info("Reading input from PAF file")
        try:
            coords = readPAF(coordsfin)
        except Exception as e:
            logger.error("Error in reading the alignment file" + e)
            # sys.exit()
            return pd.DataFrame(), {}
    else:
        logger.error("Incorrect alignment file type specified.")
        # sys.exit()
        return pd.DataFrame(), {}
    
    # 如果坐标为空，返回空值
    if len(list(coords)) == 0:
        return coords, chrlink

    if not cigar:
        if coords.shape[1] >= 12:
            coords = coords.iloc[:, 0:11]
        coords.columns = ["aStart","aEnd","bStart","bEnd","aLen","bLen","iden","aDir","bDir","aChr","bChr"]
    else:
        if coords.shape[1] > 12:
            coords = coords.iloc[:, 0:12]
        coords.columns = ["aStart","aEnd","bStart","bEnd","aLen","bLen","iden","aDir","bDir","aChr","bChr", 'cigar']

    # Sanity check input file
    try:
        coords.aStart = coords.aStart.astype('int')
    except ValueError:
        logger.error('astart is not int')
        # sys.exit()
        return pd.DataFrame(), {}

    try:
        coords.aEnd = coords.aEnd.astype('int')
    except ValueError:
        logger.error('aend is not int')
        # sys.exit()
        return pd.DataFrame(), {}

    try:
        coords.bStart = coords.bStart.astype('int')
    except ValueError:
        logger.error('bstart is not int')
        # sys.exit()
        return pd.DataFrame(), {}

    try:
        coords.bEnd = coords.bEnd.astype('int')
    except ValueError:
        logger.error('abend is not int')
        # sys.exit()
        return pd.DataFrame(), {}

    try:
        coords.aLen = coords.aLen.astype('int')
    except ValueError:
        logger.error('alen is not int')
        # sys.exit()
        return pd.DataFrame(), {}

    try:
        coords.bLen = coords.bLen.astype('int')
    except ValueError:
        logger.error('blen is not int')
        # sys.exit()
        return pd.DataFrame(), {}

    try:
        coords.iden = coords.iden.astype('float')
    except ValueError:
        logger.error('iden is not float')
        # sys.exit()
        return pd.DataFrame(), {}

    try:
        coords.aDir = coords.aDir.astype('int')
    except ValueError:
        logger.error('aDir is not int')
        # sys.exit()
        return pd.DataFrame(), {}

    if any(coords.aDir != 1):
        logger.error('aDir can only have values 1')
        # sys.exit()
        return pd.DataFrame(), {}

    try:
        coords.bDir = coords.bDir.astype('int')
    except ValueError:
        logger.error('bDir is not int')
        # sys.exit()
        return pd.DataFrame(), {}

    for i in coords.bDir:
        if i not in [1,-1]:
            logger.error('bDir can only have values 1/-1')
            # sys.exit()
            return pd.DataFrame(), {}

    try:
        coords.aChr = coords.aChr.astype(str)
    except:
        logger.error('aChr is not string')
        # sys.exit()
        return pd.DataFrame(), {}

    try:
        coords.bChr = coords.bChr.astype(str)
    except:
        logger.error('bChr is not string')
        # sys.exit()
        return pd.DataFrame(), {}

    # Filter small alignments
    if args.f:
        logger.info('Filtering alignments')
        logger.debug('Number of alignments before filtering: {}'.format(coords.shape[0]))
        coords = coords.loc[coords.iden > 90]
        coords = coords.loc[(coords.aLen>100) & (coords.bLen>100)]
        logger.debug('Number of alignments after filtering: {}'.format(coords.shape[0]))

    ## check for bstart > bend when bdir is -1
    check = np.unique(coords.loc[coords.bDir == -1, 'bStart'] > coords.loc[coords.bDir == -1, 'bEnd'])
    if len(check) > 1:
        logger.error('Inconsistent start and end position for inverted alignment in query genome. For inverted alignments, either all bstart < bend or all bend > bstart')
        # sys.exit()
        return pd.DataFrame(), {}
    elif len(check) == 0:
        logger.info('No Inverted alignments present.')
    elif check[0] == True:
        pass
    else:
        logger.info('For inverted alignments, bstart was less than bend. Swapping them.')
        coords.loc[coords.bDir == -1, 'bStart'] = coords.loc[coords.bDir == -1, 'bStart'] + coords.loc[coords.bDir == -1, 'bEnd']
        coords.loc[coords.bDir == -1, 'bEnd'] = coords.loc[coords.bDir == -1, 'bStart'] - coords.loc[coords.bDir == -1, 'bEnd']
        coords.loc[coords.bDir == -1, 'bStart'] = coords.loc[coords.bDir == -1, 'bStart'] - coords.loc[coords.bDir == -1, 'bEnd']

    coords.sort_values(['aChr', 'aStart', 'aEnd', 'bChr', 'bStart', 'bEnd'], inplace=True)

    ## Ensure that chromosome IDs are same for the two genomes.
    ## Either find best query match for every reference genome.
    ## Or if --no-chrmatch is set then remove non-matching chromosomes.
    if np.unique(coords.aChr).tolist() != np.unique(coords.bChr).tolist():
        logger.warning('Chromosomes IDs do not match.')
        if not args.chrmatch:
            if len(np.unique(coords.aChr)) != len(np.unique(coords.bChr)):
                logger.error("Error: unequal number of chromosomes in the genomes. Exiting")
                # sys.exit(1)
                return pd.DataFrame(), {}
            else:
                logger.warning("Matching them automatically. For each reference genome, most similar query genome will be selected. Check mapids.txt for mapping used.")
                chromMaps = defaultdict(dict)
                for i in np.unique(coords.bChr):
                    for j in np.unique(coords.aChr):
                        a = np.array(coords.loc[(coords.bChr == i) & (coords.aChr == j), ["aStart", "aEnd"]])
                        a = mergeRanges(a)
                        chromMaps[j][i] = len(a) + (a[:, 1] - a[:, 0]).sum()

                assigned = []
                fout = open(args.dir+args.prefix+"mapids.txt", "w")
                for chrom in np.unique(coords.aChr):
                    maxid = max(chromMaps[chrom].items(), key=lambda x: x[1])[0]
                    if maxid in assigned:
                        logger.error("{} in genome B is best match for two chromosomes in genome A. Cannot assign chromosomes automatically.".format(maxid))
                        fout.close()
                        os.remove(args.dir+args.prefix+"mapids.txt")
                        # sys.exit()
                        return pd.DataFrame(), {}
                    assigned.append(maxid)
                    fout.write(chrom+"\t"+maxid+"\n")
                    logger.info("setting {} as {}".format(maxid, chrom))
                    coords.loc[coords.bChr == maxid, "bChr"] = chrom
                    chrlink[maxid] = chrom
                fout.close()
        else:
            logger.warning("--no-chrmatch is set. Not matching chromosomes automatically.")
            aChromo = set(coords["aChr"])
            bChromo = set(coords["bChr"])
            badChromo = list(aChromo - bChromo) + list(bChromo - aChromo)
            if len(badChromo) > 0:
                logger.warning(", ".join(badChromo) + " present in only one genome. Removing corresponding alignments")
            coords = coords.loc[~coords.aChr.isin(badChromo) & ~coords.bChr.isin(badChromo)]

    ## Check for presence of directed alignments
    achrs = np.unique(coords.aChr).tolist()
    for achr in achrs:
        if coords.loc[(coords.aChr==achr) & (coords.bChr==achr) & (coords.bDir == 1),].shape[0] == 0:
            hombchr = [k for k,v in chrlink.items() if v==achr]
            if len(hombchr) == 1:
                hombchr = hombchr[0]
            elif len(hombchr) == 0:
                hombchr = achr
            else:
                logger.error('Homologous chromosomes were not identified correctly. Try assigning the chromosome ids manually.')
                # sys.exit()
                return pd.DataFrame(), {}
            # logger.warning('Reference chromosome ' + achr + ' do not have any directed alignments with its homologous chromosome in the query genome (' + hombchr + '). Filtering out all corresponding alignments.')
            coords = coords.loc[~(coords.aChr == achr)]
            coords = coords.loc[~(coords.bChr == achr)]

    ## Check for presence of too many inverted alignments
    for achr in achrs:
        dir_range = mergeRanges(np.array(coords.loc[(coords.aChr==achr) & (coords.bChr==achr) & (coords.bDir==1), ["aStart", "aEnd"]]))
        dir_len = len(dir_range) + (dir_range[:, 1] - dir_range[:, 0]).sum()
        inv_range = mergeRanges(np.array(coords.loc[(coords.aChr==achr) & (coords.bChr==achr) & (coords.bDir==-1), ["aStart", "aEnd"]]))
        inv_len = len(inv_range) + (inv_range[:, 1] - inv_range[:, 0]).sum()
        if inv_len > dir_len:
            hombchr = [k for k,v in chrlink.items() if v==achr]
            if len(hombchr) == 1:
                hombchr = hombchr[0]
            elif len(hombchr) == 0:
                hombchr = achr
            else:
                logger.error('Homologous chromosomes were not identified correctly. Try assigning the chromosome ids manually.')
                # sys.exit()
                return pd.DataFrame(), {}
            # logger.warning('Reference chromosome ' + achr + ' has high fraction of inverted alignments with its homologous chromosome in the query genome (' + hombchr + '). Ensure that same chromosome-strands are being compared in the two genomes, as different strand can result in unexpected errors.')
    return coords, chrlink
# END

def mergeRanges(ranges):
    """
    Take a 2D numpy array, with each row as a range and return merged ranges i.e. ranges which are overlapping would be
    combined as well.
    :param ranges:
    :return:
    """
    from collections import deque
    import numpy as np
    if len(ranges) < 2:
        return ranges
    ranges = ranges[ranges[:, 0].argsort()]
    for i in ranges:
        if i[0] > i[1]:
            garb = i[0]
            i[0] = i[1]
            i[1] = garb
    min_value = ranges[0, 0]
    max_value = ranges[0, 1]
    out_range = deque()
    for i in ranges[1:]:
        if i[0] > max_value:
            out_range.append([min_value, max_value])
            min_value = i[0]
            max_value = i[1]
        elif i[1] > max_value:
            max_value = i[1]
    out_range.append([min_value, max_value])
    return np.array(out_range)


########################################################################################################################
def apply_TS(astart, aend, bstart, bend, threshold, mxgap = 100000000000):
    # cdef:
    i = len(astart)
    j = len(astart)
    n = len(astart)
    df = {}
    # cpp_map[long, cpp_deq[long]].iterator   mapit
    for i in range(n):
        for j in range(i+1, n):
            if (astart[j] - aend[i]) < mxgap:        # Select only alignments with small gaps
                if (astart[j] - astart[i]) > threshold:
                    if (aend[j] - aend[i]) > threshold:
                        if (bstart[j] - bend[i]) < mxgap:        # Select only alignments with small gaps
                            if (bstart[j] - bstart[i]) > threshold:
                                if (bend[j] - bend[i]) > threshold:
                                    if i not in df:
                                        df[i] = []
                                    df[i].append(j)
    out = {}
    for i in range(n):
        if i in df:
            out[i] = [df[i][j] for j in range(len(df[i]))]
        else:
            out[i] = []
    return out


########################################################################################################################
def getSynPath(blocks):
    synPath = []
    scores = [block.score for block in blocks]
    lastBlock = scores.index(max(scores))
    while blocks[lastBlock].bestParentID != -1:
        synPath.append(lastBlock)
        lastBlock = blocks[lastBlock].bestParentID        
    synPath.append(lastBlock)
    return(synPath[::-1])


################################################# 读取文件，获取coords坐标 #################################################
def getCoords(
    args, 
    alignMentFileName, 
    refFileName = "", 
    qryFileName = ""
):
    logger = logging.getLogger('getCoords')

    # chrlink is a dict with query genome ID as key and matching reference genome as values
    coords, chrlink = readCoords(alignMentFileName, args)

    # 如果坐标为空，返回空值
    if len(list(coords)) == 0:
        return coords, chrlink

    achrs = np.unique(coords.aChr).tolist()
    bchrs = np.unique(coords.bChr).tolist()
    for chrid in achrs + bchrs:
        try:
            i = int(chrid)
            logger.error('Numerical chromosome id: {} found in alignments. Please change it to string. Example: for chromosome 1, ">Chr1" is valid but ">1" is invalid. Exiting.'.format(chrid))
            # sys.exit()
            return pd.DataFrame(), {}
        except ValueError as e:
            pass

    achr_size = {}
    bchr_size = {}
    for achr in achrs:
        achr_size[achr] = coords.loc[coords.aChr == achr, ['aStart', 'aEnd']].max().max()
        bchr_size[achr] = coords.loc[coords.bChr == achr, ['bStart', 'bEnd']].max().max()

    key_found = []
    if refFileName:
        for chrid, seq in readfasta(refFileName).items():
            if chrid in achrs:
                try:
                    i = int(chrid)
                    logger.error('Numerical chromosome id found: {} in reference genome. Please change it to string. Example: for chromosome 1, ">Chr1" is valid but ">1" is invalid. Exiting.'.format(chrid))
                    # sys.exit()
                    return pd.DataFrame(), {}
                except ValueError as e:
                    pass
                key_found = key_found + [chrid]
                if len(seq) < achr_size[chrid]:
                    logger.error('Length of reference sequence of ' + chrid + ' is less than the maximum coordinate of its aligned regions. Exiting.')
                    # sys.exit()
                    return pd.DataFrame(), {}
        for achr in achrs:
            if achr not in key_found:
                logger.error('Chromosome ID ' + achr + ' is present in alignments but not in reference genome fasta. Exiting.')
                # sys.exit()
                return pd.DataFrame(), {}
    
    key_found = []
    if qryFileName:
        if len(chrlink) > 0:
            for chrid, seq in readfasta(qryFileName).items():
                if chrid in list(chrlink.keys()):
                    try:
                        i = int(chrid)
                        logger.error('Numerical chromosome id found: {} in query genome. Please change it to string. Example: for chromosome 1, ">Chr1" is valid but ">1" is invalid. Exiting.'.format(chrid))
                        # sys.exit()
                        return pd.DataFrame(), {}
                    except ValueError as e:
                        pass
                    key_found = key_found + [chrid]
                    if len(seq) < bchr_size[chrlink[chrid]]:
                        logger.error('Length of query sequence of ' + chrid + ' is less than the maximum coordinate of its aligned regions. Exiting.')
                        # sys.exit()
                        return pd.DataFrame(), {}
            for bchr in list(chrlink.keys()):
                if bchr not in key_found:
                    logger.error('Chromosome ID ' + bchr + ' is present in alignments but not in query genome fasta. Exiting.')
                    # sys.exit()
                    return pd.DataFrame(), {}
        else:
            for chrid, seq in readfasta(args.qry.name).items():
                if chrid in list(bchr_size.keys()):
                    try:
                        i = int(chrid)
                        logger.error('Numerical chromosome id found: {} in query genome. Please change it to string. Example: for chromosome 1, ">Chr1" is valid but ">1" is invalid. Exiting.'.format(chrid))
                        # sys.exit()
                        return pd.DataFrame(), {}
                    except ValueError as e:
                        pass
                    key_found = key_found + [chrid]
                    if len(seq) < bchr_size[chrid]:
                        logger.error('Length of query sequence of ' + chrid + ' is less than the maximum coordinate of its aligned regions. Exiting.')
                        # sys.exit()
                        return pd.DataFrame(), {}
            for bchr in list(bchr_size.keys()):
                if bchr not in key_found:
                    logger.error('Chromosome ID ' + bchr + ' is available in alignments but not in query genome fasta. Exiting.')
                    # sys.exit()
                    return pd.DataFrame(), {}

    return coords, chrlink



################################################# 获取共线性坐标 #################################################
def getSyn(coords):
    # 保存所有共线性坐标
    synDatas = pd.DataFrame()

    # 如果坐标为空，返回空值
    if len(list(coords)) == 0:
        return synDatas

    threshold = 50  ##args.threshold

    logger = logging.getLogger("syri")
    logger.info("starting")
    logger.debug("memory usage: " + str(psutil.Process(os.getpid()).memory_info()[0]/2.**30))

    uniChromo = list(np.unique(coords.aChr))
    logger.info('Analysing chromosomes: {}'.format(uniChromo))

    # 每个染色体对都遍历
    for chromo in uniChromo:
        logger = logging.getLogger("syri."+chromo)
        coordsData = coords[(coords.aChr == chromo) & (coords.bChr == chromo) & (coords.bDir == 1)]
        logger.info(chromo + " " + str(coordsData.shape))
        logger.info("Identifying Synteny for chromosome " + chromo)
        df = apply_TS(
            coordsData.aStart.values, 
            coordsData.aEnd.values, 
            coordsData.bStart.values, 
            coordsData.bEnd.values, 
            threshold
        )
        blocks = [alignmentBlock(i, df[i], coordsData.iloc[i]) for i in df.keys()]
        for block in blocks:
            i = 0
            while(i < len(block.children)):
                block.children = list(set(block.children) - set(blocks[block.children[i]].children))
                i+=1
            block.children.sort()

            for child in block.children:
                blocks[child].addParent(block.id)

            scores = [blocks[parent].score for parent in block.parents]
            if len(scores) > 0:
                block.bestParent(block.parents[scores.index(max(scores))], max(scores))
        synPath = getSynPath(blocks)
        synData = coordsData.iloc[synPath].copy()
        synDatas = pd.concat([synDatas, synData], ignore_index=True)
        del(coordsData, blocks, df)
        collect()

    return synDatas


########################################################################################################################
def main(args, alignMentFileName):
    # ###### 读取比对文件，获取所有坐标对 ###### #
    coords, chrlink = getCoords(args, alignMentFileName)

    # 如果坐标为空，返回空值
    if len(list(coords)) == 0:
        return pd.DataFrame()

    # ###### 获取共线性坐标 ###### #
    synDatas = getSyn(coords)

    return synDatas
