'''
@File    :   Bam2CuttingBW.py
@Time    :   2020/8/15 4:30 下午
@Author  :   KeeeeepGoing(Yuechao Wu)
@Version :   
@Contact :   wuyuechao@zhangtaolab.org
@Desc    :   copy from https://github.com/forrestzhang/bagatelle.git
'''

import pysam
import sys
import pyBigWig



def dhstobw(bamfile, samplename, excludechr='', library='Duke'):
    # Washington is under processing
    """

    :param bamfile:
    :param bwfile:
    :param library:Duke or Washington

        Duke: |=====>
                        <=====|

        Washington: |===========|

        Out put cutting site '|'
    :return:
    """
    if excludechr == None:
        excludechr = ''

    bwfile = samplename + '_' + 'cutting.bw'

    bamfor = Baminfo(bamfile)

    bw = pyBigWig.open(bwfile, "w")

    excludechrs = excludechr.split(',')

    countchrs = list()

    for chrom in list(bamfor.chrlen.items()):

        if chrom[0] not in excludechrs:
            countchrs.append(chrom)

    bw.addHeader(countchrs)

    print(countchrs)

    for chromosome in bamfor.chrlen:

        if chromosome not in excludechrs:


            end = bamfor.chrlen[chromosome]

            dhscut = dhcutcount(bamfile=bamfile, chromosome=chromosome, start=1,
                                       end=end, library=library)

            if dhscut:

                starts = list()

                values = list()

                for start in sorted(dhscut):
                    starts.append(start)

                    values.append(float(dhscut[start]))

                bw.addEntries(chromosome, starts=starts, values=values,
                              span=1, step=1)

    bw.close()



def openBam(bamFile):
    try:
        bam = pysam.Samfile(bamFile, 'rb')

    except IOError:

        sys.exit("The file {} does not exist".format(bamFile))

    except:

        sys.exit("The file {} does not have BAM format ".format(bamFile))

    try:

        if 'check_index' in dir(bam):

            assert (bam.check_index())

        else:
            # The proper check_index() function wasn't implemented until pysam 0.8.4!
            assert (bam._hasIndex())

    except:

        sys.exit("{} does not appear to have an index. You MUST index the file first!".format(bamFile))

    if bam.mapped == 0:
        sys.exit("Samtools reports that the number of mapped "
                 "reads is zero for the file {}. Please check "
                 "that the file is properly indexed and that "
                 "it contains mapped reads.".format(bamFile))

    return bam



class Baminfo:


    def __init__(self, bamfile):
        self.bamfile = bamfile

        self.samfile = openBam(self.bamfile)

        self.chrlen = self.getchrlen()

    def getchrlen(self):
        ref_lengths = self.samfile.lengths

        sam_ref = self.samfile.references

        refere_ncenumber = self.samfile.nreferences

        chrlen = dict()

        for i in range(refere_ncenumber):
            chr_length = ref_lengths[i]

            chrlen[sam_ref[i]] = chr_length

        return chrlen


def dhcutcount(bamfile, chromosome, start, end, library='Duke'):
    """

    :param bamfile: bamfile
    :param chromosome:
    :param start:
    :param end:
    :param library: Duke or Washington

        Duke: |=====>
                        <=====|

        Washington: |===========|

        Out put cutting site '|'

    :return: dictionary of cutting site count
    """

    samfile = openBam(bamfile)

    readscount = dict()

    if library == 'Duke':

        for aligned_read in samfile.fetch(reference=str(chromosome), start=start, end=end):

            if aligned_read.is_reverse:

                site = aligned_read.aend

            else:

                site = aligned_read.pos

            site = site + 1

            if site in readscount:

                readscount[site] = readscount[site] + 1

            else:

                readscount[site] = 1

    elif library == 'Washington':

        pass

    else:

        pass

    return readscount


