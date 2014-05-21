from __future__ import print_function
from __future__ import division
import io
from Sampleinfor import *
from countreads import *
from Hotspot import *
from bgcount import *
from multiprocessing import Pool
from kernelsmooth import *
from kernel import *
import numpy as np

class KeyboardInterruptError(Exception):

    pass

def normalizeratio(sampleinfors):

    normailziedratio = dict()

    adjreads = list()

    for sampleinfor in sampleinfors:

        adjreads.append(sampleinfor.fregion.adjreads)

    # minreadscount = min(adjreads)

    for sampleinfor in sampleinfors:

        normailziedratio[sampleinfor.samplename] = sampleinfor.fregion.adjreads/sampleinfor.fregion.countgenomelength

    return normailziedratio


def hotspotswriter(samplename, hotspots):

    bedfilename =samplename+ '_' + 'hotspots' + ".bed"

    open_bed = io.FileIO(bedfilename, 'w')

    for hotspot in hotspots:

        #bedlist = [str(hotspot.chromosome), str(hotspot.start), str(hotspot.end), hotspot.hotspotid]
        bedlist = [str(hotspot.chromosome), str(hotspot.start), str(hotspot.end)]

        linker = "\t"

        outstring = linker.join(bedlist) + "\n"

        open_bed.write(outstring)

    open_bed.close()


def mergedhotsportscountwrite(sampleinfors, mergedhotspots, nthreads, outfilename="merged_hotsopts_count.txt"):

    try:

        outfile = open(outfilename,'w')

        headerstring = "location"

        for sampleinfor in sampleinfors:

            headerstring = headerstring + "\t" + sampleinfor.samplename

        print(headerstring, file=outfile)

        normalizedratio = normalizeratio(sampleinfors=sampleinfors)

        pars = list()

        # for hotspot in mergedhotspots:
        #
        #     outstr = hotspot.region
        #
        #     for sampleinfor in sampleinfors:
        #
        #         nowreads = dhsingleregioncounter(bamfile=sampleinfor.datafile, region=hotspot.region)
        #
        #         normailziedcount = nowreads/(hotspot.end-hotspot.start+1)/normalizedratio[sampleinfor.samplename]
        #
        #         outstr = outstr + "\t" + str(normailziedcount)
        #
        #     print(outstr, file=outfile)
        for hotspot in mergedhotspots:

            par = dict()

            par['hotspot'] = hotspot

            par['sampleinfors'] = sampleinfors

            par['normalizedratio'] = normalizedratio

            pars.append(par)

        pool = Pool(nthreads)



        countinthreads = pool.map(hotspotscounter, pars)

        for countstr in countinthreads:

            print (countstr, file=outfile)

        pool.close()

        outfile.close()

    except KeyboardInterrupt:

        pool.terminate()

        print ("You cancelled the program!")

        sys.exit(1)

    except Exception, e:

        print ('got exception: %r, terminating the pool' % (e,))

        pool.terminate()

        print ('pool is terminated')

    finally:
        #     print ('joining pool processes')
        pool.join()


def hotspotscounter(par):

    try:

        hotspot = par['hotspot']

        sampleinfors = par['sampleinfors']

        outstr = hotspot.region

        normalizedratio = par['normalizedratio']

        for sampleinfor in sampleinfors:

            nowreads = dhsingleregioncounter(bamfile=sampleinfor.datafile, region=hotspot.region)

            normailziedcount = nowreads/(hotspot.end-hotspot.start+1)/normalizedratio[sampleinfor.samplename]

            outstr = outstr + "\t" + str(normailziedcount)

        return outstr

    except KeyboardInterrupt:

        sys.stderr.write("User interrupt\n")

        sys.exit(0)


def wigwritte(sampleinfors, kernellength, nthreads):

    try:

        normalizedratio = normalizeratio(sampleinfors=sampleinfors)

        pars = list()

        for sampleinfor in sampleinfors:

            uniqreate = dhuniquerate(fregion=sampleinfor.fregion)

            maxscore = dhnoncontrol(uniqueratio=uniqreate, threshold=50, kernellength=kernellength, nthreads=nthreads)

            samplenormalizedratio = normalizedratio[sampleinfor.samplename]

            par = dict()

            # par['uniqreate'] = uniqreate

            par['maxscore'] = maxscore

            par['sampleinfor'] = sampleinfor

            par['samplenormalizedratio'] = samplenormalizedratio

            par['kernellength'] = kernellength

            par['maxscore'] = maxscore/samplenormalizedratio

            pars.append(par)

        pool = Pool(nthreads)

        pool.map(wigwritter, pars)

        pool.close()

    except KeyboardInterrupt:

        pool.terminate()

        print ("You cancelled the program!")

        sys.exit(1)

    except Exception, e:

        print ('got exception: %r, terminating the pool' % (e,))

        pool.terminate()

        print ('pool is terminated')

    finally:
        #     print ('joining pool processes')
        pool.join()


def wigwritter(par):

    try:
        sampleinfor = par['sampleinfor']

        kernellength = par['kernellength']

        samplenormalizedratio = par['samplenormalizedratio']

        kernellength = par['kernellength']

        maxscore = par['maxscore']

        bamfile = sampleinfor.datafile

        wigfilename = sampleinfor.samplename+".wig"

        wigio = open(wigfilename, 'w')

        kernel = smooth_kernel(kernellength)

        kernel_score = list()

        for w in sorted(kernel):

            kernel_score.append(kernel[w])

        for chromosome in sampleinfor.fregion.count_chr:

            print('track type=wiggle_0 name="',sampleinfor.samplename,'" description="', sampleinfor.samplename,'"',sep='',file=wigio)

            print('variableStep	chrom=', chromosome, sep='', file=wigio)

            # chrregion = chromosome+":"+str(1)+"-"+str(sampleinfor.fregion.chrs_length[chromosome])
            for scare in range(0, int(sampleinfor.fregion.chrs_length[chromosome]/1000000)+1):

                startsite = scare * 1000000 + 1

                endsite = (scare + 1) * 1000000

                if endsite > sampleinfor.fregion.chrs_length[chromosome]:

                    endsite = sampleinfor.fregion.chrs_length[chromosome]

                regionnow = chromosome+":" + str(startsite) + "-" + str(endsite)

                smoothedscore = regionsmooth(bamfile=bamfile, region=regionnow,
                                             chr_length=sampleinfor.fregion.chrs_length[chromosome],
                                             kernelsize=kernellength)

                for site in sorted(smoothedscore['score'].keys()):

                    score = smoothedscore['score'][site]/samplenormalizedratio

                    if score > maxscore:

                        score = maxscore

                    score = round(score, 3)

                    print (str(site)+"\t"+str(score), file=wigio)

        wigio.close()

    except KeyboardInterrupt:

        print ("You cancelled the program!")

        sys.exit(1)