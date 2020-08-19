'''
@File    :   FDR_calu_write.py
@Time    :   2020/8/15 4:10 下午
@Author  :   KeeeeepGoing(Yuechao Wu)
@Version :   
@Contact :   wuyuechao@zhangtaolab.org
@Desc    :   
'''
from typing import Dict, Any

import numpy as np
from decimal import Decimal
import pyBigWig
import scipy.stats as stats
from statsmodels.stats import multitest
from multiprocessing import Pool


def calculate_FDR(hotspots, CuttingBwFile, nthreads, boundary_width=50):
    """

    :param hotspots:
    :param CuttingBwFile:
    :param countchr:
    :param boundary_width:
    :return: 返回全部peak中每个碱基的FDR
    """

    hotspot_list = []
    CuttingBwObject = pyBigWig.open(CuttingBwFile)

    for hotspot in hotspots:
        chrom = hotspot.chromosome
        start = hotspot.start
        end = hotspot.end

        hotspot_list.append({'chrom':chrom,'start':start, 'end':end, 'chromlength':CuttingBwObject.chroms(chrom), 'CuttingBwFile':CuttingBwFile, 'boundary_width':boundary_width})

    pool = Pool(nthreads)

    hotspot_pvalue_dict = pool.map(calculate_FDR_threads, hotspot_list)

    pool.close()

    return  hotspot_pvalue_dict





def calculate_FDR_threads(hotspot_list):

    chrom = hotspot_list['chrom']
    chromlength = hotspot_list['chromlength']
    start = hotspot_list['start']
    end = hotspot_list['end']
    CuttingBwFile = hotspot_list['CuttingBwFile']
    boundary_width = hotspot_list['boundary_width']
    CuttingBwObject = pyBigWig.open(CuttingBwFile)

    hotspot_pvalue_dict = {}

    hotspot_id = '_'.join([chrom, str(start), str(end)])

    hotspot_pvalue_dict[hotspot_id] = []

    hotspots_localmu = GetPeakLocalMU(CuttingBwObject=CuttingBwObject, chrom=chrom, start=start, end=end,
                                      extendsize=50000)

    for base in range(start, end):
        boundary_left = int(base - boundary_width / 2)
        boundary_right = int(base + boundary_width / 2)

        boundary_left = 1 if boundary_left < 1 else boundary_left
        boundary_right = chromlength if boundary_right > chromlength else boundary_right

        regioncount = GetBaseLocalCount(CuttingBwObject=CuttingBwObject, chrom=chrom, boundary_left=boundary_left,
                                        boundary_right=boundary_right)

        base_p_value = poissonpvalue(regioncount, hotspots_localmu * boundary_width)

        hotspot_pvalue_dict[hotspot_id].append(base_p_value)

    hotspot_pvalue_dict[hotspot_id] = list(multitest.multipletests(hotspot_pvalue_dict[hotspot_id], method='fdr_bh')[1])

    return hotspot_pvalue_dict






# def calculate_FDR(hotspots, CuttingBwFile, countchr, boundary_width=50):
#     """
#
#     :param hotspots:
#     :param CuttingBwFile:
#     :param countchr:
#     :param boundary_width:
#     :return: 返回全部peak中每个碱基的FDR
#     """
#     base_FDR_dict = { i:{} for i in countchr}
#
#     CuttingBwObject = pyBigWig.open(CuttingBwFile)
#
#     hotspot_pvalue_dict = {}
#
#     for hotspot in hotspots:
#
#         chrom = hotspot.chromosome
#         chromlength = CuttingBwObject.chroms(chrom)
#         start = hotspot.start
#         end = hotspot.end
#         hotspot_id = '_'.join([chrom, str(start), str(end)])
#
#         hotspot_pvalue_dict[hotspot_id] = []
#
#         hotspots_localmu = GetPeakLocalMU(CuttingBwObject=CuttingBwObject, chrom=chrom, start=start, end=end, extendsize=50000)
#
#         for base in range(start, end):
#
#             boundary_left = int(base - boundary_width/2)
#             boundary_right = int(base + boundary_width/2)
#
#             boundary_left = 1 if boundary_left < 1 else boundary_left
#             boundary_right = chromlength if boundary_right > chromlength else boundary_right
#
#             regioncount = GetBaseLocalCount(CuttingBwObject=CuttingBwObject, chrom=chrom, boundary_left=boundary_left, boundary_right=boundary_right)
#
#             base_p_value = poissonpvalue(regioncount, hotspots_localmu*boundary_width)
#
#             hotspot_pvalue_dict[hotspot_id].append(base_p_value)
#
#
#     for hotspot_id in hotspot_pvalue_dict :
#
#         hotspot_pvalue_list = hotspot_pvalue_dict[hotspot_id]
#         hotspot_pvalue_dict[hotspot_id] = list(multitest.multipletests(hotspot_pvalue_list, method='fdr_bh')[1])
#
#     return hotspot_pvalue_dict


def write_hotspot_FDR(hotspot_pvalue_dict, baseFDRFile):

    outfile = open(baseFDRFile, 'w')

    for hotspot_dict in hotspot_pvalue_dict:

        hotspot_chr = ''.join(hotspot_dict.keys()).split('_')[0]
        hotspot_start = int(''.join(hotspot_dict.keys()).split('_')[1])
        hotspot_end = int(''.join(hotspot_dict.keys()).split('_')[2])

        base_FDR_list = list(hotspot_dict.values())[0]

        count = 0

        for base in range(hotspot_start, hotspot_end):
            print(hotspot_chr, str(base), str(base+1), 'l', str(base_FDR_list[count]), sep='\t', file=outfile)
            count += 1

    outfile.close()






def GetPeakLocalMU(CuttingBwObject, chrom, start, end, extendsize=50000):
    """

    :param CuttingBwObject:
    :param chrom:
    :param start:
    :param end:
    :param extendsize:
    :return: 返回指定区域上下游50k切点位置的均值
    """
    chrlengh = CuttingBwObject.chroms(chrom)

    start = start - extendsize

    end = end + extendsize

    start = 1 if start < 1 else start

    end = chrlengh if end > chrlengh else end

    localmu = np.mean(np.nan_to_num(CuttingBwObject.values(chrom, start, end, numpy=True)))

    return localmu


def GetBaseLocalCount(CuttingBwObject, chrom, boundary_left, boundary_right):
    """

    :param CuttingBwObject:
    :param chrom:
    :param boundary_left:
    :param boundary_right:
    :return: 返回位点上下游boundary_left的切点总和
    """
    regioncout = np.sum(np.nan_to_num(CuttingBwObject.values(chrom, boundary_left, boundary_right, numpy=True)))

    return regioncout


def poissonpvalue(x, mu):
    """

    :param x:
    :param mu:
    :return: 返回x在以mu为均数泊松分布的概率
    """
    poissonpvalue = 1 - stats.poisson.cdf(x, mu)

    return poissonpvalue



