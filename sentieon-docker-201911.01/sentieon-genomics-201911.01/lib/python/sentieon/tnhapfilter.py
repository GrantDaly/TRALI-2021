#!/usr/bin/env python
from __future__ import print_function
import argparse
import math
import sys
import vcflib
import os
from vcflib.compat import *

class TNhapFilter:

    params = {
        'contamination' : [ 0., "Contamination fraction to filter", ""],
        'max_alt_cnt' : [ 1, "Maximum alt allele count", ""],
        'max_event_cnt' : [ 2, "Maximum events in region", ""],
        'min_median_base_qual' : [ 20, "Minimum median base quality", ""],
        'min_median_mapq' : [ 30, "Minimum median map", ""],
        'max_diff_fraglen' : [ 10000, "Maximum median difference in fragment length", ""],
        'min_pir_median' : [ 5, "Minimum median read position", ""],
        'max_normal_art_lod' : [ 0., "Maximum normal artifact LOD", ""],
        'min_tumor_lod' : [ 5.3, "Minimum tumorLOD", ""],
        'max_strand_prob' : [ 0.99, "Maximum strand artifact probability", ""],
        'max_germline_prob' : [ 0.025, "Maximum germline probability", ""],
        'min_strand_af' : [ 0.01, "Minimum strand artifact allele fraction", ""],
    }

    @classmethod
    def add_arguments(cls, parser):
        for k,v in cls.params.items():
            parser.add_argument('--'+k, default=v[0], type=type(v[0]), help=v[1] + ' (default: '+str(v[0])+')', metavar=v[2])

    def __init__(self, args, t_smid, n_smid):
        self.args = args
        self.t_smid = t_smid
        self.n_smid = n_smid

    def applyFilters(self, v):
        filters = []
        tumor_ann = v.samples[self.t_smid]

        tlods = v.info.get('TLOD', [])
        max_tlod_idx = -1
        if tlods and isinstance(tlods, list):
            max_tlod_idx = tlods.index(max(tlods))
        
            # insufficientEvidenceFilter 
            if tlods[max_tlod_idx] < self.args.min_tumor_lod:
                filters.append('t_lod')

        # clusteredEventFilter
        if 'ECNT' in v.info and v.info['ECNT'] > self.args.max_event_cnt:
            filters.append('clustered_events')
        
        # duplicatedAltReadFilter 
        
        # triallelicFilter
        numPass = 0
        if isinstance(tlods, list):
            for tlod in tlods:
                if tlod > self.args.min_tumor_lod:
                    numPass += 1
        if numPass > self.args.max_alt_cnt:
            filters.append('multiallelic')
          
        # PanelOfNormalsFilter
        if 'IN_PON' in v.info:
            filters.append('panel_of_normals')

        # germlineVariantFilter         
        p_germline = v.info.get('P_GERMLINE', [])
        if p_germline and max_tlod_idx != -1 and len(p_germline) > max_tlod_idx:
            log10GermlinePost = p_germline[max_tlod_idx]
            if log10GermlinePost > math.log10(self.args.max_germline_prob):
                filters.append('germline_risk')
         
        # artifactInNormalFilter
        if 'N_ART_LOD' in v.info:
            n_art_lods = v.info['N_ART_LOD']
            if max_tlod_idx != -1 and len(n_art_lods) > max_tlod_idx and n_art_lods[max_tlod_idx] > self.args.max_normal_art_lod:
                filters.append('artifact_in_normal')
        
        # STRFilter
        if 'RPA' in v.info and 'RU' in v.info:
            rpa = v.info['RPA']
            ru = v.info['RU']
            if len(ru) > 1 and isinstance(rpa, list) and len(rpa) > 1:
                refCnt = rpa[0]
                altCnt = rpa[1]
                if (refCnt - altCnt == 1):
                    filters.append('str_contraction')
        
        # strandArtifactFilter
        if 'SA_POST_PROB' in tumor_ann and 'SA_MAP_AF' in tumor_ann:
            sa_post_prob = tumor_ann['SA_POST_PROB']            
            sa_map_af = tumor_ann['SA_MAP_AF']
            if sa_post_prob and isinstance(sa_post_prob, list) and sa_map_af:
                max_sa_idx = sa_post_prob.index(max(sa_post_prob))                
                if max_sa_idx != 2:
                    if sa_post_prob[max_sa_idx] > self.args.max_strand_prob and sa_map_af[max_sa_idx] < self.args.min_strand_af:
                        filters.append('strand_artifact')  
        
        #contaminationFilter
        if 'AF' in tumor_ann:
            af_tumor = tumor_ann['AF']
            if isinstance(af_tumor, list) and af_tumor:
                max_af = max(af_tumor)
                if (max_af < self.args.contamination):
                    filters.append('contamination')
            
        # medianBaseQualityDifferenceFilter
        if 'MBQ' in tumor_ann:
            abq = tumor_ann['MBQ']
            if abq and isinstance(abq, list) and isinstance(abq[0], int) and abq[0] < self.args.min_median_base_qual: 
                filters.append('base_quality')

        # medianMappingQualityDifferenceFilter
        if 'MMQ' in tumor_ann:
            amq = tumor_ann['MMQ']
            if amq and isinstance(amq, list) and isinstance(amq[0], int) and amq[0] < self.args.min_median_mapq: 
                filters.append('mapping_quality')

        # medianFragmentLengthDifferenceFilter
        if 'MFRL' in tumor_ann:
            afl = tumor_ann['MFRL']
            if afl and isinstance(afl, list) and len(afl)>1 and math.fabs(afl[1] - afl[0]) > self.args.max_diff_fraglen:
                filters.append('fragment_length')
    
        # readPositionFilter
        if 'MPOS' in tumor_ann:
            t_gt = []
            if 'GT' in tumor_ann:
                t_gt = tumor_ann.get('GT').split('/')
            n_gt = []
            if 'GT' in v.samples[self.n_smid]:
                n_gt = v.samples[self.n_smid].get('GT').split('/')
            count = {}
            for gt in t_gt + n_gt:
                if gt != '0':
                    count[gt] = count.get(gt, 0) + 1 
            if count:
                max_count_gt = max(count, key=count.get)
                max_t_gt_idx = -1
                if max_count_gt in t_gt:
                    max_t_gt_idx = t_gt.index(max_count_gt)
                if max_t_gt_idx >= 1:
                    arp = tumor_ann['MPOS']
                    insertSize = max(0, len(v.alt[max_t_gt_idx-1]) - len(v.ref))            
                    if arp and ( insertSize + arp[0] < self.args.min_pir_median ):
                        filters.append('read_position')

        flds = v.line.split('\t')
        flds[6] = filters and ';'.join(sorted(set(filters))) or 'PASS'
        v.line = '\t'.join(flds)
        return v

extra_headers = (
    '##FILTER=<ID=PASS,Description="Accept as a confident somatic mutation">',
    '##FILTER=<ID=artifact_in_normal,Description="artifact_in_normal">',
    '##FILTER=<ID=base_quality,Description="alt median base quality">',
    '##FILTER=<ID=clustered_events,Description="Clustered events observed in the tumor">',
    '##FILTER=<ID=contamination,Description="contamination">',
    '##FILTER=<ID=duplicate_evidence,Description="evidence for alt allele is overrepresented by apparent duplicates">',
    '##FILTER=<ID=fragment_length,Description="abs(ref - alt) median fragment length">',
    '##FILTER=<ID=germline_risk,Description="Evidence indicates this site is germline, not somatic">',
    '##FILTER=<ID=mapping_quality,Description="ref - alt median mapping quality">',
    '##FILTER=<ID=multiallelic,Description="Site filtered because too many alt alleles pass tumor LOD">',
    '##FILTER=<ID=panel_of_normals,Description="Blacklisted site in panel of normals">',
    '##FILTER=<ID=read_position,Description="median distance of alt variants from end of reads">',
    '##FILTER=<ID=str_contraction,Description="Site filtered due to contraction of short tandem repeat region">',
    '##FILTER=<ID=strand_artifact,Description="Evidence for alt allele comes from one read direction only">',
    '##FILTER=<ID=t_lod,Description="Tumor does not meet likelihood threshold">',
    '##filtering_status=These calls have been filtered by FilterMutectCalls to label false positives with a list of failed filters and true positives with PASS.',
)

remove_headers = ('##filtering_status=*', '##source=*')

expect_types = {
    'AF': {'Number': 'A', 'Type': 'Float' },
    'GT': {'Number': '1', 'Type': 'String' },
    'MBQ': {'Number': 'A', 'Type': 'Integer' },
    'MFRL': {'Number': 'R', 'Type': 'Integer' },
    'MMQ': {'Number': 'A', 'Type': 'Integer' },
    'MPOS': {'Number': 'A', 'Type': 'Integer' },
    'SA_MAP_AF': {'Number': '3' , 'Type': 'Float'},
    'SA_POST_PROB': {'Number': '3', 'Type': 'Float'},
    'ECNT': {'Number': '1', 'Type': 'Integer' },
    'IN_PON': {'Number': '0', 'Type': 'Flag' },
    'N_ART_LOD': {'Number': 'A', 'Type': 'Float' },
    'P_GERMLINE': {'Number': 'A', 'Type': 'Float' },
    'RPA': {'Number': '.', 'Type': 'Integer' },
    'RU': {'Number': '1', 'Type': 'String' },
    'TLOD': {'Number': 'A', 'Type': 'Float' },
}

def check_header_types(vcf):
    vcf_types = {}
    for k, v in iteritems(vcf.infos):
        if k in expect_types: 
            dd = {} 
            for kk, vv in iteritems(v):
                if kk in ('Number', 'Type'):
                    dd[kk] = vv
                vcf_types[k] = dd     
    for k, v in iteritems(vcf.formats):
        if k in expect_types:
            dd = {}   
            for kk, vv in iteritems(v):                 
                if kk in ('Number', 'Type'):
                    dd[kk] = vv
                vcf_types[k] = dd                    
    return vcf_types == expect_types

def main(args):
    if not os.path.exists(args.vcf):
        print('Error: input file %s does not exist' %  args.vcf)
        return -1

    invcf = vcflib.VCF(args.vcf, 'r')
        
    if not check_header_types(invcf):
        print('Error: vcf format is not expected')
        return -1  
    
    try:
        t_smid = invcf.samples.index(args.tumor_sample)
    except ValueError:
        print('Error: tumor sample "%s" not in input file %s' %
            (args.tumor_sample, args.vcf), file=sys.stderr)
        return -1
    try:
        if args.normal_sample:
            n_smid = invcf.samples.index(args.normal_sample)
        else:
            n_smid = -1
    except ValueError:
        print('Error: normal sample "%s" not in input file %s' %
            (args.normal_sample, args.vcf), file=sys.stderr)
        return -1
    filter = TNhapFilter(args, t_smid, n_smid)

    outvcf = vcflib.VCF(args.output, 'w')
    outvcf.copy_header(invcf, extra_headers, remove_headers)
    outvcf.emit_header()
    for v in invcf:
        outvcf.emit(filter.applyFilters(v))
    outvcf.close()
    invcf.close()
    return 0

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='sentieon tnhapfilter', usage='%(prog)s [options] -v VCF --tumor_sample TUMOR_SAMPLE output')
    parser.add_argument('output', help='Output vcf file name')
    parser.add_argument('-v','--vcf',required=True, help='Input vcf file name')
    parser.add_argument('--tumor_sample', required=True, help='Tumor sample name', metavar="")
    parser.add_argument('--normal_sample', help='Normal sample name', metavar="")
    TNhapFilter.add_arguments(parser)
    sys.exit(main(parser.parse_args()))

# vim: ts=4 sw=4 expandtab
