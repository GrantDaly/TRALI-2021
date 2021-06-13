import sys
import gzip
import json
import pprint
import pdb
import pyranges
import re
import pandas as pd
import pyranges as pr
from jsonpath_ng import jsonpath, parse

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


def calculateHaplotypes(query, subject, qStart, contig, contigStart, contigEnd):

    alignLength = contigEnd - contigStart
    outInsertions = []
    outDeletions = []
    outMutations = []
    #find insertions (deletions in the mito)
    p = re.compile('-+')

    # until I can get correct blast outputs, these intervals won't be right, because I'm expecting an exact match to work back coordinates but indels have "-" character. Going to output the sequences in the mean time.
    for match in p.finditer(query):
        tempStart = match.span()[0]
        tempEnd = match.span()[1]
        numInsLength = tempEnd - tempStart

        tempStart = tempStart + qStart
        outInsertions.append((tempStart, numInsLength, query, subject))

    #find deletions (deletions in numt)
    p = re.compile('-+')
    for match in p.finditer(subject):
        tempStart = match.span()[0]
        tempEnd = match.span()[1]
        numInsLength = tempEnd - tempStart

        tempStart = tempStart + qStart
        outDeletions.append((tempStart, numInsLength, query, subject))

    #mismatches
    #print("New Pair")
    qPos = qStart
    for allele in zip(query, subject):
        qBase, sBase = allele
        qBase = qBase.upper()
        sBase = sBase.upper()
        #print("Query %s | Subject %s Pos %d" % (qBase, sBase, qPos))
        if(qBase == "-"):
            #print("Query -")
            continue
        elif(sBase == "-"):
            #print("Subject -")
            qPos += 1
            continue
        elif(qBase != sBase):
            #print("Mismatch")
            #print("%s\t%s\t%s\t%s" % (qStart, qPos, qBase, sBase))
            outMutations.append((qPos, qBase, sBase))
        qPos += 1
    return {'chrom': contig, 'start': contigStart, 'end':contigEnd,
            'length': alignLength, 'mismatches': outMutations, 'insertions': outInsertions,
            'deletions': outDeletions }

    
if __name__ == '__main__':
    recordList = []
    bedDFList = []
    haplotypeList = []
    pp = pprint.PrettyPrinter(indent=4, depth=4)
    filename = sys.argv[1]
    outPrefix = sys.argv[2]
    with gzip.open(filename, 'rt', encoding='ascii') as zipfile:
        blast_file = json.load(zipfile)
    
    counter = 1
    mitoReference = SeqIO.to_dict(SeqIO.parse("mitochondria.oneline.fa", "fasta"))
    
    for hit in blast_file['BlastOutput2'][0]['report']['results']['search']['hits']: 
        accession = hit['description'][0]['accession']

        if(accession != 'chrM'):
            print(accession)
            for hsp in hit['hsps']:
                #print(hsp.keys())

                #get subject strand
                if (hsp['hit_strand'] == "Plus"):
                    tempStrand = "+"
                elif(hsp['hit_strand'] == "Minus"):
                    tempStrand = "-"
                else:
                    print("invalid strand")

                tempStartSub = hsp['hit_from']
                tempEndSub = hsp['hit_to']

                sortedPositions = sorted((tempStartSub, tempEndSub))
                tempStartSub = sortedPositions[0]
                tempEndSub = sortedPositions[1]
                if(tempStartSub > tempEndSub):
                    print("#######invalid##########")
                    print("%d\t%d % tempStartSub tempEndSub")
                #pdb.set_trace()
                
                #tempStartQuery = hsp['query_from']
                #tempEndQuery = hsp['query_to']
                # blast had been outputting incorrect query coordinates, so I am manually finding them
                tempLength = hsp['align_len']
                percIdentity = ( hsp['identity'] / tempLength ) * 100
                tempScore = hsp['score']
                #pp.pprint(hsp)

                tempSub = hsp['hseq']
                tempQuery = hsp['qseq']

                

                #split the query up by - character and query first segment.
                #not ideal because it may not be the longest segment
                splitList = tempQuery.split("-")
                #longestEntry = sorted(test, key= lambda x : len(x), reverse=True)
                joinedEntry = "".join(splitList)
                # seq.find seems to have upper length limit so using first 50 nucleotides
                queryStart = mitoReference['chrM'].seq.find(joinedEntry[:25]) + 1
                
                #if(queryStart == 0):
                #    pdb.set_trace()
                numtName = "Numt_" + str(counter)
                tempBedDict = {'Chromosome': accession, 'Start': tempStartSub,
                               'End': tempEndSub, 'Name': numtName,
                               'Score': percIdentity, 'Strand': tempStrand}
                bedDFList.append(tempBedDict)
                record = SeqRecord(
                Seq(tempSub),
                    id=numtName,
                    description=accession + ":" + str(tempStartSub) + "-" + str(tempEndSub))
                
                #pdb.set_trace()
                recordList.append(record)
                counter += 1

                correctedStart = min(hsp['hit_from'], hsp['hit_to'])
                correctedEnd = max(hsp['hit_from'], hsp['hit_to'])
                tempHaplotypeDict = calculateHaplotypes(tempQuery, tempSub, queryStart,
                                accession, correctedStart, correctedEnd)
                haplotypeList.append(tempHaplotypeDict)
    
    with open(outPrefix +'.haplotypes.json', 'w') as hapFile:
        json.dump({'haplotypes': haplotypeList}, hapFile, indent=4)

    with open(outPrefix + ".fasta", "w") as output_handle:
        SeqIO.write(recordList, output_handle, "fasta")
    #pdb.set_trace()
    bedDF = pd.DataFrame(bedDFList)
    #correct to 0-index exclusive
    bedDF['Start'] = bedDF['Start'] - 1
    outBed = pr.PyRanges(bedDF)
    outBed.to_bed(outPrefix + ".bed")
        
