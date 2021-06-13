import sys
from googleapiclient.discovery import build
import requests
import json
from pprint import pprint
import time
import glob
import pdb
import functools
import itertools
import pyranges as pr
import pysam
import re
from google.cloud import storage

from google.oauth2 import service_account

#credentials
credentials = service_account.Credentials.from_service_account_file(
    '/home/gtd1521_jagmail_southalabama_edu/huvecs/development/google-cloud/alignments-65005-fdd2a449ed54.json')

#could use this for web app
#credentials = GoogleCredentials.get_application_default()


scoped_credentials = credentials.with_scopes(
    ['https://www.googleapis.com/auth/cloud-platform'])
#build the service from api
service = build('lifesciences', 'v2beta', credentials=scoped_credentials)

#storage client
storageClient = storage.Client(project="alignments-65005", credentials=credentials)
        
inFileName = sys.argv[1]
with open(inFileName) as f:
    sampleDict = json.load(f)

with open("alignment.pipeline") as alignPipeline:
    align_request_body = json.load(alignPipeline)
with open("merge.pipeline") as mergePipeline:
    merge_bams_body = json.load(mergePipeline)
parent = 'projects/alignments-65005/locations/us-central1'
prefixDirBlob =  sampleDict['blob-prefix']
prefixDir = "gs://" + sampleDict['bucket'] + "/" + prefixDirBlob

environment = align_request_body['pipeline']['environment']
environment['OUTPREFIX'] = prefixDir
environment['REFERENCEDIR'] = sampleDict['reference-dir']
environment['REFERENCE'] = sampleDict['reference']


sampleList = sampleDict['sampleList']
for sample in sampleList:
    

    sampleName = sample['sample']
    print("Sample %s" % sampleName)

    # check for existing files
    skipSample = False
    #if there are already sorted bams, delete them
    blobPrefix = prefixDirBlob + "/" + sampleName
    blobs = storageClient.list_blobs("south-al-genomics", prefix=blobPrefix)
    prog = re.compile(r'.sorted.bam')
    for blob in blobs:
        result = prog.search(blob.name)
        if(result):
            print("deleting %s" % blob.name)
            blob.delete()
    blobs = storageClient.list_blobs("south-al-genomics", prefix=blobPrefix)
    prog = re.compile(r'.sorted.merged.bam')
    #see if the merged entry is already there. If so, skip this sample
    for blob in blobs:
        result = prog.search(blob.name)
        if(result):
            print("existing merged bam %s , skipping to next sample" % blob.name)
            skipSample = True

    #don't want to repeat work so skipping if we already have the merged entry
    if(skipSample == True):
        continue

    
    numberLanes = len(sample['flowcellList'])
    for flowcell in sample['flowcellList']:
        flowcellID = flowcell['id']
        fastq1 = flowcell['fastq1']
        fastq2 = flowcell['fastq2']

        print("Aligning Flowcell %s" % flowcellID)
        pipeline = align_request_body['pipeline']['environment']
        pipeline['FASTQ1'] = fastq1
        pipeline['FASTQ2'] = fastq2
        pipeline['SAMPLE'] = sampleName
        pipeline['ID'] = flowcellID
        
        pprint(pipeline)
        workRequest = service.projects().locations().pipelines().run(parent=parent, body=align_request_body)
        workResponse = workRequest.execute()

        statusParent = workResponse['name']

        print("status url %s" % statusParent)
        notDone = True
        while notDone:
            statusRequest = service.projects().locations().operations().get(name=statusParent)
            statusResponse = statusRequest.execute()

            #if("events" in statusResponse['metadata'].keys()):
            #    machineList = statusResponse['metadata']['events']
            #    pprint(machineList)
            
            if("done" in statusResponse.keys()):
                print("Done")
                notDone = False
            print("sleep 3 minutes")
            time.sleep(180)
            
    if(numberLanes > 1):
        #merge alignments and delete unmerged
        mergeEnvs = merge_bams_body['pipeline']['environment']
        mergeEnvs['BAMDIR'] = prefixDir
        mergeEnvs['SAMPLE'] = sampleName
        #pdb.set_trace()
        workRequest = service.projects().locations().pipelines().run(parent=parent, body=merge_bams_body)
        workResponse = workRequest.execute()

        statusParent = workResponse['name']
        notDone = True
        while notDone:
            statusRequest = service.projects().locations().operations().get(name=statusParent)
            statusResponse = statusRequest.execute()

            if("done" in statusResponse.keys()):
                print("Done")
                notDone = False
            print("sleep 5 minutes")
            time.sleep(300)
        
        
        blobPrefix = prefixDirBlob + "/" + sampleName
        blobs = storageClient.list_blobs("south-al-genomics", prefix=blobPrefix)
        prog = re.compile(r'.sorted.bam')
        for blob in blobs:
            result = prog.search(blob.name)
            if(result):
                blob.delete()