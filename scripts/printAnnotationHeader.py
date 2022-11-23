import re

# POSSIBLE_TARGETS = sorted(["molMol01", "takFla02", "takRub01", "tetNig2", "synSco01", "hipCom01", "hipEre01", "cynSem", "monAlb01"])
# POSSIBLE_OUTGROUPS = sorted(["ampCit01", "oreNil02", "neoBri01", "hapBur01", "punNye01", "mayZeb03", "notFur02", "kryMar01", "ausLim01", 
#                     "xipMac02", "xipHel01", "xipCou01", "poeRet02", "poeFor01", "funHet01", "latCal01", "parOli02", "serQui01", 
#                     "serDum01", "priCar1", "gasAcu14", "punPun02", "larCro01", "miiMii01", "dicLab01", "labBer01", "bolPec01"])
# TARGET_CLADES = sorted(["clade:tetraodontiformes", "clade:syngnathids", "clade:sole", "clade:eel"])
# TARGET_CLADE_ASSIGNMENTS = {"molMol01": "tetraodontiformes", "takFla02": "tetraodontiformes", "takRub01": "tetraodontiformes", 
#                             "tetNig2": "tetraodontiformes", "synSco01": "syngnathids", "hipCom01": "syngnathids",
#                             "hipEre01": "syngnathids", "cynSem": "sole", "monAlb01": "eel"}


# header = re.sub(',', '\t', "chrom,start,end,length,ID,assc'dGene,slidingWindowSizesUsed (10|25|50|100bp)," +
#                        "1stUSTSS_name,1stUSTSS_geneID,1stUSTSS_transcriptID,1stUSTSS_description,1stUSTSS_distance,1stDSTSS_name," +
#                      "1stDSTSS_geneID,1stDSTSS_transcriptID,1stDSTSS_description,1stDSTSS_distance,assc'dDel,ixCodingExons," +
#                        "ix5UTRs,ix3UTRs,ixIntrons,ixNonCodingExons,numMappedSpecies,numTargetViolators,numTargetCladeViolators," +
#                       "numTargetsOk,numMappedTargets,numOutgroupViolators,numOutgroupsConserved,numMappedOutgroups")
header = re.sub(',', '\t', "chrom,start,end,length,ID,assc'dGene,slidingWindowSizesUsed (10|25|50|100bp)," +
                       "1stUSTSS_name,1stUSTSS_geneID,1stUSTSS_transcriptID,1stUSTSS_description,1stUSTSS_distance,1stDSTSS_name," +
                     "1stDSTSS_geneID,1stDSTSS_transcriptID,1stDSTSS_description,1stDSTSS_distance,assc'dDel,ixCodingExons," +
                       "ix5UTRs,ix3UTRs,ixIntrons,ixNonCodingExons")
# header = header + "\t" + "\t".join(POSSIBLE_TARGETS)
# header = header + "\t" + "\t".join(POSSIBLE_OUTGROUPS)
# header = header + "\t" + "\t".join(TARGET_CLADES)

print(header)
