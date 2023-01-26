
"""
Script that takes as input an excel/csv file with all the original info on FDA approved drugs (file coming from a paper)
and outputs a tsv file only with the DiseaseID - ENSG_IDs pairs. Where the DiseasIDs are unique and the ENSG are grouped
under the same DiseaseID.
"""

import pandas as pd
import numpy as np
import logging

logging.basicConfig(level=logging.INFO)

# HOME
input_f = "/home/ricard/Downloads/23246866.xlsx"
output_path = "/home/ricard/PycharmProjects/vcf_annotation_tool/output_vcf/reformatted_excel.tsv"

#LAB
# input_f = "/home/ricard/Downloads/23246866.xlsx"
# output_path = "/home/ricard/PycharmProjects/vcf_annotation_tool/output_vcf/reformatted_excel.tsv"


excel = pd.read_excel(input_f)
pairs_df = excel[['diseaseId', 'Human target Ids']]

ensg_unique = []
for ensg in pairs_df['Human target Ids']: # for each ENSG item (one or more than one each time) return only unique IDs, or '-' if there are none
    try:
        if not isinstance(ensg, (float, int)) and ensg is not None:
            ensg_unique.append(set(ensg.split(';')))
        else:
            ensg_unique.append('-') # if it is NaN
    except TypeError as e:
        logging.error(e, ensg)

out_df = pd.DataFrame(columns=['diseaseId', 'ensembl_geneId'])

for disease, ensg_set in zip(pairs_df['diseaseId'], ensg_unique):
    if disease is np.nan:
        disease = 'None'
    elif 'DOID' in disease:
        disease = 'DO_' + str(disease.split('_')[1])
    elif 'HP' in disease:
        disease = 'HPO_HP:' + str(disease.split('_')[1])
    elif 'Orphanet' in disease:
        disease = 'ORDO_' + str(disease.split('_')[1])
    else:
        disease = disease

    if ensg_set != '-':
        ensg_id = ",".join(ensg_set)
    else:
        ensg_id = 'None'

    new_row = {'diseaseId': disease, 'ensembl_geneId': ensg_id}
    out_df = out_df.append(new_row, ignore_index=True)


# print(out_df)
out_df.to_csv(output_path, sep="\t", index=False)

logging.info("Output file was created in {}".format(output_path))