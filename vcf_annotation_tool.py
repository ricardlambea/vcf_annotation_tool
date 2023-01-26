import pandas as pd
import numpy as np
import requests
import json
from cyvcf2 import VCF
import time
import logging
import urllib3


urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning) # disables InsecureRequestWarning that appeared every time I did a request to DISGENET PLUS

def disgenet_api_request_string(items_to_query: str, endpoint: str) -> dict:
    """
    Function that performs a query (a genomic variant) to a DISGENET API specified endpoint, and returns
    a dictionary with data related to that variant.

    :param items_to_query: string with the format chromosome:position, required to make a query to the API
    :param endpoint: Endpoint string.
    :return: python dictionary
    """

    ##### API REQUEST to DISGENET PLUS
    params = {}
    API_KEY = "1eb0b8dc-b975-4fc4-b41c-03dcb26a1117"  # ricard premium api key.     "6b2d3b03-7f44-49e1-8801-d9bbc948afe7" # janet developer api key
    params['page_number'] = "0"  # page_number default value is 0. Shows 100 top variant-disease (VD) associations, ordered by descending VD score
    # params['skey'] = "pkdgnpswd" # test new endpoint parameter. Must be commented for the 'variant' endpoint, otherwise will raise an error
    # params['page_size'] = '10'  # test new endpoint parameter. Must be commented for the 'variant' endpoint, otherwise will raise an error
    # params['source'] = 'ALL'  # test new endpoint parameter. Must be commented for the 'variant' endpoint, otherwise will raise an error
    # params['variant'] = items_to_query
    params['chromcoord'] = items_to_query

    HTTPheadersdict = {}
    HTTPheadersdict['Authorization'] = API_KEY # disgenet user API key
    HTTPheadersdict['accept'] = 'application/json'

    ## get request using rsIDs
    response = requests.get(endpoint, params = params, headers = HTTPheadersdict, verify = False)

    ## get request using chr_position. TEST SERVER URL
    # response = requests.get("http://18.156.26.238:9007/api/vda/evidence", params=params, headers=HTTPheadersdict, verify=False)

    if not response.ok:
        try:
            if response.status_code == 429:
                while response.ok is False:
                    print("You have reached a query limit for your user. Please wait {} seconds until next query".format(response.headers['x-rate-limit-retry-after-seconds']))
                    logging.warning("You have reached a query limit for your user. Please wait {} seconds until next query".format(response.headers['x-rate-limit-retry-after-seconds']))
                    time.sleep(int(response.headers['x-rate-limit-retry-after-seconds']))
                    print("Your rate limit is now restored")
                    logging.info("Your rate limit is now restored")
                    # Repeat your query
                    response = requests.get(endpoint, params=params,headers=HTTPheadersdict, verify=False)
                    if response.ok is True:
                        break
                    else:
                        continue
            elif response.status_code != 429:
                logging.warning("HTTP ERROR WITH REQUEST: {}".format(response.url))  # r.raise_for_status()
                logging.warning("HTTP STATUS CODE: {}".format(response.status_code))
                print("HTTP ERROR WITH REQUEST: {}".format(response.url))  # r.raise_for_status()
                print("HTTP STATUS CODE: {}".format(response.status_code))

                while response.ok is False:
                    logging.info("... Sleeping 20 seconds from now ...")
                    print("... Sleeping 20 seconds from now ...")
                    time.sleep(20)
                    response = requests.get(endpoint, params=params, headers=HTTPheadersdict, verify=False)
                    if response.ok is True:
                        break
                    else:
                        continue
        except Exception:
            pass

    response_parsed = json.loads(response.text)

    if response_parsed["status"] == 'OK' and response_parsed["paging"]["totalElements"] == 0 and response_parsed["paging"]["totalElementsInPage"] == 0 and 'payload' not in response_parsed.keys():
        return None

    elif len(response_parsed['payload']) == 0: # If Response_parsed['payload'] is empty for query item
        return None

    JSON = response.json()
    # if endpoint.endswith('variant'): # only for variant endpoint
    #     print(JSON['payload'])
        # return JSON

    processed_json = {}
    pmid_set = set()
    disease_umls = []
    disease_name_set = set()
    disease_umls_set = set()
    ncbi_set = set()
    disease_type_set = set()
    score_set = set()
    risk_allele_set = set()

    for variant_dict in JSON['payload']:  # JSON is a list of dicts, where each dict corresponds to a variant, but these variants can be repeated
        if isinstance(variant_dict, dict) and len(variant_dict) > 0 and 'variantStrID' in variant_dict.keys():
            c_rsid = variant_dict['variantStrID']
        else:
            logging.warning("Problem with the rsID. variant_dict['variantStrID'] is {}, and JSON has keys {}".format(variant_dict['variantStrID'], JSON['payload']))
            return None
        
        pmid = variant_dict['pmid'] if 'pmid' in variant_dict.keys() else None
        if pmid is not None: pmid_set.add(pmid)
        disease_umls = [id for id in variant_dict['diseaseVocabularies'] if 'diseaseVocabularies' in variant_dict.keys() and id.startswith('UMLS')][0]
        disease_umls_set.add(disease_umls)
        disease_name = variant_dict['diseaseName'].upper()
        disease_name_set.add(disease_name)
        num_pmids = str(len(list(pmid_set)))
        num_diseases = str(len(disease_name_set))
        if 'riskAllele' in variant_dict.keys() and variant_dict['riskAllele'] is not None: risk_allele_set.add(variant_dict['riskAllele'])
        if 'score' in variant_dict.keys(): score = score_set.add(str(variant_dict['score']))
        if 'geneNcbiIDs' in variant_dict.keys(): [ncbi_set.add(str(ncbi)) for ncbi in variant_dict['geneNcbiIDs']] # list of integers (IDs)
        chrom = variant_dict['chromosome'] if 'chromosome' in variant_dict.keys() else None
        coords = variant_dict['coord'] if 'coord' in variant_dict.keys() else None
        if 'diseaseType' in variant_dict.keys(): disease_type_set.add(variant_dict['diseaseType'])

        if c_rsid not in processed_json.keys():  # if rsID not in dictionary create new entry
            processed_json[c_rsid] = {'disease_name': disease_name_set, 'diseaseID': disease_umls_set, 'num_diseases': num_diseases, 'pmids': pmid_set, 'num_pmids': str(num_pmids), 'ncbiIDs': ncbi_set, 'chromosome': chrom, 'coordinates': coords, 'diseaseType': disease_type_set, 'score':score_set, 'riskAllele': risk_allele_set}

        elif c_rsid in processed_json.keys():  # Otherwise, update it
            processed_json[c_rsid]['disease_name'] = processed_json[c_rsid]['disease_name'].union(disease_name_set)
            processed_json[c_rsid]['diseaseID'] = processed_json[c_rsid]['diseaseID'].union(disease_umls_set)
            processed_json[c_rsid]['num_diseases'] = str(len(list(processed_json[c_rsid]['diseaseID'])))
            processed_json[c_rsid]['pmids'] = processed_json[c_rsid]['pmids'].union(pmid_set)
            processed_json[c_rsid]['num_pmids'] = str(len(list(processed_json[c_rsid]['pmids'])))
            processed_json[c_rsid]['score'] = processed_json[c_rsid]['score'].union(score_set)
            processed_json[c_rsid]['riskAllele'] = processed_json[c_rsid]['riskAllele'].union(risk_allele_set)
            processed_json[c_rsid]['ncbiIDs'] = processed_json[c_rsid]['ncbiIDs'].union(ncbi_set)
            processed_json[c_rsid]['diseaseType'] = processed_json[c_rsid]['diseaseType'].union(disease_type_set)
            if chrom != processed_json[c_rsid]['chromosome']: # these values should be the same for the same variant
                raise ValueError("THE CURRENT VALUE OF CHROMOSOME DOES NOT MATCH THE PREVIOUS VALUE FOR THE SAME rsID !!!")
            elif coords != processed_json[c_rsid]['coordinates']:
                raise ValueError("THE CURRENT VALUE OF COORDINATES DOES NOT MATCH THE PREVIOUS VALUE FOR THE SAME rsID !!!")

        # Once previous data is pushed to dict, reset variables
        pmid_set = set()
        disease_name_set = set()
        disease_umls_set = set()
        disease_type_set = set()
        ncbi_set = set()
        score_set = set()
        risk_allele_set = set()

        # ---- FOR VARIANT ENDPOINT ONLY
        # c_rsid = variant_dict['strID']
        # msc = str(variant_dict['mostSevereConsequences']) if variant_dict['mostSevereConsequences'] is not None else ''
        # threeletID = str(variant_dict['threeletterID']) if variant_dict['threeletterID'] is not None else ''
        # variant2genes = str(variant_dict['variantToGenes']) if variant_dict['variantToGenes'] is not None else ''
        # alleles = str(variant_dict['alleles']) if variant_dict['alleles'] is not None else ''
        # processed_json[c_rsid] = {'mostSevereConsequences':msc, 'threeletterID':threeletID, 'variant2Genes':variant2genes, 'alleles':alleles}


    for rsid, value in processed_json.items():
        value['disease_name'] = "; ".join(value['disease_name']) if len(value['disease_name']) > 1 else str(list(value['disease_name'])[0]) if len(value['disease_name']) > 0 else 'None'
        value['diseaseID'] = "; ".join(value['diseaseID']) if len(value['diseaseID']) > 1 else str(list(value['diseaseID'])[0]) if len(value['diseaseID']) > 0 else 'None'
        pmid_lst = [str(x) for x in value['pmids']]
        value['pmids'] = "; ".join(pmid_lst) if len(pmid_lst) > 1 else str(pmid_lst[0]) if len(pmid_lst) > 0 else 'None'
        value['diseaseType'] = "; ".join(value['diseaseType']) if len(value['diseaseType']) > 1 else str(list(value['diseaseType'])[0]) if len(value['diseaseType']) > 0 else 'None'
        value['ncbiIDs'] = "; ".join(value['ncbiIDs']) if len(value['ncbiIDs']) > 1 else str(list(value['ncbiIDs'])[0]) if len(value['ncbiIDs']) > 0 else 'None'
        value['riskAllele'] = "; ".join(value['riskAllele']) if len(value['riskAllele']) > 1 else str(list(value['riskAllele'])[0]) if len(value['riskAllele']) > 0 else 'None'
        value['score'] = "; ".join(value['score']) if len(value['score']) > 1 else str(list(value['score'])[0]) if len(value['score']) > 0 else 'None'

    return processed_json


def append_alternative_allele(request_json, alt_allele):
    """Append alternative allele from input VCF file to JSON object returned from DISGENET request"""
    if len(request_json) > 1:  # In some cases there are several rsIDs associated to the same position
        for rsid, value in request_json.items():
            riskAllele = value['riskAllele']
            c_rsid = rsid
            request_json[c_rsid].update({'altAllele': alt_allele})

    else:  # If there is only one rsID associated to that position
        riskAllele = list(request_json.values())[0]['riskAllele']
        c_rsid = list(request_json.keys())[0]
        request_json[c_rsid].update({'altAllele': alt_allele})

    return request_json


def create_out_df(processed_json, out_filename):
    df = pd.DataFrame.from_dict(processed_json, orient='index',columns=['disease_name', 'diseaseID', 'pmids', 'num_pmids', 'ncbiIDs', 'chromosome', 'coordinates', 'diseaseType','score','riskAllele','altAllele'])
    df.reset_index(level=0, inplace=True)
    df.rename(columns={'index': 'rsIDs'}, inplace=True)
    df.to_csv(out_filename, sep="\t", header=True, index=False)



if __name__ == "__main__":
    # HOME
    # maf = '/home/ricard/maf_files/5289feed-dd60-4b72-8a5d-3b27262ef9e9/368ff99a-a718-47c9-a773-7d04adcd6da9.wxs.aliquot_ensemble_masked.maf'
    # out_filename = '/home/ricard/maf_files/vda.tsv'
    #vcf_path = "/home/ricard/PycharmProjects/vcf_annotation_tool/liftover/lifted_50k_GCAT_clean.vcf"
    #out_filename = "/home/ricard/PycharmProjects/vcf_annotation_tool/output_vcf/testing_vdas.out"

    # LAB
    # maf_path = '/home/ricard/Documents/vcf_annotation_tool/maf_files/5289feed-dd60-4b72-8a5d-3b27262ef9e9/368ff99a-a718-47c9-a773-7d04adcd6da9.wxs.aliquot_ensemble_masked.maf'
    # out_filename = '/home/ricard/Documents/vcf_annotation_tool/maf_files/vda.tsv'
    vcf_path = "/home/ricard/Documents/vcf_annotation_tool/liftover/lifted_50k_GCAT_clean.vcf"
    out_filename = "/home/ricard/Documents/vcf_annotation_tool/output_vcf/testing_vdas.tsv"

    start_time = time.time()

    if 'maf_path' in dir():
        # read input maf file and discard header lines
        with open(maf_path, 'r') as fh:
            i_file = [f.strip() for f in fh.readlines()]
            i_file = i_file[7:]

        # reformat the data so a dataframe can be created
        new_list = []
        for line in i_file:
            new_line = line.split("\t")
            for x, i in enumerate(new_line):
                if i == '':
                    new_line[x] = i.replace('', 'nan')
            new_list.append(new_line)
        headers = new_list.pop(0)

        # create dataframe from data
        df = pd.DataFrame(new_list, columns=headers)
        df.replace('nan', np.nan, inplace=True)
        rsids = df['rsIDs']
        rsids_list = [i for i in rsids.to_list() if i is not np.nan and str(i).startswith('rs')]


    elif 'vcf_path' in dir():

        logging.basicConfig(level=logging.INFO)
        logging.info(">>>> START TIME: {}".format(start_time))
        print(">>>> START TIME: {}".format(start_time))
        logging.info(">>>> INPUT FILE: {}".format(vcf_path))
        print(">>>> INPUT FILE: {}".format(vcf_path))
        logging.info(">>>> OUTPUT FILE: {}".format(out_filename))
        print(">>>> OUTPUT FILE: {}".format(out_filename))
                
        endpoint = "https://api.disgenetplus.com/api/v1/vda/evidence"
        logging.info(">>>> ENDPOINT USED: {}".format(endpoint))
        i = 0
        final_dict = {}

        for variant in VCF(vcf_path):
            query_item = str(variant.CHROM).lstrip('chr') + ':' + str(variant.POS)
            alt_allele = str(variant.ALT[0])

            try:
                request_json = disgenet_api_request_string(query_item, endpoint)

                if request_json is not None:
                    request_json_mod = append_alternative_allele(request_json, alt_allele)

                    for rsid in request_json_mod:
                        if rsid in final_dict.keys():
                            logging.info("rsID {} already in output dictionary".format(rsid))
                        else:
                            logging.info("Added rsID {} to output file".format(rsid))
                            final_dict.update(request_json_mod)
                else:
                    logging.warning("No information available in DISGENET PLUS for entry: {}".format(query_item))

            except json.JSONDecodeError:
                logging.error("JSONDecodeError arised by query: {}".format(query_item))

            i += 1
            if i % 100 == 0:
                print("Iteration number {}".format(str(i)))

        print(">>>> Creating output tsv file...")
        create_out_df(final_dict, out_filename)
        print(">>>> File created at '" + str(out_filename) + "'")
        logging.info(">>>> File created at '" + str(out_filename) + "'")

    logging.info("--- EXECUTION TIME: {:.2f} minutes, {:.2f} seconds ---".format((time.time() - start_time)/60, (time.time() - start_time)%60))
    print("--- EXECUTION TIME: {:.2f} minutes, {:.2f} seconds ---".format((time.time() - start_time)/60, (time.time() - start_time)%60))
