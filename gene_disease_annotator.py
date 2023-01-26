import pandas as pd
import numpy as np
import requests
import json
import time
import logging
import urllib3

urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)  # disables InsecureRequestWarning that appeared every time I did a request to DISGENET PLUS


def disgenet_api_request_string(items_to_query: tuple, endpoint: str, page_num: str) -> dict:
    """Gene-Disease API request to DISGENET PLUS, specifying endpoint and which page_number to be returned."""

    json_dict = {}
    params = {}
    params['page_number'] = page_num  # page_number default value is 0. Shows 100 top variant-disease (VD) associations, ordered by descending VD score
    # params['skey'] = "pkdgnpswd"
    # params['page_size'] = '10'
    # params['source'] = 'ALL'
    params['disease'] = items_to_query[0]  # test new endpoint parameter
    params['gene_ensembl_id'] = items_to_query[1]

    HTTPheadersdict = {}
    HTTPheadersdict['Authorization'] = API_KEY # disgenet user API key
    HTTPheadersdict['accept'] = 'application/json'

    ## get request using rsIDs
    response = requests.get(endpoint, params=params, headers=HTTPheadersdict, verify=False)

    ## get request using chr_position. TEST SERVER URL
    # response = requests.get("http://18.156.26.238:9007/api/vda/evidence", params=params, headers=HTTPheadersdict, verify=False)

    if not response.ok:
        try:
            if response.status_code == 429:
                while response.ok is False:
                    print(
                        "You have reached a query limit for your user. Please wait {} seconds until next query".format(
                            response.headers['x-rate-limit-retry-after-seconds']))
                    logging.warning(
                        "You have reached a query limit for your user. Please wait {} seconds until next query".format(
                            response.headers['x-rate-limit-retry-after-seconds']))
                    time.sleep(int(response.headers['x-rate-limit-retry-after-seconds']))
                    print("Your rate limit is now restored")
                    logging.info("Your rate limit is now restored")
                    # Repeat your query
                    response = requests.get(endpoint, params=params, headers=HTTPheadersdict, verify=False)
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

    # print(">>>>>>>>>> Response URL: {}".format(response.url))
    # print('>>>>>>>>>> STATUS: {}'.format(response_parsed["status"]))
    # print('>>>>>>>>>> TOTAL NUMBER OF RESULTS (rsIDs): {}'.format(response_parsed["paging"]["totalElements"]))
    # print('>>>>>>>>>> NUMBER OF RESULTS RETRIEVED BY CURRENT CALL (PAGE NUMBER {}): {}'.format(response_parsed["paging"]["currentPageNumber"], response_parsed["paging"]["totalElementsInPage"]))

    if response_parsed["status"] == 'OK' and response_parsed["paging"]["totalElements"] == 0 and \
            response_parsed["paging"]["totalElementsInPage"] == 0 and 'payload' not in response_parsed.keys():
        logging.info(">>>> Good request, but there was no payload object for query item: {}!!!".format(items_to_query))
        return None

    elif len(response_parsed['payload']) == 0:
        # logging.info(">>>> Response_parsed['payload'] is empty for query item: {}".format(items_to_query))
        return None

    JSON = response.json()
    # if endpoint.endswith('variant'): # only for variant endpoint
    #     print(JSON['payload'])
    # return JSON

    processed_json = {}
    pmid_set = set()
    disease_vocabs_set = set()
    disease_name_set = set()
    disease_umls_set = set()
    ensg_id_set = set()
    disease_type_set = set()
    score_set = set()
    gene_symbol_set = set()
    ncbi_id_set = set()

    for gda_dict in JSON['payload']:  # JSON is a list of dicts, where each dict corresponds to a variant, but these variants can be repeated
        if isinstance(gda_dict, dict) and len(gda_dict) > 0 and 'geneEnsemblIDs' in gda_dict.keys() and len(gda_dict['geneEnsemblIDs']) == 1:
            c_ensg_id = items_to_query[0]
        else:
            logging.warning("Problem with the geneEnsemblIDs. There is more than one ID: {}".format(gda_dict['geneEnsemblIDs']))
            return None

        pmid = gda_dict['pmid'] if 'pmid' in gda_dict.keys() else None
        if pmid is not None: pmid_set.add(pmid)
        disease_umls = [id for id in gda_dict['diseaseVocabularies'] if 'diseaseVocabularies' in gda_dict.keys() and id.startswith('UMLS')][0]
        disease_umls_set.add(disease_umls)
        disease_name = gda_dict['diseaseName'].upper()
        disease_name_set.add(disease_name)
        [disease_vocabs_set.add(id) for id in gda_dict['diseaseVocabularies'] if
         'diseaseVocabularies' in gda_dict.keys() and id is not None]
        num_pmids = str(len(list(pmid_set)))
        [ensg_id_set.add(id) for id in gda_dict['geneEnsemblIDs'] if 'geneEnsemblIDs' in gda_dict.keys() and id is not None]
        if 'score' in gda_dict.keys() and gda_dict['score'] is not None: score = score_set.add(str(gda_dict['score']))
        if 'diseaseType' in gda_dict.keys() and gda_dict['diseaseType'] is not None: disease_type_set.add(
            gda_dict['diseaseType'])
        if 'symbolOfGene' in gda_dict.keys() and gda_dict['symbolOfGene'] is not None: gene_symbol_set.add(
            gda_dict['symbolOfGene'])
        if 'geneNcbiID' in gda_dict.keys() and gda_dict['geneNcbiID'] is not None: ncbi_id_set.add(
            str(gda_dict['geneNcbiID']))

        if c_ensg_id not in processed_json.keys():  # if rsID not in dictionary create new entry
            processed_json[c_ensg_id] = {'disease_name': disease_name_set, 'diseaseID': disease_umls_set,
                                            'diseaseVocabularies': disease_vocabs_set, 'diseaseType': disease_type_set,
                                            'pmids': pmid_set, 'num_pmids': str(num_pmids),
                                            'geneSymbol': gene_symbol_set, 'score': score_set,
                                            'ensemblGeneIDs': ensg_id_set, 'geneNcbiID': ncbi_id_set}

        elif c_ensg_id in processed_json.keys():  # if rsID in dictionary, update it
            processed_json[c_ensg_id]['disease_name'] = processed_json[c_ensg_id]['disease_name'].union(
                disease_name_set)
            processed_json[c_ensg_id]['diseaseID'] = processed_json[c_ensg_id]['diseaseID'].union(
                disease_umls_set)
            processed_json[c_ensg_id]['diseaseType'] = processed_json[c_ensg_id]['diseaseType'].union(
                disease_type_set)
            processed_json[c_ensg_id]['diseaseVocabularies'] = processed_json[c_ensg_id][
                'diseaseVocabularies'].union(disease_vocabs_set)
            processed_json[c_ensg_id]['pmids'] = processed_json[c_ensg_id]['pmids'].union(pmid_set)
            processed_json[c_ensg_id]['num_pmids'] = str(len(list(processed_json[c_ensg_id]['pmids'])))
            processed_json[c_ensg_id]['score'] = processed_json[c_ensg_id]['score'].union(score_set)
            processed_json[c_ensg_id]['ensemblGeneIDs'] = processed_json[c_ensg_id]['ensemblGeneIDs'].union(
                ensg_id_set)
            processed_json[c_ensg_id]['geneSymbol'] = processed_json[c_ensg_id]['geneSymbol'].union(
                gene_symbol_set)
            processed_json[c_ensg_id]['geneNcbiID'] = processed_json[c_ensg_id]['geneNcbiID'].union(ncbi_id_set)

        # once previous data is pushed to dict, reset variables
        pmid_set = set()
        disease_name_set = set()
        disease_umls_set = set()
        disease_type_set = set()
        disease_vocabs_set = set()
        score_set = set()
        ensg_id_set = set()
        gene_symbol_set = set()
        ncbi_id_set = set()

    for ensg_id, value in processed_json.items():
        value['disease_name'] = "; ".join(value['disease_name']) if len(value['disease_name']) > 1 else str(
            list(value['disease_name'])[0]) if len(value['disease_name']) > 0 else 'None'
        value['diseaseID'] = "; ".join(value['diseaseID']) if len(value['diseaseID']) > 1 else str(
            list(value['diseaseID'])[0]) if len(value['diseaseID']) > 0 else 'None'
        value['diseaseType'] = "; ".join(value['diseaseType']) if len(value['diseaseType']) > 1 else str(
            list(value['diseaseType'])[0]) if len(value['diseaseType']) > 0 else 'None'
        value['diseaseVocabularies'] = "; ".join(value['diseaseVocabularies']) if len(
            value['diseaseVocabularies']) > 1 else str(list(value['diseaseVocabularies'])[0]) if len(
            value['diseaseVocabularies']) > 0 else 'None'
        pmid_lst = [str(x) for x in value['pmids']]
        value['pmids'] = "; ".join(pmid_lst) if len(pmid_lst) > 1 else str(pmid_lst[0]) if len(pmid_lst) > 0 else 'None'
        value['score'] = "; ".join(value['score']) if len(value['score']) > 1 else str(list(value['score'])[0]) if len(
            value['score']) > 0 else 'None'
        value['ensemblGeneIDs'] = "; ".join(value['ensemblGeneIDs']) if len(value['ensemblGeneIDs']) > 1 else str(
            list(value['ensemblGeneIDs'])[0]) if len(value['ensemblGeneIDs']) > 0 else 'None'
        value['geneSymbol'] = "; ".join(value['geneSymbol']) if len(value['geneSymbol']) > 1 else str(
            list(value['geneSymbol'])[0]) if len(value['geneSymbol']) > 0 else 'None'
        value['geneNcbiID'] = "; ".join(value['geneNcbiID']) if len(value['geneNcbiID']) > 1 else str(
            list(value['geneNcbiID'])[0]) if len(value['geneNcbiID']) > 0 else 'None'

    return processed_json


def append_alternative_allele(request_json, alt_allele):
    """Append alternative allele from input VCF file to JSON object returned from DISGENET request"""

    if len(request_json) > 1:  # in some cases there are several rsIDs associated to the same position
        for disease_id, value in request_json.items():
            riskAllele = value['riskAllele']
            c_disease_id = disease_id
            request_json[c_disease_id].update({'altAllele': alt_allele})

    else:  # if there is only one rsID associated to that position
        riskAllele = list(request_json.values())[0]['riskAllele']
        c_disease_id = list(request_json.keys())[0]
        request_json[c_disease_id].update({'altAllele': alt_allele})

    return request_json


def create_out_df(processed_json, out_filename):
    df = pd.DataFrame.from_dict(processed_json, orient='index',
                                columns=['disease_name', 'diseaseID', 'diseaseType', 'diseaseVocabularies', 'pmids',
                                         'num_pmids', 'ensemblGeneIDs', 'geneNcbiID', 'geneSymbol', 'score'])
    df.reset_index(level=0, inplace=True)
    df.rename(columns={'index': 'inputDiseaseID'}, inplace=True)
    df.to_csv(out_filename, sep="\t", header=True, index=False)


def multiple_input_string(disease_ids, ensg_ids):
    """ Create input string of multiple items for DISGENET PLUS API."""
    disease_str, ensg_str = '', ''
    max_length = 0

    for disease_id, ensg_id in zip(disease_ids, ensg_ids):
        if disease_id != 'None' and ensg_id != 'None':
            disease_str += str(disease_id) + ','
            ensg_str += str(ensg_id) + ','
            max_length += 1

            if max_length >= 100:
                disease_str = disease_str.rstrip(',')
                ensg_str = ensg_str.rstrip(',')

                return (disease_str, ensg_str)

    else:
        disease_str = disease_str.rstrip(',')
        ensg_str = ensg_str.rstrip(',')

        return (disease_str, ensg_str)



if __name__ == "__main__":
    # HOME
    vcf_path = "/home/ricard/PycharmProjects/vcf_annotation_tool/output_vcf/reformatted_excel.tsv"
    out_filename = "/home/ricard/PycharmProjects/vcf_annotation_tool/output_vcf/testing_gdas_singlepairs.tsv"

    # LAB
    # vcf_path = "/home/ricard/Documents/vcf_annotation_tool/liftover/lifted_50k_GCAT_clean.vcf"
    # out_filename = "/home/ricard/Documents/vcf_annotation_tool/output_vcf/lifted_50k_GCAT_clean.out"

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

        df = pd.DataFrame(new_list, columns=headers)
        df.replace('nan', np.nan, inplace=True)
        rsids = df['rsIDs']
        rsids_list = [i for i in rsids.to_list() if i is not np.nan and str(i).startswith('rs')]


    elif 'vcf_path' in dir():

        logging.basicConfig(level=logging.INFO)
        logging.info(">>>> START TIME: {}".format(start_time))
        # print(">>>> START TIME: {}".format(start_time))
        logging.info(">>>> INPUT FILE: {}".format(vcf_path))
        # print(">>>> INPUT FILE: {}".format(vcf_path))
        logging.info(">>>> OUTPUT FILE: {}".format(out_filename))
        # print(">>>> OUTPUT FILE: {}".format(out_filename))

        endpoint = "https://api.disgenetplus.com/api/v1/gda/summary"
        logging.info(">>>> ENDPOINT USED: {}".format(endpoint))
        i = 0
        final_dict = {}

        df = pd.read_csv(vcf_path, sep="\t", index_col=None)

        for tuple_disease_ensg in zip(df['diseaseId'], df['ensembl_geneId']):
            i += 1
            print("Iteration number {}".format(str(i)))

            # Create completely unique disease-ensg (1-1) pairs
            if len(tuple_disease_ensg[1].split(',')) > 1:
                break_down_pairs = [(tuple_disease_ensg[0], ensg) for ensg in tuple_disease_ensg[1].split(',')]
                for pairs in break_down_pairs:
                    request_json = disgenet_api_request_string(pairs, endpoint, '0')
                    if request_json is not None:
                        if len(request_json) > 1:  # does not happen in our test file
                            print(request_json)
                            print("JSON OBJECT HAS LENGTH GREATER THAN 1")
                            exit()
                        final_dict.update(request_json)
                    else:
                        logging.warning(
                            "No information available in DISGENET PLUS for entry: {}".format(pairs))

            else:
                request_json = disgenet_api_request_string(tuple_disease_ensg, endpoint, '0')

                if request_json is not None:
                    if len(request_json) > 1:  # does not happen in our test file
                        print(request_json)
                        print("JSON OBJECT HAS LENGTH GREATER THAN 1")
                        exit()
                    final_dict.update(request_json)
                else:
                    logging.warning("No information available in DISGENET PLUS for entry: {}".format(tuple_disease_ensg))

        # print(">>>> Creating output tsv file...")
        create_out_df(final_dict, out_filename)
        logging.info(">>>> File created at '" + str(out_filename) + "'")
        # print(">>>> File created at '" + str(out_filename) + "'")

    logging.info("--- EXECUTION TIME: {:.2f} minutes, {:.2f} seconds ---".format((time.time() - start_time) / 60, (time.time() - start_time) % 60))
    # print("--- EXECUTION TIME: {:.2f} minutes, {:.2f} seconds ---".format((time.time() - start_time) / 60, (time.time() - start_time) % 60))

