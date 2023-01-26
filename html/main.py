import json
import time
from jinja2 import Environment, FileSystemLoader, StrictUndefined, Template
import pandas as pd
import re


"""
First this script has to be executed alone.
Then, once outputFile has been created, execute the following:
    $ python3 -m http.server -d html_report
Then open a web browser an type on the URL: 'localhost:8000', and the created page should appear.
"""


def reader_generator(file):
    while True:
        line = file.readline().strip()
        if not line:
            break
        # time.sleep(0.5)
        yield line


# home
variants_disgenet = "/home/ricard/PycharmProjects/vcf_annotation_tool/output_vcf/testing_vdas.tsv"
log_file = "/home/ricard/PycharmProjects/vcf_annotation_tool/output_vcf/vdas_logs.log"
outputFile = "index.html"

vcf = pd.read_csv(variants_disgenet, sep="\t")

counter_found_items = 0
not_found_counter = 0
error_counter = 0
with open(log_file, 'r') as i_file:
    for line in reader_generator(i_file):
        if 'No information' in line:
            not_found_counter += 1
        elif 'INPUT' in line:
            input_file = re.findall('\/.+', line)[0]
        elif 'OUTPUT' in line:
            output_file = re.findall('\/.+', line)[0]
        elif 'ENDPOINT' in line:
            endpoint = re.findall('https:\/\/.+', line)[0]
        elif 'ERROR' in line:
            error_counter += 1
        elif 'EXECUTION TIME' in line:
            matches = re.findall('\d+\.?\d*', line)
            minutes = matches[0]
            seconds = matches[1]
        elif 'Added rsID' in line:
            counter_found_items += 1

# vcf = pd.DataFrame(data=vcf[1:], columns=vcf[0])
# vcf_json = vcf.to_json(orient='records')
vcf_dict = vcf.to_dict(orient='records')

fileLoader = FileSystemLoader("templates")
env = Environment(loader=fileLoader, undefined=StrictUndefined)

rendered = env.get_template("report.html").render(vcf_dict = vcf_dict, title = "List of variants returned by DISGENET plus API", counterFoundItems = counter_found_items, notFoundCounter = not_found_counter, errorCounter = error_counter, inputFile = input_file, outputFile = output_file, endpoint = endpoint, minutes = minutes, seconds = seconds) # inside render there can be all the variables that we need. The names have to correspond to the ones present in the 'get_template()' file

with open("./html_report/{}".format(outputFile), "w") as out_f:
    out_f.write(rendered)

