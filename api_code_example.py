import requests
import json
import time

# Provide your API key
API_KEY = "YOUR-API-KEY"

# Specify query parameters
params = {}
# Retrieve the first page of results
params['page_number'] = '0'
# Retrieve disease associated to variant with rsID equal to rs2500281
params['variant'] = 'rs2500281'

# Create a dictionary with the HTTP headers of your API call
HTTPheadersdict = {}
# Set the 'Authorization' HTTP header equal to API_KEY
HTTPheadersdict['Authorization'] = API_KEY
# Set the 'accept' HTTP header to specify the response format: one of 'application/json', 'application/xml', 'application/csv'
HTTPheadersdict['accept'] = 'application/json'

# Query the vda evidence endpoint
response = requests.get('https://api.disgenetplus.com/api/v1/vda/evidence', params=params, headers=HTTPheadersdict, verify=False)

# If the status code of response is 429, it means you have reached one of your query limits
# you can retrieve the time you need to wait until doing a new query in the response headers
if not response.ok:
    if response.status_code == 429:
        while response.ok is False:
            print("You have reached a query limit for your user. Please wait {} seconds until next query".format(response.headers['x-rate-limit-retry-after-seconds']))
            time.sleep(int(response.headers['x-rate-limit-retry-after-seconds']))
            print("Your rate limit is now restored")
            # Repeat your query
            response = requests.get('https://api.disgenetplus.com/api/v1/vda/evidence', params=params,headers=HTTPheadersdict, verify=False)
            if response.ok is True:
                break
            else:
                continue

# Parse the response content in JSON format since we set 'accept:application/json' as HTTP header
response_parsed = json.loads(response.text)
print('>>>>>>>>>> STATUS: {}'.format(response_parsed["status"]))
print('>>>>>>>>>> TOTAL NUMBER OF RESULTS (rsIDs): {}'.format(response_parsed["paging"]["totalElements"]))
print('>>>>>>>>>> NUMBER OF RESULTS RETRIEVED BY CURRENT CALL (PAGE NUMBER {}): {}'.format(response_parsed["paging"]["currentPageNumber"], response_parsed["paging"]["totalElementsInPage"]))






