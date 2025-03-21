import pandas as pd
import requests
import re
from io import StringIO

INPUT_FILENAME = 'input.csv' ###REPLACE WITH INPUT CSV FILE
OUTPUT_FILENAME = 'output.csv' ###REPLACE WITH DESIRED OUTPUT CSV FILE NAME
COL_1 = 'V_CDR3_J_Sequence'
COL_2 = 'VL3'
COL_3 = 'VSL2'
COL_4 = 'Filename'
COL_5 = 'Read ID'

df = pd.read_csv(INPUT_FILENAME, usecols=['Filename','Read ID','V_CDR3_J_Sequence'])

#Access site using Google Chrome, open console, select desired settings (statistics), and run
#Network tab: find pondr.com (type = document) (this is the main post request), right click, and copy as cuRL
#Convert to matching python code at this website https://www.scrapingbee.com/curl-converter/python/
headers = {
    'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.7',
    'Accept-Language': 'en-US,en;q=0.9,zh-CN;q=0.8,zh;q=0.7',
    'Cache-Control': 'max-age=0',
    'Connection': 'keep-alive',
    'Content-Type': 'multipart/form-data; boundary=----WebKitFormBoundarynRyrQvFDHnBYrvSw',
    'Origin': 'https://pondr.com',
    'Referer': 'https://pondr.com/',
    'Sec-Fetch-Dest': 'document',
    'Sec-Fetch-Mode': 'navigate',
    'Sec-Fetch-Site': 'same-origin',
    'Sec-Fetch-User': '?1',
    'Upgrade-Insecure-Requests': '1',
    'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/131.0.0.0 Safari/537.36',
    'dnt': '1',
    'sec-ch-ua': '"Google Chrome";v="131", "Chromium";v="131", "Not_A Brand";v="24"',
    'sec-ch-ua-mobile': '?0',
    'sec-ch-ua-platform': '"macOS"',
}

final_df = pd.DataFrame(columns=[COL_4, COL_5, COL_1, COL_2, COL_3])
count = 0

for idx, row in df.iterrows():
    filename = row['Filename']
    samplename = row['Read ID']
    sequence = row['V_CDR3_J_Sequence'].replace('*', '')
    #Insert data using {sequence}
    data = f'------WebKitFormBoundarynRyrQvFDHnBYrvSw\r\nContent-Disposition: form-data; name="VL3"\r\n\r\non\r\n------WebKitFormBoundarynRyrQvFDHnBYrvSw\r\nContent-Disposition: form-data; name="VSL2"\r\n\r\non\r\n------WebKitFormBoundarynRyrQvFDHnBYrvSw\r\nContent-Disposition: form-data; name="CHStart"\r\n\r\n\r\n------WebKitFormBoundarynRyrQvFDHnBYrvSw\r\nContent-Disposition: form-data; name="CHEnd"\r\n\r\n\r\n------WebKitFormBoundarynRyrQvFDHnBYrvSw\r\nContent-Disposition: form-data; name="ProteinName"\r\n\r\nname\r\n------WebKitFormBoundarynRyrQvFDHnBYrvSw\r\nContent-Disposition: form-data; name="AccessionCode"\r\n\r\n\r\n------WebKitFormBoundarynRyrQvFDHnBYrvSw\r\nContent-Disposition: form-data; name="Sequence"\r\n\r\n{sequence}\r\n------WebKitFormBoundarynRyrQvFDHnBYrvSw\r\nContent-Disposition: form-data; name="stats"\r\n\r\non\r\n------WebKitFormBoundarynRyrQvFDHnBYrvSw\r\nContent-Disposition: form-data; name="submit_result"\r\n\r\nSubmit Query\r\n------WebKitFormBoundarynRyrQvFDHnBYrvSw\r\nContent-Disposition: form-data; name=".cgifields"\r\n\r\nCAN\r\n------WebKitFormBoundarynRyrQvFDHnBYrvSw\r\nContent-Disposition: form-data; name=".cgifields"\r\n\r\nCDF\r\n------WebKitFormBoundarynRyrQvFDHnBYrvSw\r\nContent-Disposition: form-data; name=".cgifields"\r\n\r\nCH\r\n------WebKitFormBoundarynRyrQvFDHnBYrvSw\r\nContent-Disposition: form-data; name=".cgifields"\r\n\r\nVL3\r\n------WebKitFormBoundarynRyrQvFDHnBYrvSw\r\nContent-Disposition: form-data; name=".cgifields"\r\n\r\nVLXT\r\n------WebKitFormBoundarynRyrQvFDHnBYrvSw\r\nContent-Disposition: form-data; name=".cgifields"\r\n\r\nVSL2\r\n------WebKitFormBoundarynRyrQvFDHnBYrvSw\r\nContent-Disposition: form-data; name=".cgifields"\r\n\r\nXL1\r\n------WebKitFormBoundarynRyrQvFDHnBYrvSw\r\nContent-Disposition: form-data; name=".cgifields"\r\n\r\ngraphic\r\n------WebKitFormBoundarynRyrQvFDHnBYrvSw\r\nContent-Disposition: form-data; name=".cgifields"\r\n\r\nseq\r\n------WebKitFormBoundarynRyrQvFDHnBYrvSw\r\nContent-Disposition: form-data; name=".cgifields"\r\n\r\nstats\r\n------WebKitFormBoundarynRyrQvFDHnBYrvSw\r\nContent-Disposition: form-data; name=".cgifields"\r\n\r\nwcwraw\r\n------WebKitFormBoundarynRyrQvFDHnBYrvSw--\r\n'
    response = requests.post('https://pondr.com/', headers=headers, data=data, verify=False)
    cleaned = response.content.decode("utf-8")
    #Regex pattern to find "Average Prediction Score: " followed by a number
    pattern = r"Average Prediction Score:\s*([\d\.]+)"
    #Find all matches
    matches = re.findall(pattern, cleaned)
    #Save to variables
    try:
        VL3_mean = float(matches[0])
    except:
        VL3_mean = 'NA'
    try:
        VSL2_mean = float(matches[1])
    except:
        VSL2_mean = 'NA'

    row = pd.DataFrame({COL_4: filename, COL_5: samplename, COL_1: sequence, COL_2: VL3_mean, COL_3: VSL2_mean}, index=[0])
    print(row)
    final_df = pd.concat([final_df, row])
    if count % 100 == 0:
        final_df.to_csv(OUTPUT_FILENAME, header=True, index=False)
    count += 1
