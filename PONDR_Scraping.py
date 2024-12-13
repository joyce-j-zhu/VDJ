import pandas as pd
import requests
from io import StringIO

INPUT_FILENAME = 'input.csv'
OUTPUT_FILENAME = 'output.csv'
COL_1 = 'V_CDR3_J_Sequence'
COL_2 = 'VL3'
COL_3 = 'VSL2'
COL_4 = 'Filename'
COL_5 = 'Read ID'

df = pd.read_csv(INPUT_FILENAME, usecols=['Filename','Read ID','V_CDR3_J_Sequence'])

#Find POST and copy as cURL and find matching python code at this website https://www.scrapingbee.com/curl-converter/python/
#EXAMPLE
headers = {
    'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10.15; rv:121.0) Gecko/20100101 Firefox/121.0',
    'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,*/*;q=0.8',
    'Accept-Language': 'en-US,en;q=0.5',
    'Referer': 'http://www.pondr.com/',
    'Content-Type': 'multipart/form-data; boundary=---------------------------13776033352742610101080519968',
    'Origin': 'http://www.pondr.com',
    'Connection': 'keep-alive',
    'Upgrade-Insecure-Requests': '1',
}

#Create empty dataframe
final_df = pd.DataFrame(columns=[COL_4, COL_5, COL_1, COL_2, COL_3])

#Iterate through each row to obtain V-CDR3-J sequence
for idx, row in df.iterrows():
    filename = row['Filename']
    samplename = row['Read ID']
    sequence = row['V_CDR3_J_Sequence']
    data = f'-----------------------------13776033352742610101080519968\r\nContent-Disposition: form-data; name="VL3"\r\n\r\non\r\n-----------------------------13776033352742610101080519968\r\nContent-Disposition: form-data; name="VSL2"\r\n\r\non\r\n-----------------------------13776033352742610101080519968\r\nContent-Disposition: form-data; name="CHStart"\r\n\r\n\r\n-----------------------------13776033352742610101080519968\r\nContent-Disposition: form-data; name="CHEnd"\r\n\r\n\r\n-----------------------------13776033352742610101080519968\r\nContent-Disposition: form-data; name="ProteinName"\r\n\r\nTest\r\n-----------------------------13776033352742610101080519968\r\nContent-Disposition: form-data; name="AccessionCode"\r\n\r\n\r\n-----------------------------13776033352742610101080519968\r\nContent-Disposition: form-data; name="Sequence"\r\n\r\n{sequence}\r\n-----------------------------13776033352742610101080519968\r\nContent-Disposition: form-data; name="wcwraw"\r\n\r\non\r\n-----------------------------13776033352742610101080519968\r\nContent-Disposition: form-data; name="submit_result"\r\n\r\nSubmit Query\r\n-----------------------------13776033352742610101080519968\r\nContent-Disposition: form-data; name=".cgifields"\r\n\r\nCH\r\n-----------------------------13776033352742610101080519968\r\nContent-Disposition: form-data; name=".cgifields"\r\n\r\nVL3\r\n-----------------------------13776033352742610101080519968\r\nContent-Disposition: form-data; name=".cgifields"\r\n\r\nVSL2\r\n-----------------------------13776033352742610101080519968\r\nContent-Disposition: form-data; name=".cgifields"\r\n\r\nXL1\r\n-----------------------------13776033352742610101080519968\r\nContent-Disposition: form-data; name=".cgifields"\r\n\r\nCDF\r\n-----------------------------13776033352742610101080519968\r\nContent-Disposition: form-data; name=".cgifields"\r\n\r\nCAN\r\n-----------------------------13776033352742610101080519968\r\nContent-Disposition: form-data; name=".cgifields"\r\n\r\nseq\r\n-----------------------------13776033352742610101080519968\r\nContent-Disposition: form-data; name=".cgifields"\r\n\r\nVLXT\r\n-----------------------------13776033352742610101080519968\r\nContent-Disposition: form-data; name=".cgifields"\r\n\r\nwcwraw\r\n-----------------------------13776033352742610101080519968\r\nContent-Disposition: form-data; name=".cgifields"\r\n\r\nstats\r\n-----------------------------13776033352742610101080519968\r\nContent-Disposition: form-data; name=".cgifields"\r\n\r\ngraphic\r\n-----------------------------13776033352742610101080519968--\r\n'
    response = requests.post('http://www.pondr.com/', headers=headers, data=data)
    cleaned = response.content.decode("utf-8")

    #find start and end index between last <PRE> and </PRE> tags
    i = len(cleaned)
    flag = False
    while(i > 0):
        if cleaned[i-4:i] == "/PRE":
            end = i-5
            flag = True
            i -= 1
            continue
        if flag and cleaned[i-3:i] == "PRE":
            start = i+1
            break
        i -= 1

    intermediate_df = pd.read_csv(StringIO(cleaned[start:end]), sep="\t")

    #Find the average of the VL3/VSL2 individual AA scores in a V-CDR3-J AA sequence
    try:
        VL3_mean = intermediate_df['VL3'].mean()
    except:
        VL3_mean = 'NA'
    try:
        VSL2_mean = intermediate_df['Unnamed: 3'].mean()
    except:
        VL3_mean = 'NA'
    #Create row of desired information
    row = pd.DataFrame({COL_4: filename, COL_5: samplename, COL_1: sequence, COL_2: VL3_mean, COL_3: VSL2_mean}, index=[0])
    #Print output for visualization
    print(row)
    #Add row to existing dataframe
    final_df = pd.concat([final_df, row], ignore_index = True)

#OUTPUT TO CSV
final_df.to_csv(OUTPUT_FILENAME, header=True, index=False)
