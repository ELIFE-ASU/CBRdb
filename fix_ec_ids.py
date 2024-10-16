import os
import time
import pandas as pd
import requests
from requests import Session
from requests.adapters import HTTPAdapter
from urllib3.util import Retry


def prepare_session():
    # Make the session
    s = Session()
    # Add retries
    retries = Retry(
        total=5,
        backoff_factor=0.1,
        status_forcelist=[502, 503, 504],
        allowed_methods={'POST'},
    )
    # Mount the session
    s.mount('https://', HTTPAdapter(max_retries=retries))
    return s


def get_ec_ids(session, kegg_website=r"https://rest.kegg.jp/link/enzyme/reaction", request_sleep=0.2):
    # Get the data
    try:
        # Get the response
        response = session.get(f"{kegg_website}", timeout=10.0)
        # Limit the number of requests
        time.sleep(request_sleep)
    except requests.exceptions.RequestException as e:
        # Some error in the connection
        print(f"Error connection exception {e}", flush=True)
        exit()
    # Check if the response is ok
    if response.ok:
        # Split the response into lines
        lines = response.text.replace("rn:", "").replace("ec:", "").split("\n")[:-1]
        # Split each string by the tab character and create a list of lists
        split_data = [item.split('\t') for item in lines]
        # Create a DataFrame from the list of lists
        df = pd.DataFrame(split_data, columns=['Reaction', 'Enzyme'])
        return df

    else:
        # Some error in the response
        print(f"Error response {response.status_code}")
        exit()


if __name__ == "__main__":
    print("Program started", flush=True)
    file_ec_ids = "Data/ec_ids.csv.zip"

    if not os.path.exists(file_ec_ids):
        session = prepare_session()
        data = get_ec_ids(session)
        # write the data to a file
        data.to_csv(file_ec_ids, compression='zip', encoding='utf-8')

    # Read the data from the file
    data = pd.read_csv(file_ec_ids, compression='zip')
    print(f"Data: {data}", flush=True)

    # read the merged data
    file_merged = "Data/full_processed_merged.csv.zip"
    data_merged = pd.read_csv(file_merged, compression='zip')