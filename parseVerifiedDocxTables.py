#!/usr/local/bin/python3.6
from checks.utils import remove_white_space
import logging
import traceback

import boto3
import numpy as np
import pandas as pd

s3 = boto3.client('s3')
global_found_tables = []
analytes = []


def get_base_path(filename):
    file_folder = '/'.join(filename.replace(' ', '').split('/')[:-1])
    return file_folder


def parseVerifiedDocxTables(file, bucket_name, verified_tables):
    analysis_type = file[2]
    analysis_file = ''.join(analysis_type.split()) + '.csv'
    base_path = get_base_path(file[0])
    # analysis_file_path = os.path.join(base_path, analysis_file)
    analysis_file_path = base_path + "/" + analysis_file

    logging.info(f'Reading file: {analysis_file_path}, from bucket: {bucket_name}')
    obj = s3.get_object(Bucket=bucket_name, Key=analysis_file_path)

    metadata = pd.read_csv(obj['Body'])
    tables_with_dataframes = []
    metadata = metadata.astype(str)
    for verified_table in verified_tables:
        table_index = int(verified_table["table_index"])
        table_title = verified_table["table_title"]
        analyte_name = verified_table["analyte_name"]
        table_type = verified_table["table_type"]
        table_subtype = verified_table.get("table_subtype", "")
        title = verified_table["tb_title"]
        try:
            csv_path = metadata[metadata['index'] == str(table_index)].iloc[0]['file_name']
            table_obj = s3.get_object(Bucket=bucket_name, Key=csv_path)
            df_table = pd.read_csv(table_obj['Body'], keep_default_na=False, skip_blank_lines=True)
            df_table.replace("", np.nan, inplace=True)
            # logging.info(f'Table Title {table_title} and df columns are  {list(df_table.columns)}')
            df_table = df_table.dropna(how='all')
            df_table = df_table.fillna("")
            df_table = df_table.reset_index(drop=True)
            df_table = df_table.applymap(str)
            df_table = remove_white_space(df_table)
            # df_json = df_table.to_json()
            # logging.info(f'Table JSON data {df_json}')

            table_info = {"table_title": table_title, "analyte_name": analyte_name, "table_type": table_type,
                          "tb_title": title,
                          "table_rows": df_table, "table_subtype": table_subtype}
            tables_with_dataframes.append(table_info)
        except Exception as e:
            tb = traceback.format_exc()
            logging.error(f'Error while reading table {table_title}. Error is  {tb}')
            continue
    return tables_with_dataframes
