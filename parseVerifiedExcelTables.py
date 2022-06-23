#!/usr/local/bin/python3.6

import numpy as np
import pandas as pd
from checks.utils import remove_white_space

global_found_tables = []
analytes = []


def removeRepeatedValRows(df_table):
    try:
        if len(df_table.columns) > 1:
            df_table = df_table.dropna(axis=0, how='all')
            repeatedValRows = df_table.eq(df_table.iloc[:, 0], axis=0).all(axis=1)
            df_table = df_table[~repeatedValRows]
    except:
        pass

    return df_table


def parseVerifiedExcelTables(verified_tables):
    tables_with_dataframes = []
    parsed_indexes = []

    for verified_table in verified_tables:
        try:
            if not str(verified_table['table_index']).isdigit():
                table_index = str(verified_table["table_index"])
                table_title = verified_table["table_title"]
                analyte_name = verified_table["analyte_name"]
                table_type = verified_table["table_type"]
                table_subtype = verified_table.get("table_subtype", "")
                title = verified_table["tb_title"]
                if table_index not in parsed_indexes:
                    parsed_indexes.append(table_index)
                    file_location = table_index.split('*****')[0]
                    file_location = file_location.rstrip()
                    file_location = "/tmp/" + str(file_location.split('/')[-1])
                    if ".csv" in str(file_location):
                        raw_table = pd.read_csv(file_location, header=None)
                    elif ".xls" in str(file_location):
                        sheet = table_index.split('*****')[-1]
                        try:
                            excel_file = pd.ExcelFile(file_location)
                            raw_table = excel_file.parse(sheet, header=None, keep_default_na=False, skip_blank_lines=True)
                            raw_table.replace("", np.nan, inplace=True)
                            raw_table.dropna(how="all", inplace=True)
                        except Exception as e:
                            pass
                    raw_table = raw_table.applymap(lambda x: x.strip() if isinstance(x, str) else x)
                    raw_table = raw_table.reset_index(drop=True)
                    # raw_table = removeRepeatedValRows(raw_table)
                    raw_table = raw_table.fillna("")
                    raw_table = raw_table.applymap(str)
                    raw_table = remove_white_space(raw_table)
                    if raw_table.empty:
                        continue
                    else:
                        # df_json = raw_table.to_json()

                        table_info = {"table_title": table_title, "analyte_name": analyte_name, "table_type": table_type,
                                      "tb_title": title,
                                      "table_rows": raw_table, "table_subtype": table_subtype}
                        tables_with_dataframes.append(table_info)
        except:
            continue

    return tables_with_dataframes
