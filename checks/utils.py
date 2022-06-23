import configparser
import datetime
import decimal
import json
import re
import statistics
from decimal import Decimal

import dateutil.parser
import numpy as np
import pandas as pd

from checks.check_result import red


def get_duplicate_columns(table):
    columns = list(table.columns)
    duplicate_column = [False]*len(columns)
    prev_col = columns[0]
    for i in range(len(columns) - 1):
        if prev_col == columns[i+1]:
            duplicate_column[i] = True

        prev_col = columns[i+1]
    return np.array(duplicate_column)


def check_orientation(df):
    column = list(df.columns)
    concentration = extract_nominal_conc(column)
    if concentration is None:
        return "VR"
    else:
        return "HZ"


def check_batches_column(test_rows):
    firsts = list(set(list(test_rows['batch_number_first'].dropna())))
    lasts = list(set(list(test_rows['batch_number_last'].dropna())))

    if len(lasts) > len(firsts):
        return test_rows['batch_number_last']
    elif len(firsts) > len(lasts):
        return test_rows['batch_number_first']
    else:
        return test_rows['batch_number_last']


def isfloat(value):
    try:
        float(value)
        return True
    except ValueError:
        return False


def parse_int_last(text):
    try:
        parsed = re.findall(r"[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", str(text))
        return abs(int(parsed[-1])) if parsed else None
    except ValueError:
        return None


def remove_duplicate_results(raw_results):
    overall_results = []

    if raw_results:
        overall_results.append(raw_results[0])

        for item in raw_results:
            num = 0
            for check_item in overall_results:
                if item.message == check_item.message and item.table_type == check_item.table_type and item.context == check_item.context:
                    num += 1
            if num == 0:
                overall_results.append(item)

    return overall_results


def get_sheet(table_type):
    tokens = table_type.split(" ")
    sheet_name = '.'.join([x[0] for x in tokens])
    return sheet_name


def write_to_excel(writer, notation_worksheet, table_obj, **kwargs):
    """
    writer: xlsx writer
    notation_worksheet: Notation sheet
    table: DataFrame
    kwargs
    count: int
    style: bold/italic
    color: blue/red
    """
    c = kwargs["count"]
    style = kwargs["style"]
    color = kwargs["color"]

    table = table_obj.table_df

    if not table.empty:
        table_type = table_obj.table_type
        table_title = table_obj.table_title

        sheet = get_sheet(table_type)
        sheet += f"_{c - 1}"

        notation_worksheet.write(f'A1', 'Sheet Name', style)
        notation_worksheet.write(f'B1', 'Table Type', style)
        notation_worksheet.write(f'C1', 'Table Title', style)
        notation_worksheet.write(f'C{c}', table_title)
        notation_worksheet.write(f'B{c}', table_type)
        notation_worksheet.write_url(f'A{c}', f"internal:'{sheet}'!A1")
        notation_worksheet.write(f'A{c}', sheet, color)
        c += 1

        table.to_excel(writer, sheet_name=sheet, index=False)
        worksheet = writer.sheets[sheet]
        (max_row, max_col) = table.shape
        column_settings = [{'header': column} for column in table.columns]
        worksheet.add_table(0, 0, max_row, max_col - 1, {'columns': column_settings})
        worksheet.set_column(0, max_col - 1, 12)

    return writer, notation_worksheet, c


def ppd_merge_split_rows(table):
    def merge_rows(_pre_row, _next_row):
        _row = list(map(lambda x, y: f"{x} {y}" if str(x).strip() != "" else y, _pre_row, _next_row))
        return _row

    values = table.values.tolist()
    column = table.columns
    pre_row = []
    df_values = []
    is_merged = False
    for row in values:
        _pre_row = list(filter(lambda x: "theoretical" in str(x).lower() or "concentration" in str(x).lower() or (
                "%" in str(x).lower() and "diff" in str(x).lower()), pre_row))
        if _pre_row:
            next_row = list(filter(
                lambda x: "concentration" in str(x).lower() or "theoretical" in str(x).lower() or "dilution" in str(
                    x).lower(), row))
            if next_row:
                merged_row = merge_rows(pre_row, row)
                df_values.append(merged_row)
                is_merged = True
        pre_row = row
        if not is_merged:
            df_values.append(row)
        else:
            is_merged = False
            if len(df_values) > 2: df_values.pop(-2)
    final_table = pd.DataFrame(data=df_values, columns=column)
    return final_table


def find_units(header: list):
    found_unit = ""
    for i, cell in enumerate(header):
        cell_string = str(cell)
        cell_array = cell_string.split(" ")
        for word in cell_array:
            if "ml" in word.lower() and "g" in word.lower() and "/" in word.lower():
                if len(word) == 5 or len(word) == 4:
                    found_unit = str(word)
                elif len(word) == 6:
                    found_unit = word.split(")")[0]
                elif len(word) == 7:
                    found_unit = word
    return found_unit


def extract_nominal_conc(raw_headers):
    nominal_conc = dict()
    c = 1
    try:
        unit = find_units(raw_headers)
        raw_headers = [str(x).split(unit)[0].strip() for x in raw_headers if len(str(x).split(unit)) == 2]
        try:
            nominal_values = [Decimal(x.split(" ")[-1].strip()) for x in raw_headers]
        except (decimal.InvalidOperation, IndexError):
            nominal_values = [Decimal(x.split("(")[-1].strip()) for x in raw_headers]
        for i, val in enumerate(nominal_values):
            search_header_name = raw_headers[i].lower()
            val = round(float(val), 3)
            if "lloq" in search_header_name:
                nominal_conc["lloq_conc"] = val
            elif "uhqc" in search_header_name:
                nominal_conc["uhqc_conc"] = val
            elif "hqc" in search_header_name:
                nominal_conc["hqc_conc"] = val
            else:
                nominal_conc[f"conc_{c}"] = val
                c += 1
    except:
        return None
    return nominal_conc


def read_json_object(table_type):
    result_dict = dict()
    json_file_name = '/tmp/ObjectNotation.json'
    with open(json_file_name, 'r') as r:
        table_dict = json.load(r)

    for object_dict in table_dict:
        if object_dict.get('type') == table_type:
            result_dict = object_dict
            break
    return result_dict


def generalize_column_names(col_names, dictionary_name):
    found_cols = []
    missing_cols = []
    for col_name in col_names:
        match = False
        for i in dictionary_name['columns']:
            values = [str(x).lower().strip() for x in i["value"]]
            c = str(col_name).lower().strip().replace("\n", "")
            if str(col_name).lower().strip().replace("\n", "") in values:
                col_name = i['name']
                found_cols.append(col_name)
                match = True
                break
        if not match:
            found_cols.append(col_name)
    columns = dictionary_name['columns']
    missing_cols_dict = [cols for cols in columns if cols["name"] not in found_cols]

    column_iter = iter(missing_cols_dict)
    count = 0
    while True:
        try:
            name = next(column_iter)["name"]
            if name == "mean" or name == "cv" or name == "re" or name == "sd" or name == "n" or name == "te" or name == "ar" or name == "cycle_id":
                del missing_cols_dict[count]
                column_iter = iter(missing_cols_dict)
                count = 0
                continue
            else:
                count += 1
        except StopIteration:
            break

    for i in missing_cols_dict:
        missing_cols.append(i["desc"])
    return found_cols, missing_cols


def round_df_values(df, digit):
    df_list = df.values.tolist()
    process_list = []
    for data in df_list:
        data = [f"%.{digit}f" % (round(float(x), digit)) if str(x).isdigit() else x for x in data]
        process_list.append(data)
    column = df.columns
    process_df = pd.DataFrame(data=process_list, columns=column)
    return process_df


def format_table_data_frame(table: pd.DataFrame, table_type: str):
    table_dict = read_json_object(table_type)
    table_cols, missing_cols = generalize_column_names(table.columns, table_dict)
    table.columns = table_cols
    return table, missing_cols


def remove_header(df: pd.DataFrame, template_type=None) -> pd.DataFrame:
    df = remove_numeric_header(df)
    if str(template_type) != "ppd":
        df = remove_duplicate_columns(df)
    df = df.reset_index(drop=True)
    return df


def remove_numeric_header(df: pd.DataFrame) -> pd.DataFrame:
    for i in range(2):
        found_digit = False
        df_column = df.columns
        for col in df_column:
            if str(col).replace(".", "").isdigit():
                found_digit = True
                break
        if found_digit:
            new_header = df.iloc[0]
            df = df[1:].reset_index(drop=True)
            df.columns = new_header
    return df


def get_analyte(analytes: list, table_title: str):
    table_title = table_title.lower()
    for analyte in analytes:
        if analyte.lower() in table_title:
            return analyte
    return None


def check_split_batch_analyte(analyte: str):
    analyte = analyte.split("/")
    return True if len(analyte) > 1 else False


def split_batch_analyte(analyte: str):
    analyte = analyte.split("/")
    analyte = [x.strip() for x in analyte]
    return analyte


def remove_footer_sb_table(df):
    row, cols = df.shape
    if row >= 3:
        for i in range(row - 1, 0, -1):
            df_row = df.loc[i]
            df_row = [str(x).strip() for x in df_row]
            first_cell = df_row[0]
            nan_count = df_row.count("")
            if cols <= 7:
                if nan_count > (cols - 2) and first_cell != "":

                    df = df[:-1].reset_index(drop=True)

                else:
                    break
            elif cols > 7:
                if nan_count > (cols // 2 + 2) and first_cell != "":

                    df = df[:-1].reset_index(drop=True)
                else:
                    break
    return df


def remove_footer(df):
    row, cols = df.shape
    if row >= 3:
        for i in range(row - 1, 0, -1):
            df_row = df.loc[i]
            df_row = [str(x).strip() for x in df_row]
            first_cell = df_row[0]
            nan_count = df_row.count("")
            if cols <= 7:
                if nan_count > (cols - 2):

                    df = df[:-1].reset_index(drop=True)

                else:
                    break
            elif cols > 7:
                if nan_count >= (cols // 2 + 2):

                    df = df[:-1].reset_index(drop=True)
                else:
                    break
    return df


def validate_unit(header, table_title, table_type):
    found_unit = []
    result = []
    for i, cell in enumerate(header):
        cell_string = str(cell)
        cell_array = cell_string.split(" ")

        for word in cell_array:
            if "ml" in word.lower() and "g" in word.lower() and "/" in word.lower():
                if len(word) == 5 or len(word) == 4:
                    found_unit.append(str(word))
                elif len(word) == 6:
                    found_unit.append(word.split(")")[0])
                elif len(word) == 7:
                    found_unit.append(word)
    if len(set(found_unit)) > 1:
        result.append(
            red(f"Unable to process table due to inconsistent concentration units in table column headers", None,
                table_title, table_type))

    return result


def find_nominal_conc(raw_headers):
    nominal_conc = dict()
    try:
        unit = find_units(raw_headers)
        if unit:
            column = list(filter(lambda x: unit in str(x), raw_headers))
            raw_headers = [str(x).split(unit)[0].strip() for x in raw_headers if len(str(x).split(unit)) == 2]

            try:
                nominal_values = [Decimal(x.split(" ")[-1].strip().replace(",", "")) for x in raw_headers]
            except decimal.InvalidOperation:
                nominal_values = [Decimal(x.split("(")[-1].strip()) for x in raw_headers]

            for i, val in enumerate(nominal_values):
                col_prefix = column[i]
                nominal_conc[col_prefix] = val

    except decimal.InvalidOperation as e:
        print("Error", e)
        return None
    return nominal_conc


def contains_split_col(x):
    return "mean" in str(x).lower() or "average" in str(x).lower() or str(x).lower() == "n"


def remove_last_element(values: list):
    if len(values) >= 2:
        if str(values[0]).lower() != "n":
            values = [x for x in values if str(x).lower() != 'n']
        else:
            values.pop(1)
    return values


def split_static_table(table_df, round="one"):
    data = pd.DataFrame()
    try:
        column_no = 0
        column = list(table_df.columns)
        column_values = table_df[column[0]].to_list()
        column_values = list(filter(contains_split_col, column_values))
        if not column_values:
            column_values = table_df[column[1]].to_list()
            column_values = list(filter(contains_split_col, column_values))
            column_no = 1
        if not column_values:
            column_values = table_df[column[2]].to_list()
            column_values = list(filter(contains_split_col, column_values))
            column_no = 2
        if not column_values:
            column_values = table_df[column[3]].to_list()
            column_values = list(filter(contains_split_col, column_values))
            column_no = 3
        column_values = remove_last_element(column_values)
        if round == "one":
            value = column_values[0]
        elif round == "two":
            if len(column_values) >= 2:
                value = column_values[-1]
            else:
                return table_df, data

        split_index = table_df[table_df[column[column_no]] == value].index

        if round == "one":
            split_index = split_index[0]
        elif round == "two":
            split_index = split_index[-1]
        data = table_df[0:split_index].reset_index(drop=True)
        static = table_df[split_index:].reset_index(drop=True)

    except Exception as e:
        return table_df, data

    return data, static


def split_ppd_selectivity(table_df):

    def ppd_contains_split_col(x):
        return "result" in str(x).lower() or "n" == str(x).lower() or ("%" in str(x).lower() and "difference" in str(x).lower())

    data_df = pd.DataFrame()
    static_df = pd.DataFrame()
    column = list(table_df.columns)
    column_values = table_df[column[0]].to_list()
    column_values = list(filter(ppd_contains_split_col, column_values))
    value = column_values[0] if column_values else None
    if value:
        split_index = table_df[table_df[column[0]] == value].index
        split_index = split_index[0]
        data_df = table_df[0:split_index].reset_index(drop=True)
        static_df = table_df[split_index:].reset_index(drop=True)

    return data_df, static_df


def get_num(index, run_id):
    run_id = run_id[index:]
    num = re.findall("[0-9]+", run_id)
    return abs(int(num[0])) if num else None


def get_unique_ids(run_ids):
    difference = lambda x, y: list(set(x).difference(set(y)))
    unique_val = [difference(run_ids[i], run_ids[i + 1]) for i in range(len(run_ids) - 1)]
    filter_val = [x for y in unique_val for x in y if str(x).isdigit()]
    if len(filter_val) and str(filter_val[0]).isdigit():
        idx = [run_ids[i].index(unique_val[i][0]) for i in range(len(run_ids) - 1) if len(unique_val[i])]
        idx = min(idx)
        unique_ids = [get_num(idx, run_id) for run_id in run_ids]
        return unique_ids
    else:
        for x in range(0, len(run_ids)):
            try:
                run_ids[x] = abs(int(re.findall(r"[-+]?\d*\.*\d+", run_ids[x])[-1]))
            except IndexError:
                run_ids[x] = ""
        return run_ids


def get_cycles(header):
    cycles = [str(x).split(" ")[1] for x in header if len(str(x).split(" ")) >= 2]
    cycles = list(filter(lambda x: "x" in str(x).lower(), cycles))
    cycles = set(cycles)
    if cycles:
        return cycles
    else:
        return None


def parse_decimal_2(text: str):
    try:
        parsed = re.findall(r"[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", text)
        if parsed:
            return decimal.Decimal(parsed[0])
        else:
            pattern = re.compile(r"aql|bql|blq|alq|<lloq|>uloq")
            return text if pattern.search(text.lower()) else None
    except decimal.InvalidOperation:
        return None


def parse_decimal_2_1(text):
    if isinstance(text, decimal.Decimal):
        return text
    else:
        try:
            parsed = re.findall(r"[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", str(text))
            return decimal.Decimal(parsed[0]) if parsed else None
        except decimal.InvalidOperation:
            return None


def parse_date(str_date):
    try:
        if isinstance(str_date, datetime.date):
            return str_date
        else:
            return dateutil.parser.parse(str_date).date()
    except Exception as e:
        return None


def next_table_title(tables, index):
    try:
        table_title = tables[index]["table_title"]
        return table_title
    except IndexError:
        return ""


def parse_decimal(text):
    try:
        text = re.sub(",", "", str(text))
        parsed = re.findall(r"[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", str(text))
        return decimal.Decimal(str(parsed[0]).replace(",", "")) if parsed else None
    except decimal.InvalidOperation:
        return None


def parse_last_decimal(text):
    try:
        parsed = re.findall(r"[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", str(text))
        return decimal.Decimal(str(parsed[-1]).replace(",", "")) if parsed else None
    except decimal.InvalidOperation:
        return None


def parse_signed_int(text):
    try:
        parsed = re.findall(r"[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", str(text))
        return int(float(str(parsed[0]).replace(",", ""))) if parsed else None
    except ValueError:
        return None


def check_stability_orientation(df):
    column = list(df.columns)
    cycles = get_cycles(column)
    if cycles is None:
        return "VR"
    else:
        return "HZ"


def filter_conc_col(columns, conc_unit):
    return [x for x in columns if conc_unit in x]


def ppd_merge_split_header(table):
    prev_col = []
    table = remove_numeric_header(table)
    table = remove_duplicate_columns(table)
    for index, row in table.iterrows():
        columns = list(table.columns)
        row = row.values.tolist()
        if str(columns[0]).strip() == "" and str(row[0].strip()) != "":
            prev_col = columns.copy()
            columns = [f"{x} {y}" for x, y in zip(columns, row)]
            table.columns = columns
            table = table[1:].reset_index(drop=True)
            break
        elif str(columns[0]).strip() == "":
            table.columns = row
            table = table[1:].reset_index(drop=True)
            continue
        else:
            break
    return table, prev_col


def remove_duplicate_columns(table: pd.DataFrame) -> pd.DataFrame:
    for _ in range(table.shape[0]):
        tab_col = table.columns
        columns = {col for col in tab_col if str(col).strip() != ""}
        if len(columns) == 1 and len(columns) != len(tab_col):
            table.columns = table.iloc[0]
            table = table[1:].reset_index(drop=True)
        else:
            break
    return table


def ppd_merge_concentration_row_to_header(table):
    def merge(row_val, header_val):
        val = header_val
        if unit in header_val:
            val = header_val.replace(unit, f"{row_val} {unit}")
        return val

    columns = list(table.columns)
    unit = find_units(columns)
    column = columns[0]
    if unit and str(column).strip():
        table = ppd_merge_split_rows(table)
        rows = ['concentration after dilution', 'theoretical concentration']
        for i, row in enumerate(rows):
            columns = list(table.columns)
            nominal_table = table[table[column].str.strip().str.lower() == row].reset_index(drop=True)
            if not nominal_table.empty:
                nominal_values = nominal_table.iloc[0].to_list()
                final_columns = list(map(merge, nominal_values, columns))
                table = table[~(table[column].str.strip().str.lower() == row)].reset_index(drop=True)
                table.columns = final_columns

    return table


def format_df_calc_stats(df):
    column = list(df.columns)
    found_re = find_re_column(df)
    found_cycle = find_cycle_column(df)
    final_df = pd.DataFrame()
    if found_cycle:
        cycle_column = list(filter(lambda x: "cycle" in str(x), column))
    count = 0
    run_column = list(filter(lambda x: "run_id" in str(x), column))
    date_column = list(filter(lambda x: "assay_date" in str(x), column))
    conc = find_nominal_conc(column)
    sub_static_df = pd.DataFrame()
    for i in conc:
        if found_re:
            sub_df = pd.DataFrame(columns=['assay_date', 'run_id', 'column', 'nominal', 'conc', 're'])
        else:
            sub_df = pd.DataFrame(columns=['assay_date', 'run_id', 'column', 'nominal', 'conc'])
        sub_df['run_id'] = df[run_column[0]]
        sub_df['nominal'] = conc[i]
        sub_df['conc'] = df[i].apply(parse_decimal)
        sub_df['column'] = i
        if date_column:
            sub_df['assay_date'] = df[date_column[0]]

        if found_cycle:
            sub_df.insert(loc=0, column="cycle_id", value=df[cycle_column[0]])

        if found_re:
            column_index = df.columns.get_loc(i)
            if str(column[column_index + 1]).strip() == "re" or "re_" in str(column[column_index + 1]).strip():
                re_col = f"re_{count}"
                df.columns.values[column_index + 1] = re_col
                try:
                    sub_df["re"] = df[re_col].apply(parse_decimal)
                except KeyError:
                    sub_df["re"] = df["re"].apply(parse_decimal)
            else:
                sub_df["re"] = ""

        values = sub_df["conc"].to_list()
        sub_static_df = pd.concat([sub_static_df, build_static_df(values, conc[i])]).reset_index(drop=True)
        final_df = pd.concat([final_df, sub_df]).reset_index(drop=True)

        count += 1
    return final_df, sub_static_df, found_re, found_cycle


def format_calibration_table(df, nominal_conc):
    re_cols = list(filter(lambda x: "re" in str(x), df.columns))
    found_re = find_re_column(df)
    final_df = pd.DataFrame()
    count = 0
    column = list(df.columns)
    run_column = list(filter(lambda x: "run_id" in str(x), column))
    date_column = list(filter(lambda x: "assay_date" in str(x), column))
    if date_column:
        date_found = True
    else:
        date_found = False

    for i in nominal_conc:
        if found_re:
            sub_df = pd.DataFrame(columns=['assay_date', 'run_id', 'sample', 'column', 'nominal', 'conc', 're'])
        else:
            sub_df = pd.DataFrame(columns=['assay_date', 'run_id', 'sample', 'column', 'nominal', 'conc'])
        if date_found:
            sub_df['assay_date'] = df[date_column[0]]
        sub_df['run_id'] = df[run_column[0]]
        sub_df['sample'] = i
        sub_df['nominal'] = nominal_conc[i]
        sub_df['conc'] = df[i]
        sub_df['column'] = i
        if found_re:
            if len(re_cols) > 1:
                column_index = df.columns.get_loc(i)
                df.columns.values[column_index + 1] = f"re_{count}"
                sub_df["re"] = df[f"re_{count}"].apply(parse_decimal)
            else:
                sub_df["re"] = df["re"].apply(parse_decimal)
        final_df = pd.concat([final_df, sub_df]).reset_index(drop=True)
        count += 1

    return final_df, found_re


def fill_val(df: pd.DataFrame) -> pd.DataFrame:
    columns = ['run_id', 'cycle_id', 'assay_date', 'mean_analyte_area', 'mean_is_area', 'assay_date', 'qc_level',
               'mean_solution_peak_area', 'time', 'std_analyte_area', 'batch_qc_is_mean_area']
    for col in columns:
        try:
            df[col] = df[col].replace([''], np.nan)
        except KeyError:
            continue
    df.ffill(inplace=True)
    return df


def merge_tables(annotated_tables, unique_ids):
    concat_tables = []
    for run_id in unique_ids:
        run_tables = []
        r_id = ""
        for tables in annotated_tables:
            if run_id == tables["run"]:
                run_tables.append(tables["table"])
                r_id = tables["run"]
        df = pd.concat(run_tables).reset_index(drop=True)
        concat_tables.append({"table": df, "run": r_id})

    return concat_tables


def string_formatter(data):
    ls = []
    for val in data:
        if "&" in str(val):
            val = val.replace("&", "&amp;")
            ls.append(val)
        elif "datetime.datetime" in str(type(val)):
            ls.append(val.strftime('%d-%b-%Y'))

        elif "timestamps.Timestamp" in str(type(val)):
            val = val.to_pydatetime()
            ls.append(val.strftime('%d-%b-%Y'))

        elif issubclass(type(val), type(pd.NaT)):
            ls.append("")
        else:
            ls.append(val)
    return ls


def remove_duplicate(data):
    new_data = []
    for val in data:
        if val not in new_data:
            new_data.append(val)
    return new_data


def split_run_df(df):
    all_tables = []
    duplicate_id = []
    df_column = list(df.columns)
    run_column = list(filter(lambda x: "cycle" in str(x).lower(), df_column))
    if not run_column:
        run_column = list(filter(lambda x: "run_id" in str(x).lower(), df_column))

    try:
        ids = df[run_column[0]].to_list()
        ids = [x for x in ids if not re.sub(r"[-%\s()/_.]", "", str(x)).isalpha()]
    except Exception as e:
        all_tables.append({"table": df, "run": "", "date": ""})
        return all_tables, duplicate_id
    ids = [x for x in ids if str(x).strip() != ""]
    unique_ids = remove_duplicate(ids)
    duplicate_id = set([x for x in ids if ids.count(x) > 1])

    c = 0
    if len(ids) > 1:
        for i, val in enumerate(ids):
            if i != 0:
                split_index = df[df[run_column[0]] == val].index[0]
            if i != 0 and i == 1:
                data = df[0:split_index].reset_index(drop=True)
                all_tables.append({"table": data, "run": ids[c]})
                df = df[split_index:].reset_index(drop=True)
                c += 1
            if i != 1 and i != 0:
                data = df[0:split_index].reset_index(drop=True)
                df = df[split_index:].reset_index(drop=True)
                all_tables.append({"table": data, "run": ids[c]})
                c += 1
        all_tables.append({"table": df, "run": ids[c]})
        c += 1
    else:
        all_tables.append({"table": df, "run": ids[c]})
    all_tables = merge_tables(all_tables, unique_ids)
    return all_tables, duplicate_id


def build_static_df(values, nominal_value):
    values = [x for x in values if x is not None]
    sub_static_df = pd.DataFrame(columns=["calc_mean", "calc_sd", "calc_cv", "calc_re", "calc_n", "nominal"])
    try:
        calc_mean = statistics.mean(values)
    except (statistics.StatisticsError, TypeError):
        calc_mean = Decimal(0)

    try:
        calc_re = 100 * (calc_mean - nominal_value) / nominal_value if calc_mean else Decimal(0)
    except (decimal.DivisionByZero, TypeError):
        calc_re = Decimal(0)

    try:
        std_dev = statistics.stdev(values)
    except (statistics.StatisticsError, TypeError):
        std_dev = Decimal(0)

    try:
        calc_cv = 100 * (std_dev / calc_mean) if calc_mean else Decimal(0)
    except decimal.DivisionByZero:
        calc_cv = Decimal(0)

    try:
        sub_static_df["calc_mean"] = [calc_mean]
        sub_static_df["calc_sd"] = std_dev
        sub_static_df["calc_n"] = Decimal(len(values))
        sub_static_df["calc_re"] = calc_re
        sub_static_df["calc_cv"] = calc_cv
        sub_static_df["nominal"] = nominal_value
    except Exception as e:
        pass

    return sub_static_df


def find_re_column(df: pd.DataFrame):
    column = df.columns
    found = False
    for col in column:
        if str(col).strip() == "re":
            found = True
            break
        elif str(col).strip() == "%RE":
            found = True
            break
        elif str(col).strip() == "% RE":
            found = True
            break
        elif "re_" in str(col).lower():
            found = True
            break
    return found


def find_te_column(df):
    column = df.columns
    found = False
    for col in column:
        if "te" in str(col):
            found = True
            break
        elif "%TE" in str(col):
            found = True
            break
        elif "% TE" in str(col):
            found = True
            break
    return found


def find_ar_column(df):
    column = df.columns
    found = False
    for col in column:
        if "ar" in str(col):
            found = True
            break
        elif "%AR" in str(col):
            found = True
            break
        elif "% AR" in str(col):
            found = True
            break
    return found


def find_cycle_column(df: pd.DataFrame):
    column = df.columns
    found = False
    for col in column:
        if "cycle_id" in str(col):
            found = True
            break
        elif "cycle" in str(col):
            found = True
            break
    return found


def concat_static_df(static_df, calc_static_df):
    static_df, found_sd = convert_static_df_to_decimal(static_df)
    static_df = pd.concat([static_df, calc_static_df], axis=1)
    column = ["mean", "sd", "cv", "re", "n"]
    for col in column:
        try:
            static_df[f"per_diff_{col}"] = calculate_per_diff(static_df[col], static_df[f"calc_{col}"])
        except KeyError:
            continue
    return static_df.fillna(""), found_sd


def parse_float(text):
    try:
        parsed = re.findall(r"[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", str(text))
        return float(parsed[0]) if parsed else None
    except ValueError:
        return None


def parse_str_float(text):
    if str(text).isalpha():
        return text
    else:
        try:
            parsed = re.findall(r"[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", str(text))
            return float(parsed[0]) if parsed else None
        except ValueError:
            return None


def convert_static_df_to_decimal(static_df):
    found_sd = True
    column = ['mean', 'sd', 'cv', 're', 'n']
    for col in column:
        try:
            static_df[col] = static_df[col].apply(parse_decimal)
        except KeyError:
            if col == 'sd':
                found_sd = False
            continue
    return static_df, found_sd


def calculate_per_diff(reported_values, calculated_values):
    per_diff_res = list()
    for i, val in enumerate(calculated_values.fillna("")):
        try:
            if val != Decimal(0):
                per_diff_res.append(
                    Decimal(round(abs((round(float(val), 5) - float(reported_values[i])) / round(float(val), 5)), 3)))
            elif reported_values[i] is not None or str(reported_values[i]) != "" or reported_values[i] != Decimal(0):
                per_diff_res.append(Decimal(abs((float(reported_values[i]) - float(val)) / float(reported_values[i]))))
            else:
                per_diff_res.append(Decimal(0))
        except:
            per_diff_res.append(Decimal(0))

    return per_diff_res


def drop_blank_col(df):
    try:
        df["assay_date"] = df["assay_date"].replace("", np.nan)
    except KeyError:
        pass
    try:
        df["run_id"] = df["run_id"].replace("", np.nan)
    except KeyError:
        pass
    try:
        df["cycle_id"] = df["cycle_id"].replace("", np.nan)
    except KeyError:
        pass

    re_column = list(df.columns)
    re_column = list(filter(lambda x: "re_" in str(x) or str(x).strip() == "re", re_column))
    if re_column:
        for cols in re_column:
            try:
                df[cols] = df[cols].replace("", np.nan)
                df[cols] = df[cols].replace("NA", np.nan)
            except KeyError:
                df["re"] = df[cols].replace("", np.nan)
                df["re"] = df[cols].replace("NA", np.nan)

    df = df.dropna(axis=1, how='all')
    return df


def get_stability_type(table_title: str):
    table_title = table_title.lower()
    if ("freeze" in table_title) and ("thaw" in table_title) or ("f/t" in table_title) or (
            "freeze/thaw" in table_title):
        return "FT"
    elif "frozen" in table_title:
        return "FR"
    elif "long" in table_title and "term" in table_title or "long-term" in table_title:
        return "LT"
    elif "hemolysis" in table_title:
        return "HM"
    elif "hyperlipidemic" in table_title:
        return "HP"
    elif "processed" in table_title or "extract" in table_title:
        return "PS"
    else:
        return None


def split_time_point(text):
    text = str(text)
    text = text.replace(":", "")
    values = re.findall(r"[a-zA-Z]+.?\d+|\d+|\d?\d?\d?.?[a-zA-Z]+.?[a-zA-Z]+", text)
    values = [str(x).strip() for x in values]
    if len(values) > 2:
        values[1] = " ".join(values[1:])
        del values[2:]
    return values


def format_peak_area_table(table):
    table = fill_val(table)
    try:
        table["analyte_peak_area"] = table["analyte_peak_area"].apply(parse_decimal)
        table["mean_analyte_area"] = table["mean_analyte_area"].apply(parse_decimal)
        table["is_peak_area"] = table["is_peak_area"].apply(parse_decimal)
        table["mean_is_area"] = table["mean_is_area"].apply(parse_decimal)
        table["per_analyte_response"] = table["per_analyte_response"].apply(parse_decimal)
        table["per_is_response"] = table["per_is_response"].apply(parse_decimal)
        table["calc_per_analyte_response"] = (table["analyte_peak_area"] / table["mean_analyte_area"]) * 100
        table["calc_per_is_response"] = (table["is_peak_area"] / table["mean_is_area"]) * 100
    except:
        raise Exception
    return table


def parse_int(text):
    if isinstance(text, float) or isinstance(text, int):
        try:
            return int(text)
        except ValueError:
            return None
    else:
        try:
            parsed = re.findall(r"[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", str(text))
            return abs(int(parsed[0])) if parsed else None
        except ValueError:
            return None


def split_qc_df(df):
    all_tables = []
    df_column = list(df.columns)
    run_column = list(filter(lambda x: "level" in str(x).lower(), df_column))
    if not run_column:
        run_column = list(filter(lambda x: "concentration" in str(x).lower(), df_column))

    try:
        ids = df[run_column[0]].to_list()
        ids = [x for x in ids if " qc" in str(x).lower()]
    except IndexError:
        all_tables.append({"table": df, "level": ""})
        return all_tables
    ids = [x for x in ids if str(x).strip() != ""]
    unique_ids = remove_duplicate(ids)

    c = 0
    if len(ids) > 1:
        for i, val in enumerate(ids):
            if i != 0:
                split_index = df[df[run_column[0]] == val].index[0]
            if i != 0 and i == 1:
                data = df[0:split_index].reset_index(drop=True)
                all_tables.append({"table": data, "level": ids[c]})
                df = df[split_index:].reset_index(drop=True)
                c += 1
            if i != 1 and i != 0:
                data = df[0:split_index].reset_index(drop=True)
                df = df[split_index:].reset_index(drop=True)
                all_tables.append({"table": data, "level": ids[c]})
                c += 1
        all_tables.append({"table": df, "level": ids[c]})
        c += 1
    else:
        all_tables.append({"table": df, "level": ids[c]})
    return all_tables


def format_date(ref_date):
    # Covert date YYY-MM-DD into DD-MM-YYY
    try:
        ref_date = datetime.datetime.strptime(str(ref_date), "%Y-%m-%d").strftime("%d-%b-%Y")
        return ref_date
    except:
        return None


def validate_batch_run(table, table_title, table_type):
    result = []
    run_ids = table["run_id"].to_list()
    ids = [x for x in run_ids if not re.sub(r"[-%\s()/_]", "", x).isalpha()]
    ids = sorted([parse_int(x) for x in ids if str(x).strip() != ""])
    ids = [str(x) for x in ids]
    previous_id = ids[0]
    previous_id_index = run_ids.index(previous_id)
    for r_id in ids:
        try:
            r_id_index = run_ids.index(r_id, previous_id_index + 1)
        except ValueError:
            continue
        if previous_id == r_id:
            if previous_id_index + 1 == r_id_index:
                pass
            else:
                result.append(red(f"Unable to process table due to duplicate Batch Run ID {r_id} for more than one set"
                                  f" of result", None, table_title, table_type))
                break
        else:
            previous_id = r_id
            previous_id_index = r_id_index
    return result


def format_value(calc_value, decimal_point=2):
    try:
        return f"%.{decimal_point}f" % (
            round(float(calc_value), decimal_point)) if calc_value is not None else calc_value
    except ValueError:
        return calc_value


class FlagProperties:
    def __init__(self, file_path):
        self.file_path = file_path
        self.msg_dict = dict()

    def get_table_dict(self, table_type):
        """
        Read the file passed as parameter as a properties file.
        """
        config = configparser.RawConfigParser()
        config.read(self.file_path)
        data = config[table_type]
        self.msg_dict = data

    def get_message(self, key, **kwargs):
        """
        **kwargs
        batch: batch number
        cycle: cycle number
        reported_value: table value
        calc_value: calculated value
        threshold: test threshold
        added_str: string for batch failed/passed
        column: table column
        """
        kwargs = [('{' + key + '}', kwargs[key]) for key in kwargs]
        print(*kwargs)
        message = self.msg_dict[key]
        try:
            for item in kwargs:
                message = message.replace(item[0], str(item[1]))
        except Exception as e:
            pass

        return message


def process_static_df(static_df: pd.DataFrame, table_type: str):
    static_df = static_df.fillna("")
    static_df = drop_blank_col(static_df)
    static_df = static_df.T
    index_column = list(static_df.index)
    index_column = list(filter(lambda x: find_units(index_column) in str(x), index_column))
    static_df = remove_header(static_df)
    static_df, missing_cols = format_table_data_frame(static_df, table_type)
    static_df.insert(loc=0, column="column", value=index_column)
    static_df = static_df.reset_index(drop=True)
    return static_df


def find_cv_column(df):
    column = df.columns
    if "cv" in column:
        return True
    else:
        return False


def find_sd_column(df):
    column = df.columns
    if "sd" in column:
        return True
    else:
        return False


def find_count_column(df):
    column = df.columns
    if "n" in column:
        return True
    else:
        return False


def find_mean_column(df):
    column = df.columns
    if "mean" in column:
        return True
    else:
        return False


def find_stats_re_column(df):
    column = df.columns
    if "re" in column:
        return True
    else:
        return False


def remove_extra_char(df):
    column = list(df.columns)
    for col in column:
        values = []
        if "nominal" in col:
            val = df[col].to_list()
            x = val[0]
            for v in val:
                if str(v).strip() != "":
                    try:
                        v = decimal.Decimal(v.replace(",", ""))
                    except:
                        v = parse_decimal(str(v))
                    values.append(v)
                    x = v
                else:
                    values.append(x)
            df[col] = values

        elif "dilution" in col:
            val = df[col].to_list()
            x = val[0]
            for v in val:
                if str(v).strip() != "":
                    try:
                        v = decimal.Decimal(v.replace(",", ""))
                    except:
                        v = decimal.Decimal(v)
                    values.append(v)
                    x = v
                else:
                    values.append(x)
            df[col] = values

        elif "date" in col or "_id" in col:
            val = df[col].to_list()
            x = val[0]
            for v in val:
                if str(v).strip() != "":
                    values.append(v)
                    x = v
                else:
                    values.append(x)
            df[col] = values
        elif "samples" in col:
            continue
        else:
            try:
                val = df[col].to_list()
                val = [parse_decimal(str(x)) for x in val]
                df[col] = val
            except Exception as e:
                print("Error", e)

    return df


def calculate_re(values, nominal_values):
    res = []
    for index, val in enumerate(values):
        try:
            res.append(((val - nominal_values[index]) / nominal_values[index]) * 100)
        except (decimal.InvalidOperation, ZeroDivisionError, TypeError):
            res.append(Decimal(0))
    return res


def extract_dilution_factor(table_header):
    factor = []
    concentrations = list(find_nominal_conc(table_header).values())
    all_values = []
    cell_len = 0
    for text in table_header:
        parsed = re.findall(r"[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", text)
        if parsed:
            parsed = [Decimal(x) for x in parsed]
            if cell_len < len(parsed):
                cell_len = len(parsed)
            all_values.append(parsed)
    for i, values in enumerate(all_values):
        if i == 0:
            try:
                val_index = values.index(concentrations[i])
                val = values[val_index + 1]
                flag = True
            except (ValueError, IndexError):
                flag = False
        if flag:
            if cell_len == len(values):
                factor.append(values[val_index + 1])
            else:
                factor.append(Decimal(1))
        elif not flag:
            if cell_len == len(values):
                factor.append(values[val_index - 1])
            else:
                factor.append(Decimal(1))

    return factor


def split_analyte_table_df(parsed_tables: list, user_analytes: list) -> list:
    """Separate table dataframe for multiple analyte, if single table for multiple analyte"""
    analyte_table_list = []
    for annotated_table in parsed_tables:
        table_df = annotated_table['table_rows']
        table_df = remove_header(table_df)
        table_df = remove_footer(table_df)
        analyte_col = list(filter(lambda x: "analyte" in str(x).lower(), table_df.columns))
        if not analyte_col:
            return parsed_tables
        analyte_col = analyte_col[0]
        analytes = table_df[analyte_col].unique()
        for table_analyte in analytes:
            analyte_table = table_df[table_df[analyte_col] == table_analyte].reset_index(drop=True)
            annotated_table["table_rows"] = analyte_table
            analyte = [x for x in user_analytes if x.lower() == table_analyte.lower()]
            annotated_table["analyte"] = analyte[0] if analyte else table_analyte
            table_dict = annotated_table.copy()
            analyte_table_list.append(table_dict)

    return analyte_table_list


def merge_batch_performance_tables(parsed_tables: list, analytes: list) -> dict:
    table_dict = {}
    concat_tables = pd.DataFrame()
    table_type = ""
    analysis_type = ""
    tb_title = ""
    table_title = ""
    for annotated_table in parsed_tables:
        table_title = annotated_table["table_title"]
        tb_title = annotated_table["tb_title"]
        table_type = annotated_table["table_type"]
        analysis_type = annotated_table["analysis_type"]
        analyte = get_analyte(analytes, table_title + tb_title)
        table = annotated_table["table_rows"]
        table = remove_footer(remove_header(table))
        concat_tables = pd.concat([concat_tables, table], ignore_index=True)
        table_dict[analyte] = {"table_title": table_title}
    table_dict["table_type"] = table_type
    table_dict["table_rows"] = concat_tables
    table_dict["analysis_type"] = analysis_type
    table_dict["tb_title"] = tb_title
    return table_dict


def stability_name(title_array):
    if "bench-top" in title_array:
        return "Bench-top Stability"
    elif "long-term" in title_array:
        return "Long-term Stability"
    elif "short-term" in title_array:
        return "Short-term Stability"
    elif "blood" in title_array:
        return "Blood Stability"
    else:

        try:
            st_index = title_array.index("stability")
        except ValueError:
            return None
        try:
            added_str = " ".join(title_array[(st_index - 1):(st_index + 1)])
        except IndexError:
            return None
        added_str = re.sub("[(),]", "", added_str)
        added_str = added_str.title()
        return added_str


def remove_white_space(df: pd.DataFrame):
    df_obj = df.select_dtypes(['object'])
    df[df_obj.columns] = df_obj.apply(lambda x: x.str.strip())
    return df


def check_te_column(df: pd.DataFrame, batch: str):
    try:
        is_empty = df[df["run_id"] == batch]["te"].replace("", np.nan).dropna().empty
        return False if is_empty else True
    except KeyError:
        return False


def get_reported_samples_num(annotated_table):
    table_rows = annotated_table["table_rows"]
    samples = [1]
    columns = [0, 1]

    def find_analyzed(x):
        x = str(x).lower()
        return ('analyze' in x) or ('total' in x)

    def get_max_sample(table):
        max_sample = 1
        if not table.empty:
            table = table.applymap(lambda x: parse_int(x))
            max_sample_series = table.max()
            max_value = max_sample_series.max()
            try:
                max_sample = int(max_value.item())
            except AttributeError:
                pass
        return max_sample

    for col in columns:
        try:
            found_rows = table_rows[table_rows[str(col)].apply(lambda x: find_analyzed(x))]
        except KeyError:
            found_rows = table_rows[table_rows[col].apply(lambda x: find_analyzed(x))]

        if not found_rows.empty:
            samples.append(get_max_sample(found_rows))

    return max(samples)


def parse_expiry_dates(annotated_table):
    table_rows = annotated_table["table_rows"]

    found_dates = []
    expiry_date = None

    def find_analyzed(x):
        x = str(x).lower()
        return 'expir' in x or "retest" in x

    try:
        found_columns = table_rows.loc[:, table_rows.apply(lambda x: find_analyzed(x))]
        found_columns = found_columns.values.tolist()
        if found_columns:
            for column in found_columns:
                if column:
                    for element in column:
                        try:
                            tested = parse_date(str(element))
                            if tested is not None:
                                found_dates.append(tested)
                        except:
                            continue

        found_rows = table_rows[table_rows.apply(lambda x: find_analyzed(x), axis=1)]
        found_rows = found_rows.values.tolist()
        if found_rows:
            for row in found_rows:
                if row:
                    for element in row:
                        try:
                            tested = parse_date(str(element))
                            if tested is not None:
                                found_dates.append(tested)
                        except:
                            continue

        if found_dates:
            expiry_date = min(found_dates)
    except:
        pass

    return expiry_date



