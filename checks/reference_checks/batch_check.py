import re

from checks.check_result import green, yellow, red
from checks.utils import get_unique_ids


def alternate_cross_reference_batch_numbers(parsed_table_boj, batch_summary_table, analyte: str = None):
    table_title = parsed_table_boj.table_title
    table_type = parsed_table_boj.table_type
    analyte_str = ''
    result = []

    table_rows = parsed_table_boj.data_df

    if table_rows.empty:
        return [red(f'Data not found in table {analyte_str}', None, table_title, table_type)]

    if analyte is not None:
        analyte_str = f'for analyte {analyte}'
        try:
            batch_summary_table = batch_summary_table[batch_summary_table["analyte"].str.lower() == analyte.lower()]. \
                reset_index(drop=True)
        except KeyError:
            pass

    ok_result = green(f'Batch Number Check Ok, No Findings {analyte_str}', None, table_title, table_type)

    try:
        batch_numbers = list(table_rows['run_id'].unique())
    except KeyError:
        return [red(f'No batch numbers found in table {analyte_str}', None, table_title, table_type)]

    if not isinstance(batch_numbers[0], int):
        batch_numbers = get_unique_ids(batch_numbers)

    if not batch_numbers:
        return [red(f'No batch numbers found in table {analyte_str}', None, table_title, table_type)]

    sorted_batch_numbers = sorted(batch_numbers)
    if list(sorted_batch_numbers) != list(batch_numbers):
        result.append(yellow(f'Batch numbers out of order {analyte_str}', None, table_title, table_type))

    if batch_summary_table.empty:
        result.append(red(f'No Batch Summary Table for comparison {analyte_str}', None, table_title, table_type))
        return result

    summary_batch_numbers = list(batch_summary_table['run_id'])

    if not isinstance(summary_batch_numbers[0], int):
        summary_batch_numbers = get_unique_ids(summary_batch_numbers)

    if not summary_batch_numbers:
        return [red(f'No batch numbers found in Batch Summary Table {analyte_str}', None, table_title, table_type)]

    missing = set(summary_batch_numbers).difference(set(batch_numbers))

    if missing:
        missing = sorted(list(missing))
        result.append(yellow(f"Missing batch numbers: {', '.join(map(str, missing))} {analyte_str}", None, table_title,
                             table_type))

    if not result:
        return [ok_result]
    else:
        return result


def cross_validate_batch_numbers(parsed_table_obj, calibration_table, regression_table):
    result = []
    total_missing = []

    raw_table = parsed_table_obj.data_df
    if raw_table.empty:
        return result

    table_title = parsed_table_obj.table_title
    table_type = parsed_table_obj.table_type

    try:
        raw_ids = list(raw_table["run_id"].unique())

    except:
        raw_ids = []

    if not isinstance(raw_ids[0], int):
        raw_ids = get_unique_ids(raw_ids)
    if not raw_ids:
        return [red('No batch numbers found in table', None, table_title, table_type)]

    sorted_batch_numbers = sorted(raw_ids)
    if list(sorted_batch_numbers) != list(raw_ids):
        result.append(yellow(f'Batch numbers out of order', None, table_title, table_type))

    ok_result = green('Batch Number Check Ok, No Findings', None, table_title, table_type)

    if calibration_table:
        for table in calibration_table:
            calibration_ids = list(table["run_id"].unique())

            if not isinstance(calibration_ids[0], int):
                calibration_ids = get_unique_ids(calibration_ids)

            missing = [x for x in calibration_ids if x not in raw_ids]
            total_missing.extend(missing)
    else:
        result.append(red(f"No Calibration Curve Samples Table for comparison", None, table_title, table_type))

    if regression_table:
        for table in regression_table:
            regression_ids = list(table["run_id"].unique())

            if not isinstance(regression_ids[0], int):
                regression_ids = get_unique_ids(regression_ids)

            missing = [x for x in regression_ids if x not in raw_ids]
            total_missing.extend(missing)
    else:
        result.append(red(f"No Regression Model Table for comparison", None, table_title, table_type))

    if total_missing:
        missing = sorted(list(set(total_missing)))
        result.append(yellow('Missing batch numbers: %s' % ', '.join(map(str, missing)), None, table_title, table_type))

    if not result:
        return [ok_result]
    else:
        return result
