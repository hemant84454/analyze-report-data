from checks.check_result import green, yellow, red
from checks.utils import parse_date


def alternate_cross_reference_assay_date(parsed_table_obj, batch_summary_table, analyte: str = None):
    result = []
    analyte_str = ''
    table_title = parsed_table_obj.table_title
    table_type = parsed_table_obj.table_type
    table_rows = parsed_table_obj.data_df

    if table_rows.empty:
        return result

    if analyte is not None:
        analyte_str = f'for analyte {analyte}'
        try:
            batch_summary_table = batch_summary_table[batch_summary_table['analyte'].str.lower() == analyte.lower()].reset_index(drop=True)
        except KeyError:
            pass

    ok_result = green(f'Assay Date Check Ok, No Findings {analyte_str}', None, table_title, table_type)

    batch_numbers = list(table_rows['run_id'].unique())

    try:
        table_assay_dates = table_rows["assay_date"].to_list()
        table_assay_dates = [x for x in table_assay_dates if str(x).strip() != ""]
        if table_rows['assay_date'].empty or table_rows['assay_date'].dropna().empty or not table_assay_dates:
            result.append(yellow(f'No Assay Dates Found {analyte_str}', None, table_title, table_type))
            return result
    except KeyError:
        result.append(yellow(f'No Assay Dates Found {analyte_str}', None, table_title, table_type))
        return result

    if batch_summary_table.empty:
        return result

    summary_batch_numbers = list(batch_summary_table['run_id'])
    if not summary_batch_numbers:
        return result

    try:
        summary_assay_dates = batch_summary_table["assay_date"].to_list()
        summary_assay_dates = [x for x in summary_assay_dates if str(x).strip() != ""]
        if batch_summary_table['assay_date'].empty or batch_summary_table['assay_date'].dropna().empty or not summary_assay_dates:
            result.append(yellow(f'No Assay Dates found in Batch Summary Table {analyte_str}', None, table_title, table_type))
            return result
    except KeyError:
        result.append(yellow(f'No Assay Dates found in Batch Summary Table {analyte_str}', None, table_title, table_type))
        return result

    for batch in batch_numbers:
        try:
            assay_date = table_rows[table_rows['run_id'] == batch]['assay_date'].tolist()[0]
        except IndexError:
            assay_date = None

        try:
            summary_assay_date = batch_summary_table[batch_summary_table['run_id'] == batch]['assay_date'].tolist()[0]
        except IndexError:
            summary_assay_date = None

        if assay_date is not None and summary_assay_date is not None:
            try:
                if parse_date(assay_date) != parse_date(summary_assay_date):
                    result.append(
                        red(f'[Batch {batch}] Incorrect Assay date, when cross checked with Batch Run Summary Table {analyte_str}',
                            None, table_title, table_type))
            except:
                continue

    if not result:
        return [ok_result]
    else:
        return result
