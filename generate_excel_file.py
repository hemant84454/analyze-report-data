import pandas as pd
import os
import boto3
s3_cc = boto3.client('s3')


def create_excel_file(all_source_tables, s3_path, bucket_name):
    output_file = '/tmp/multi_sheet.xlsx'
    writer = pd.ExcelWriter(output_file, engine='xlsxwriter')
    c = 2
    column_a_len = 0
    column_b_len = 0
    workbook = writer.book
    notation_worksheet = workbook.add_worksheet('Notations')
    bold = workbook.add_format({'bold': True})
    blue = workbook.add_format({'color': 'blue', 'italic': True})

    for path in all_source_tables:
        notation_worksheet.write(f'A1', 'Notation', bold)
        notation_worksheet.write(f'B1', 'Meaning', bold)
        notation_worksheet.write(f'B{c}', path["notation"])
        sheet = (os.path.basename(path["path"])).split(".csv")[0]
        notation_worksheet.write_url(f'A{c}', f"internal:'{sheet}'!A1")
        notation_worksheet.write(f'A{c}', sheet, blue)
        c += 1
        if column_a_len < len(sheet):
            column_a_len = len(sheet)

        if column_b_len < len(path["notation"]):
            column_b_len = len(path["notation"])

    for table in all_source_tables:
        df = pd.read_csv(table["path"], dtype=object, keep_default_na=False)
        df = df.fillna('')
        sheet_name = (os.path.basename(table["path"])).split(".csv")[0][:31]
        df.to_excel(writer, sheet_name=sheet_name, index=False)
        worksheet = writer.sheets[sheet_name]
        (max_row, max_col) = df.shape
        column_settings = [{'header': column} for column in df.columns]
        worksheet.add_table(0, 0, max_row, max_col - 1, {'columns': column_settings})
        worksheet.set_column(0, max_col - 1, 12)

    column_settings = [{'header': 'Notation'}, {'header': 'Meaning'}]
    notation_worksheet.add_table(0, 0, c-2, 1, {'columns': column_settings})

    notation_worksheet.set_column('A:A', column_a_len)
    notation_worksheet.set_column('B:B', column_b_len)

    writer.save()
    s3_cc.upload_file(output_file, bucket_name, s3_path)