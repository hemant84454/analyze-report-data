import json
import logging
import tempfile
from itertools import groupby
from operator import itemgetter

import pandas as pd

from checks.all import run_all_checks
from checks.utils import format_date, parse_float
from checks.utils import remove_duplicate_results, isfloat, parse_date
from parseVerifiedDocxTables import parseVerifiedDocxTables
from parseVerifiedExcelTables import parseVerifiedExcelTables

null = None
logger = logging.getLogger()
logger.setLevel(logging.INFO)


def analyze_report(event):
    # excel for to write all processed dataframes
    row_count = 2  # first row for header
    excel_file = tempfile.TemporaryFile(mode="w+b")
    excel_writer = pd.ExcelWriter(excel_file, engine="xlsxwriter")
    workbook = excel_writer.book
    notation_worksheet = workbook.add_worksheet('Notations')
    bold = workbook.add_format({'bold': True})
    blue = workbook.add_format({'color': 'blue', 'italic': True})

    logger.info('event\n' + json.dumps(event))
    analysisType = event['AnalysisDetails']['analysisType']
    analysisSubtype = event['AnalysisDetails']['analysisSubtype']
    template_type = event['AnalysisDetails'].get('template', '')
    selected_tables = json.loads(event['SelectedTables'])
    input_filename = event['file_name']
    report_type = str(analysisType) + "_" + str(analysisSubtype)
    report_type = report_type.lower()
    bucket_name = event['bucket_name']
    analyte_name = ', '.join(event["AnalysisDetails"]["analyteNames"])
    lts_days = event["AnalysisDetails"]["lts80_dialog"]
    ref_date = event["AnalysisDetails"]["sampleDate"]

    #  Map options from dropdowns.json
    with open('/tmp/dropdowns.json') as f:
        data = f.read()
    mappings = json.loads(data)
    lts80_dialog = None
    lts20_dialog = 0
    extraction_dialog = None
    benchtop_dialog = None
    freeze_thaw_dialog = None
    first_patient_dialog = None
    model_dialog = ''
    anticoagulant = event["antiCoagulant"]
    samplePrep = ''
    species = event["species"]
    sample_matrix = event["matrix"]
    if str(event['AnalysisDetails']['lts80_dialog']).isdigit():
        lts80_dialog = int(event['AnalysisDetails']['lts80_dialog'])
    if str(event['AnalysisDetails']['lts20_dialog']).isdigit():
        lts20_dialog = int(event['AnalysisDetails']['lts20_dialog'])
    if isfloat(str(event['AnalysisDetails']['extraction_dialog'])):
        extraction_dialog = float(event['AnalysisDetails']['extraction_dialog'])
    if str(event['AnalysisDetails']['benchtop_dialog']):
        benchtop_dialog = parse_float(event['AnalysisDetails']['benchtop_dialog'])
    if str(event['AnalysisDetails']['freeze_thaw_dialog']).isdigit():
        freeze_thaw_dialog = int(event['AnalysisDetails']['freeze_thaw_dialog'])
    if str(event['AnalysisDetails']['first_patient_dialog']).isdigit():
        first_patient_dialog = int(event['AnalysisDetails']['first_patient_dialog'])
    grp1 = itemgetter('key')
    for key_head, grp_head in groupby(sorted(mappings, key=grp1), grp1):
        if key_head == 'regressionModel':
            for val in grp_head:
                for x in val['values']:
                    if x['name'] == event['AnalysisDetails']['regressionModel']:
                        model_dialog = x['desc']
                        break
        elif key_head == 'samplePrep':
            for val in grp_head:
                for x in val['values']:
                    if x['name'] == event['AnalysisDetails']['samplePrep']:
                        samplePrep = x['desc']
                        break
    lts20_dialog = lts80_dialog
    ref_date = parse_date(ref_date)
    multiple_analyte = event['AnalysisDetails'].get('multipleAnalytes', "False")
    multiple_analyte = True if multiple_analyte == "True" else False
    validator_options = {
        'multiple_analyte': multiple_analyte,
        'analyte_name': analyte_name,
        'report_type': report_type,
        'extraction_stability': extraction_dialog,
        'benchtop_stability': benchtop_dialog,
        'regression_model': model_dialog,
        'lts80': lts80_dialog,
        'lts20': lts20_dialog,
        'first_patient_date': first_patient_dialog,
        'freeze_thaws': freeze_thaw_dialog,
        'sample_matrix': sample_matrix,
        'anticoagulant': anticoagulant,
        'species': species,
        'samplePrep': samplePrep,
        "ref_date": ref_date,
        "template_type": template_type
    }

    old_validator_options = validator_options

    found_summary_values = {}
    mv_found_summary_values = []
    has_mv = 'no'
    has_sa = 'no'

    temp_file_list = []
    for file in input_filename['docx']:
        if "method" in file[2].lower():
            temp_file_list.append(file)
    for file in input_filename['docx']:
        if "sample analysis" in file[2].lower():
            temp_file_list.append(file)
    input_filename['docx'] = temp_file_list
    temp_file_list = []
    for file in input_filename['xlsx']:
        if 'method' in file[2].lower():
            temp_file_list.append(file)
    for file in input_filename['xlsx']:
        if "sample analysis" in file[2].lower():
            temp_file_list.append(file)
    input_filename['xlsx'] = temp_file_list
    temp_file_list = []
    for file in input_filename['pdf']:
        if 'method' in file[2].lower():
            temp_file_list.append(file)
    for file in input_filename['pdf']:
        if "sample analysis" in file[2].lower():
            temp_file_list.append(file)
    input_filename['pdf'] = temp_file_list

    verified_tables_list = json.loads(event['SelectedTables'])
    g_overall_result = []
    for file in input_filename['docx']:
        if found_summary_values:
            validator_options['extraction_stability'] = found_summary_values['found_short_stability']
            validator_options['benchtop_stability'] = found_summary_values['found_short_stability']
            if found_summary_values['found_sample_matrix']:
                validator_options['regression_model'] = found_summary_values['found_regression_weighting'][0]
            validator_options['lts80'] = found_summary_values['found_long_stability']
            validator_options['lts20'] = found_summary_values['found_long_stability']
            validator_options['freeze_thaws'] = found_summary_values['found_freeze_thaws']
            if found_summary_values['found_sample_matrix']:
                validator_options['sample_matrix'] = found_summary_values['found_sample_matrix'][0]
            if found_summary_values['found_anticoagulant']:
                validator_options['anticoagulant'] = found_summary_values['found_anticoagulant'][0]
            if found_summary_values['found_species']:
                validator_options['species'] = found_summary_values['found_species'][0]
            if found_summary_values['found_sample_prep']:
                validator_options['samplePrep'] = found_summary_values['found_sample_prep'][0]
            validator_options['found_cs_concentrations'] = found_summary_values['found_cs_concentrations']
            validator_options['found_qc_concentrations'] = found_summary_values['found_qc_concentrations']
            validator_options['found_run_length'] = found_summary_values['found_run_length']

        sub_validated_tables_list = []
        if "sample analysis" in file[2].lower():
            validator_options['report_type'] = str(analysisType) + "_sa"
            has_sa = 'yes'
            for v_table in verified_tables_list:
                if "sample" in v_table['analysis_type'].lower():
                    if str(v_table["table_index"]).isdigit():
                        sub_validated_tables_list.append(v_table)
        elif 'method' in file[2].lower():
            validator_options['report_type'] = str(analysisType) + "_mv"
            has_mv = 'yes'
            for v_table in verified_tables_list:
                if "method" in v_table['analysis_type'].lower():
                    if str(v_table["table_index"]).isdigit():
                        sub_validated_tables_list.append(v_table)
        parsed_tables = parseVerifiedDocxTables(file, bucket_name, sub_validated_tables_list)
        overall_result, found_summary_values, excel_writer, notation_worksheet, row_count = run_all_checks(
            parsed_tables, options=validator_options,
            writer=excel_writer, notation_worksheet=notation_worksheet,
            style=bold, color=blue, count=row_count)
        if 'method' in file[2].lower():
            mv_found_summary_values = found_summary_values
        for result in overall_result:
            result.table_type = str(file[2]) + '- ' + str(result.table_type)

        g_overall_result += overall_result
    for file in input_filename['xlsx']:
        sub_validated_tables_list = []
        if "sample analysis" in file[2].lower():
            validator_options['report_type'] = str(analysisType) + "_sa"
            has_sa = 'yes'
            for v_table in verified_tables_list:
                if "sample" in v_table['analysis_type'].lower():
                    if not str(v_table["table_index"]).isdigit():
                        sub_validated_tables_list.append(v_table)
        elif 'method' in file[2].lower():
            validator_options['report_type'] = str(analysisType) + "_mv"
            has_mv = 'yes'
            for v_table in verified_tables_list:
                if "method" in v_table['analysis_type'].lower():
                    if not str(v_table["table_index"]).isdigit():
                        sub_validated_tables_list.append(v_table)
        parsed_tables = parseVerifiedExcelTables(sub_validated_tables_list)
        overall_result, found_summary_values, excel_writer, notation_worksheet, row_count = run_all_checks(
            parsed_tables, options=validator_options,
            writer=excel_writer, notation_worksheet=notation_worksheet,
            style=bold, color=blue, count=row_count)
        if 'method' in file[2].lower():
            mv_found_summary_values = found_summary_values
        for result in overall_result:
            result.table_type = str(file[2]) + '- ' + str(result.table_type)

        g_overall_result += overall_result
    for file in input_filename['pdf']:
        sub_validated_tables_list = []
        if "sample analysis" in file[2].lower():
            validator_options['report_type'] = str(analysisType) + "_sa"
            has_sa = 'yes'
            for v_table in verified_tables_list:
                if "sample" in v_table['analysis_type'].lower():
                    sub_validated_tables_list.append(v_table)
        elif 'method' in file[2].lower():
            validator_options['report_type'] = str(analysisType) + "_mv"
            has_mv = 'yes'
            for v_table in verified_tables_list:
                if "method" in v_table['analysis_type'].lower():
                    sub_validated_tables_list.append(v_table)
        parsed_tables = parseVerifiedExcelTables(sub_validated_tables_list)
        overall_result, found_summary_values, excel_writer, notation_worksheet, row_count = run_all_checks(
            parsed_tables, options=validator_options,
            writer=excel_writer, notation_worksheet=notation_worksheet,
            style=bold, color=blue, count=row_count)
        if 'method' in file[2].lower():
            mv_found_summary_values = found_summary_values
        for result in overall_result:
            result.table_type = str(file[2]) + '- ' + str(result.table_type)

        g_overall_result += overall_result
    g_overall_result = remove_duplicate_results(g_overall_result)
    overall_result = [result.to_json() for result in g_overall_result]

    grouper = itemgetter('color')
    summary_statistics = dict()
    summary_statistics['red'] = 0
    summary_statistics['yellow'] = 0
    summary_statistics['green'] = 0
    for key, grp in groupby(sorted(overall_result, key=grouper), grouper):
        summary_statistics[key.lower()] = len([item for item in grp])

    grp1 = itemgetter('heading')
    grp2 = itemgetter('heading', 'tabletype')
    grp3 = itemgetter('heading', 'tabletype', 'color')
    results = []
    for key_head, grp_head in groupby(sorted(overall_result, key=grp1), grp1):
        results.append({
            'heading': key_head,
            'resultItem': []
        })
    for key_head, grp_head in groupby(sorted(overall_result, key=grp2), grp2):
        for result in results:
            if result['heading'] == key_head[0]:
                result['resultItem'].append({'tabletype': key_head[1], 'item': []})
    for key_head, grp_head in groupby(sorted(overall_result, key=grp3), grp3):
        for result in results:
            if result['heading'] == key_head[0]:
                for item in result['resultItem']:
                    if item['tabletype'] == key_head[1]:
                        item['item'].append({'color': key_head[2], 'message': [it['message'] for it in grp_head]})

    ref_date = format_date(ref_date)
    old_validator_options["ref_date"] = ref_date
    user_input_message_items = []
    for key in old_validator_options:
        if str(key) == 'report_type' or str(key) == 'first_patient_date' or str(
                key) == 'found_cs_concentrations' or str(key) == 'found_qc_concentrations' or \
                str(key) == 'lts20' or str(key) == 'found_run_length':
            continue
        elif str(key) == 'lts80':
            key_string = 'Long-Term Stability (Days)'
        elif str(key) == 'benchtop_stability':
            key_string = 'Benchtop Stability (Hours)'
        elif str(key) == 'extraction_stability':
            key_string = 'Extraction Stability (Hours)'
        elif str(key) == 'samplePrep':
            key_string = 'Sample Prep'
        else:
            key_string = str(key).title().replace('_', ' ')

        user_input_message_items.append(str(key_string) + ': ' + str(old_validator_options[key]).replace('None', ""))

    user_input_item = [{'color': ' ', 'message': user_input_message_items}]
    user_input_resultItem = [{'tabletype': 'User Inputs', 'item': user_input_item}]
    user_input = [{'heading': 'Inputs', 'resultItem': user_input_resultItem}]

    if mv_found_summary_values and has_mv == 'yes' and has_sa == 'yes':
        user_input_message_items = []
        for key in mv_found_summary_values:
            if str(key) == 'table_title' or str(key) == 'table_type' or "inter" in str(key) or "intra" in str(key):
                continue
            elif str(key) == 'found_long_stability':
                key_string = 'Long-Term Stability (Days)'
            elif str(key) == 'found_short_stability':
                key_string = 'Short-Term Stability (Hours)'
            elif str(key) == 'samplePrep':
                key_string = 'Sample Prep'
            else:
                key_string = str(key).capitalize().replace('_', ' ')

            value_string = str(mv_found_summary_values[key]).replace("Decimal(", "").replace("[", "").replace("]", ""). \
                replace("'", "").replace(")", "")
            if value_string != '0' and value_string:
                user_input_message_items.append(str(key_string) + ': ' + str(value_string))

        user_input_item = [{'color': ' ', 'message': user_input_message_items}]
        user_input_resultItem = [{'tabletype': 'Found in Method Validation', 'item': user_input_item}]
        val_user_input = [{'heading': 'Method Validation Values', 'resultItem': user_input_resultItem}]
        user_input = user_input + val_user_input

    results = user_input + results
    for resultitem in results:
        red_count = 0
        green_count = 0
        yellow_count = 0

        for result in resultitem["resultItem"]:

            for item in result['item']:
                if item["color"].lower().strip() == "red":
                    red_count = len(item["message"])
                if item["color"].lower().strip() == "green":
                    green_count = len(item["message"])
                if item["color"].lower().strip() == "yellow":
                    yellow_count = len(item["message"])

        resultitem["green"] = green_count
        resultitem["yellow"] = yellow_count
        resultitem["red"] = red_count

    column_settings = [{'header': 'Sheet Name'}, {'header': 'Table Type'}, {'header': 'Table Title'}]
    notation_worksheet.add_table(0, 0, row_count - 2, 1, {'columns': column_settings})
    excel_writer.save()
    excel_file.seek(0)
    return [results, summary_statistics, excel_file]
