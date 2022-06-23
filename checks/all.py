from collections import defaultdict

import pandas as pd
from checks.check_across_titles import check_across_titles
from checks.check_result import green, yellow, red
from checks.reference_checks.batch_check import alternate_cross_reference_batch_numbers, cross_validate_batch_numbers
from checks.reference_checks.check_multi_analyte import check_multi_analyte_across_title
from checks.reference_checks.date_check import alternate_cross_reference_assay_date
from checks.reference_checks.units_check import cross_check_unit
from checks.utils import get_reported_samples_num, next_table_title
from checks.utils import parse_expiry_dates, write_to_excel
from checks.utils import remove_duplicate_results, split_analyte_table_df
from table.table_registration import Stability, Calibration, ExtractionRecovery, DilutionLinearity, \
    AccuracyAndPrecision, QualityControl, Selectivity, RegressionModel, Specificity, BatchPerformance, MatrixEffect, \
    CarryOver, DilutionIntegrity, ISR, Sample, HookEffect, ReanalysisAndReassay, BloodStability, StockStability, \
    Interference, ValidationSummary, Parallelism, PPDSelectivity, PPDAccuracyAndPrecision
from table.table_registration import FlagProperties


def run_all_checks(tables, options=dict, **kwargs):
    """
    kwargs
    writer: xlsx writer
    notation_worksheet: Excel sheet
    style: style of font
    color: color of font
    count: row number
    """
    error_messages = FlagProperties("Errors")
    excel_writer = kwargs["writer"]
    notation_worksheet = kwargs["notation_worksheet"]
    style = kwargs["style"]
    color = kwargs["color"]
    row_count = kwargs["count"]

    if options is None:
        options = {}
    analytes = options.get("analyte_name").split(", ")
    multiple_analyte = options.get("multiple_analyte", False)
    overall_result = []
    blood_stability_tables = []
    stock_stability_tables = []
    interference_tables = []
    batch_performance_table = pd.DataFrame()
    batch_performance_dict = {}
    qc_tables = []
    co_tables = []
    bc_tables = []
    calibration_table = []
    model_tables = []
    regression_table = []
    reassay_tables = []
    isr_tables = []
    samples_tables = []
    stability_tables = []
    selectivity_tables = []
    ap_tables = []
    dil_linear_tables = []
    dil_integrity_tables = []
    recovery_tables = []
    matrix_effect_tables = []
    val_summary_tables = []
    specificity_tables = []
    hook_tables = []
    exp_tables = []
    parallelism_tables = []
    ap_stats = {}
    if multiple_analyte:
        pass_failed_batches = defaultdict(lambda: defaultdict(list))
        total_batches = defaultdict(lambda: defaultdict(list))

    else:
        total_batches = {"bc_batches": [], "regression_batches": [], "qc_batches": [], "ap_batches": []}

        pass_failed_batches = {"bc_pass_batches": [], "bc_failed_batches": [], "regression_pass_batches": [],
                               "regression_failed_batches": [], "qc_pass_batches": [], "qc_failed_batches": [],
                               "ap_pass_batches": [], "ap_failed_batches": []}

    if multiple_analyte:
        units = defaultdict()
        all_units = defaultdict(list)
        dil_integrity_dict = defaultdict(dict)
        LLOQ = dict()
        ULOQ = dict()
        qc_header_values = defaultdict(list)
        cs_header_values = defaultdict(list)
    else:
        units = ''
        all_units = []
        dil_integrity_dict = dict()
        LLOQ = 0
        ULOQ = 0
        qc_header_values = []
        cs_header_values = []

    all_found_summary_values = {"table_title": "Validation Summary", "table_type": "Validation Summary",
                                "found_sample_matrix": [], "found_anticoagulant": [], "found_sample_prep": [],
                                "found_species": [],
                                "found_regression_weighting": [], "found_LLOQ": 0, "found_ULOQ": 0,
                                "found_freeze_thaws": 0,
                                "found_short_stability": 0, "found_long_stability": 0, "found_run_length": 0,
                                "found_inter_precisions": [],
                                "found_intra_precisions": [], "found_inter_accuracies": [],
                                "found_intra_accuracies": [],
                                "found_cs_concentrations": [], "found_qc_concentrations": []}

    possible_sample_matrices = ['Serum', 'Urine', 'Plasma', 'Tissue']
    possible_anticoagulants = ['Sodium Heparin', 'Na Heparin', 'K1EDTA', 'K2EDTA', 'K₂EDTA', 'K3EDTA', 'K1-EDTA',
                               'K1EDTA', 'K2-EDTA', 'K3-EDTA', 'K₃EDTA', 'Heparin', "Sodium Citrate"]
    possible_sample_species = ['Human', 'Monkey', 'Rabbit', 'Rat', 'Dog']
    possible_sample_preps = ['Liquid-liquid extraction', 'Solid phase extraction', 'Protein precipitation']
    possible_regressions = ['Linear', 'Quadratic', 'Quad', '4PL', '5PL']

    analysts = []
    num_analyte = 1
    total_samples = 1

    if multiple_analyte:
        multi_analyte_total_samples = defaultdict()
    else:
        multi_analyte_total_samples = None

    ft_cycle = options.get("freeze_thaws")
    report_type = options["report_type"]
    template_type = options.get('template_type', '').lower()
    if report_type is None:
        report_type = 'sm_sa'
    else:
        report_type = str(report_type).lower()

    if multiple_analyte:
        sample_table = defaultdict(pd.DataFrame)
    else:
        sample_table = pd.DataFrame()

    for table in tables:
        table_type = table['table_type']
        table["analysis_type"] = report_type

        if table_type == "Overall and Individual Batch Performance":
            batch_performance_dict = table

        elif table_type == "Expiration Information":
            exp_tables.append(table)
        elif table_type == "Calibration Curve Samples":
            bc_tables.append(table)
        elif table_type == "Quality Control Samples":
            qc_tables.append(table)
        elif table_type == "Regression Model":
            model_tables.append(table)

        if 'sa' in report_type:
            if table_type == "Reanalysis and Reassay":
                reassay_tables.append(table)
            elif table_type == "Incurred Sample Reproducibility (ISR)":
                isr_tables.append(table)
            elif table_type == "Samples Table":
                samples_tables.append(table)
            elif table_type == "Samples Analyzed":
                found_total_samples = get_reported_samples_num(table)
                if found_total_samples > total_samples:
                    total_samples = found_total_samples

        if 'mv' in report_type:
            if table_type == "Stability":
                stability_tables.append(table)
            elif table_type == "Blood Stability":
                blood_stability_tables.append(table)
            elif table_type == "Stock Solution Stability":
                stock_stability_tables.append(table)
            elif table_type == "Selectivity":
                selectivity_tables.append(table)
            elif table_type == "Accuracy and Precision":
                ap_tables.append(table)
            elif table_type == "Dilution Integrity":
                dil_integrity_tables.append(table)
            elif table_type == "Dilution Linearity":
                dil_linear_tables.append(table)
            elif table_type == "Extraction Recovery":
                recovery_tables.append(table)
            elif table_type == "Matrix Effect":
                matrix_effect_tables.append(table)
            elif table_type == "Summary":
                val_summary_tables.append(table)
            elif table_type == "Specificity":
                specificity_tables.append(table)
            elif table_type == "Carry-Over":
                co_tables.append(table)
            elif table_type == "Interference Test":
                interference_tables.append(table)
            elif table_type == "Hook Effect":
                hook_tables.append(table)
            elif table_type == "Parallelism":
                parallelism_tables.append(table)

    expiration_dates = []
    expiration_date = None
    if exp_tables:
        for table in exp_tables:
            expiration_dates.append(parse_expiry_dates(table))
        if expiration_dates:
            expiration_dates = list(filter(None, expiration_dates))
            if expiration_dates:
                expiration_date = min(expiration_dates)
            else:
                expiration_date = None

    if batch_performance_dict:
        batch_performance_obj = BatchPerformance(batch_performance_dict, options, expiration_date, analytes,
                                                 multiple_analyte, template_type)
        batch_performance_obj.validate_table_format()
        if batch_performance_obj.valid_format:
            batch_performance_obj.process_table()
            excel_writer, notation_worksheet, row_count = write_to_excel(excel_writer, notation_worksheet,
                                                                         batch_performance_obj, color=color,
                                                                         style=style, count=row_count)
            if batch_performance_obj.valid_data:
                analysts = batch_performance_obj.analysts
                num_analyte = batch_performance_obj.num_analyte
                batch_performance_table = batch_performance_obj.data_df
                overall_result += batch_performance_obj.check_ref_date_and_instrument()
                overall_result += batch_performance_obj.batch_summary_check_batch_numbers()
                overall_result += batch_performance_obj.batch_summary_check_extraction_stability()
                overall_result += batch_performance_obj.batch_summary_check_expiry()
                overall_result += batch_performance_obj.stability_date_check()

        overall_result += batch_performance_obj.result
    else:
        overall_result += [red('Could not find table', None, 'Overall and Individual Batch Performance',
                               'Overall and Individual Batch Performance')]

    if bc_tables:
        for idx, table in enumerate(bc_tables):
            calibration_obj = Calibration(table, analytes, template_type)
            calibration_obj.validate_table_format()
            if calibration_obj.valid_format:
                try:
                    calibration_obj.process_table()
                    excel_writer, notation_worksheet, row_count = write_to_excel(excel_writer, notation_worksheet,
                                                                                 calibration_obj, color=color,
                                                                                 style=style, count=row_count)
                    if calibration_obj.valid_data:
                        calibration_obj.validate_table_data()
                        calibration_table.append(calibration_obj.data_df)
                        analyte = calibration_obj.analyte
                        if multiple_analyte:
                            units[analyte] = calibration_obj.conc_unit
                            cs_header_values[analyte] += calibration_obj.nominal_val
                            overall_result += alternate_cross_reference_batch_numbers(calibration_obj,
                                                                                      batch_performance_table, analyte)
                            overall_result += alternate_cross_reference_assay_date(calibration_obj,
                                                                                   batch_performance_table, analyte)
                            ULOQ[analyte] = calibration_obj.ULOQ
                            LLOQ[analyte] = calibration_obj.LLOQ
                        else:
                            LLOQ = calibration_obj.LLOQ
                            ULOQ = calibration_obj.ULOQ
                            units = calibration_obj.conc_unit
                            cs_header_values += calibration_obj.nominal_val
                            overall_result += alternate_cross_reference_batch_numbers(calibration_obj,
                                                                                      batch_performance_table)
                            overall_result += alternate_cross_reference_assay_date(calibration_obj,
                                                                                   batch_performance_table)

                        if multiple_analyte and analyte is not None:
                            pass_failed_batches["bc_pass_batches"][analyte] += calibration_obj.pass_batches
                            pass_failed_batches["bc_failed_batches"][analyte] += calibration_obj.failed_batches
                            total_batches["bc_batches"][analyte] += calibration_obj.total_batches
                        elif num_analyte == 1:
                            pass_failed_batches["bc_pass_batches"] += calibration_obj.pass_batches
                            pass_failed_batches["bc_failed_batches"] += calibration_obj.failed_batches
                            total_batches["bc_batches"] += calibration_obj.total_batches
                except:
                    message = error_messages.get_message("data_error")
                    overall_result.append(red(message, None, calibration_obj.table_title, calibration_obj.table_type))

            overall_result += calibration_obj.result
    else:
        overall_result += [red('Could not find table', None, 'Calibration Curve Samples', 'Calibration Curve Samples')]

    reported_sample = 0
    if 'sa' in report_type:

        if samples_tables:

            if multiple_analyte and len(samples_tables) == 1:
                samples_tables = split_analyte_table_df(samples_tables, analytes)

            for table in samples_tables:
                sample_obj = Sample(table, total_samples, analytes, multiple_analyte, LLOQ, units, template_type)
                sample_obj.validate_table_format()
                if sample_obj.valid_format:
                    try:
                        sample_obj.process_table()
                        excel_writer, notation_worksheet, row_count = write_to_excel(excel_writer, notation_worksheet,
                                                                                     sample_obj, color=color,
                                                                                     style=style, count=row_count)
                        if sample_obj.valid_data:
                            sample_obj.check_samples()
                            if multiple_analyte:
                                sample_table[sample_obj.analyte] = sample_obj.data_df
                                all_units[sample_obj.analyte] += sample_obj.get_concentration_units()
                                multi_analyte_total_samples[sample_obj.analyte] = sample_obj.samples
                                analyte_str = f"for analyte {sample_obj.analyte}"
                            else:
                                analyte_str = ""
                                sample_table = sample_obj.data_df
                                all_units.extend(sample_obj.get_concentration_units())

                            if total_samples > sample_obj.samples != 0:
                                tested_samples = sample_obj.samples
                                overall_result.append(
                                    red(f"Received {total_samples} (in sample listing table) but only tested {tested_samples} (in sample result table) {analyte_str}",
                                        None, "Information", "Entire Report"))

                                total_samples = tested_samples
                            elif total_samples <= 1:
                                total_samples = sample_obj.samples

                            sample_obj.check_bloq()
                    except:
                        message = error_messages.get_message("data_error")
                        overall_result.append(red(message, None, sample_obj.table_title, sample_obj.table_type))

                overall_result += sample_obj.result
        else:
            overall_result += [red('Could not find table', None, 'Samples Table', 'Samples Table')]

    if cs_header_values:
        if multiple_analyte:
            cs_header_values = {key: sorted(list(set(val))) for key, val in cs_header_values.items()}
        else:
            cs_header_values = sorted(list(set(cs_header_values)))

    all_found_summary_values['found_cs_concentrations'] = cs_header_values

    if 'sa' in report_type:
        if isr_tables:
            for table in isr_tables:
                isr_obj = ISR(table, sample_table, total_samples, LLOQ, units, analytes, multiple_analyte,
                              multi_analyte_total_samples, template_type)
                isr_obj.validate_table_format()
                if isr_obj.valid_format:
                    try:
                        isr_obj.process_table()
                        excel_writer, notation_worksheet, row_count = write_to_excel(excel_writer, notation_worksheet,
                                                                                     isr_obj, color=color, style=style,
                                                                                     count=row_count)
                        if multiple_analyte:
                            all_units[isr_obj.analyte] += [isr_obj.conc_unit]
                        else:
                            all_units.append(isr_obj.conc_unit)
                        if isr_obj.valid_data:
                            isr_obj.validate_table_data()
                    except:
                        message = error_messages.get_message("data_error")
                        overall_result.append(red(message, None, isr_obj.table_title, isr_obj.table_type))

                overall_result += isr_obj.result

        else:
            overall_result += [red('Could not find table', None, 'Incurred Sample Reproducibility (ISR)',
                                   'Incurred Sample Reproducibility (ISR)')]

    if model_tables:
        if multiple_analyte and len(model_tables) == 1:
            model_tables = split_analyte_table_df(model_tables, analytes)

        for table in model_tables:
            regression_obj = RegressionModel(table, analytes, multiple_analyte, template_type)
            regression_obj.validate_table_format()
            if regression_obj.valid_format:
                regression_obj.process_table()
                excel_writer, notation_worksheet, row_count = write_to_excel(excel_writer, notation_worksheet,
                                                                             regression_obj, color=color, style=style,
                                                                             count=row_count)
                if regression_obj.valid_data:
                    try:
                        regression_obj.validate_table_data(options)
                        regression_table.append(regression_obj.data_df)
                        analyte = regression_obj.analyte

                        if multiple_analyte and analyte is not None:
                            pass_failed_batches["regression_failed_batches"][analyte] += regression_obj.failed_batches
                            pass_failed_batches["regression_pass_batches"][analyte] += regression_obj.pass_batches
                            total_batches["regression_batches"][analyte] += regression_obj.total_batches
                            overall_result += alternate_cross_reference_batch_numbers(regression_obj,
                                                                                      batch_performance_table, analyte)
                            overall_result += alternate_cross_reference_assay_date(regression_obj,
                                                                                   batch_performance_table, analyte)
                        elif num_analyte == 1:
                            pass_failed_batches["regression_failed_batches"] += regression_obj.failed_batches
                            pass_failed_batches["regression_pass_batches"] += regression_obj.pass_batches
                            total_batches["regression_batches"] += regression_obj.total_batches
                            overall_result += alternate_cross_reference_batch_numbers(regression_obj,
                                                                                      batch_performance_table)
                            overall_result += alternate_cross_reference_assay_date(regression_obj,
                                                                                   batch_performance_table)
                    except:
                        message = error_messages.get_message("data_error")
                        overall_result.append(red(message, None, regression_obj.table_title, regression_obj.table_type))

            overall_result += regression_obj.result
    else:
        overall_result += [red('Could not find table', None, 'Regression Model', 'Regression Model')]

    if 'mv' in report_type:
        if interference_tables:
            for table in interference_tables:
                interference_obj = Interference(table, analytes, template_type)
                interference_obj.validate_table_format()
                if interference_obj.valid_format:
                    try:
                        interference_obj.process_table()
                        excel_writer, notation_worksheet, row_count = write_to_excel(excel_writer, notation_worksheet,
                                                                                     interference_obj, color=color,
                                                                                     style=style, count=row_count)
                        if multiple_analyte:
                            all_units[interference_obj.analyte] += [interference_obj.conc_unit]
                        else:
                            all_units.append(interference_obj.conc_unit)
                        if interference_obj.valid_data:
                            interference_obj.validate_table_data()
                    except:
                        message = error_messages.get_message("data_error")
                        overall_result.append(
                            red(message, None, interference_obj.table_title, interference_obj.table_type))

                overall_result += interference_obj.result
        else:
            overall_result += [red('Could not find table', None, 'Interference Test', 'Interference Test')]

        if ap_tables:

            def validate_accuracy_and_precision_table(result, stats, header_values):
                ap_obj.data_df = pd.concat([ap_obj.data_df, data_table]).reset_index(drop=True)
                ap_obj.static_df = pd.concat([ap_obj.static_df, intra_stats]).reset_index(drop=True)
                ap_obj.final_static_df = pd.concat([ap_obj.final_static_df, inter_stats]).reset_index(drop=True)
                if not ap_obj.found_duplicate_run and ap_obj.valid_data:
                    ap_obj.validate_table_data()
                    ap_header_values = ap_obj.header_values
                    analyte = ap_obj.analyte
                    result += cross_validate_batch_numbers(ap_obj, calibration_table, regression_table)

                    if multiple_analyte and analyte is not None:
                        pass_failed_batches["ap_pass_batches"][analyte] += ap_obj.passed_batch_numbers
                        pass_failed_batches["ap_failed_batches"][analyte] += ap_obj.failed_batch_numbers
                        total_batches["ap_batches"][analyte] += ap_obj.total_batches
                        header_values[analyte] += list(ap_header_values)

                    elif num_analyte == 1:
                        pass_failed_batches["ap_pass_batches"] += ap_obj.passed_batch_numbers
                        pass_failed_batches["ap_failed_batches"] += ap_obj.failed_batch_numbers
                        total_batches["ap_batches"] += ap_obj.total_batches
                        header_values += list(ap_header_values)

                    if "reject" not in ap_obj.table_title.lower() and "fail" not in ap_obj.table_title.lower():
                        if multiple_analyte and analyte is not None:
                            stats[analyte] = {'ap_header_values': ap_header_values,
                                              'intra_precisions': ap_obj.intra_precisions,
                                              'intra_accuracies': ap_obj.intra_accuracies,
                                              'inter_precisions': ap_obj.inter_precisions,
                                              'inter_accuracies': ap_obj.inter_accuracies,
                                              'inter_tes': ap_obj.inter_tes}
                        else:
                            stats = {'ap_header_values': ap_header_values,
                                     'intra_precisions': ap_obj.intra_precisions,
                                     'intra_accuracies': ap_obj.intra_accuracies,
                                     'inter_precisions': ap_obj.inter_precisions,
                                     'inter_accuracies': ap_obj.inter_accuracies, 'inter_tes': ap_obj.inter_tes}
                if multiple_analyte:
                    all_units[ap_obj.analyte] += [ap_obj.conc_unit]
                else:
                    all_units.append(ap_obj.conc_unit)
                return stats, result, qc_header_values

            data_table = pd.DataFrame()
            intra_stats = pd.DataFrame()
            inter_stats = pd.DataFrame()
            for index, table in enumerate(ap_tables):
                next_title = next_table_title(ap_tables, index + 1)
                if template_type == "ppd":
                    ap_obj = PPDAccuracyAndPrecision(table, analytes, template_type, LLOQ, ULOQ)
                else:
                    ap_obj = AccuracyAndPrecision(table, analytes, template_type, LLOQ, ULOQ)
                ap_obj.validate_table_format()
                try:
                    if ap_obj.valid_format:
                        ap_obj.process_table()
                        excel_writer, notation_worksheet, row_count = write_to_excel(excel_writer, notation_worksheet,
                                                                                     ap_obj, color=color, style=style,
                                                                                     count=row_count)
                        if next_title != ap_obj.table_title:
                            ap_stats, overall_result, qc_header_values = validate_accuracy_and_precision_table(
                                overall_result, ap_stats,
                                qc_header_values)
                            data_table = pd.DataFrame()
                            intra_stats = pd.DataFrame()
                            inter_stats = pd.DataFrame()

                        elif len(ap_tables) > 1:
                            data_table = pd.concat([data_table, ap_obj.data_df]).reset_index(drop=True)
                            intra_stats = pd.concat([intra_stats, ap_obj.static_df]).reset_index(drop=True)
                            inter_stats = pd.concat([inter_stats, ap_obj.final_static_df]).reset_index(drop=True)
                    else:
                        if next_title != ap_obj.table_title and not data_table.empty:
                            ap_stats, overall_result, qc_header_values = validate_accuracy_and_precision_table(
                                overall_result, ap_stats,
                                qc_header_values)
                            data_table = pd.DataFrame()
                            intra_stats = pd.DataFrame()
                            inter_stats = pd.DataFrame()
                except Exception as e:
                    message = error_messages.get_message("data_error")
                    overall_result.append(red(message, None, ap_obj.table_title, ap_obj.table_type))
                overall_result += ap_obj.result

        else:
            overall_result += [red('Could not find table', None, 'Accuracy and Precision', 'Accuracy and Precision')]

        if hook_tables:
            for table in hook_tables:
                hook_obj = HookEffect(table, multiple_analyte, analytes, template_type)
                hook_obj.validate_table_format()
                if hook_obj.valid_format:
                    try:
                        hook_obj.process_table()
                        excel_writer, notation_worksheet, row_count = write_to_excel(excel_writer, notation_worksheet,
                                                                                     hook_obj, color=color, style=style,
                                                                                     count=row_count)
                        if multiple_analyte:
                            all_units[hook_obj.analyte] += [hook_obj.conc_unit]
                        else:
                            all_units.append(hook_obj.conc_unit)
                        if hook_obj.valid_data:
                            hook_obj.validate_table_data(LLOQ, ULOQ)
                    except:
                        message = error_messages.get_message("data_error")
                        overall_result.append(red(message, None, hook_obj.table_title, hook_obj.table_type))

                overall_result += hook_obj.result
        elif "lm" in report_type:
            overall_result += [red('Could not find table', None, 'Hook Effect', 'Hook Effect')]

        if co_tables:
            for table in co_tables:
                co_obj = CarryOver(table, template_type)
                co_obj.validate_table_format()
                if co_obj.valid_format:
                    try:
                        co_obj.process_table()
                        excel_writer, notation_worksheet, row_count = write_to_excel(excel_writer, notation_worksheet,
                                                                                     co_obj, color=color, style=style,
                                                                                     count=row_count)
                        if co_obj.valid_data:
                            co_obj.validate_table_data()
                    except:
                        message = error_messages.get_message("data_error")
                        overall_result.append(red(message, None, co_obj.table_title, co_obj.table_type))

                overall_result += co_obj.result
        else:
            overall_result += [red('Could not find table', None, 'Carry-Over', 'Carry-Over')]

        if parallelism_tables:
            for table in parallelism_tables:
                parallelism_obj = Parallelism(table, analytes, template_type)
                parallelism_obj.validate_table_format()
                if parallelism_obj.valid_format:
                    try:
                        parallelism_obj.process_table()
                        excel_writer, notation_worksheet, row_count = write_to_excel(excel_writer, notation_worksheet,
                                                                                     parallelism_obj, color=color,
                                                                                     style=style, count=row_count)
                        if multiple_analyte:
                            all_units[parallelism_obj.analyte] += [parallelism_obj.conc_unit]
                        else:
                            all_units.append(parallelism_obj.conc_unit)
                        if parallelism_obj.valid_data:
                            parallelism_obj.validate_table_data()
                            overall_result += alternate_cross_reference_assay_date(parallelism_obj,
                                                                                   batch_performance_table)
                    except:
                        message = error_messages.get_message("data_error")
                        overall_result.append(
                            red(message, None, parallelism_obj.table_title, parallelism_obj.table_type))

                overall_result += parallelism_obj.result

    if qc_tables:
        for table in qc_tables:
            qc_obj = QualityControl(table, analytes, template_type)
            qc_obj.validate_table_format()
            if qc_obj.valid_format:
                try:
                    qc_obj.process_table()
                    excel_writer, notation_worksheet, row_count = write_to_excel(excel_writer, notation_worksheet,
                                                                                 qc_obj, color=color, style=style,
                                                                                 count=row_count)
                    if multiple_analyte:
                        all_units[qc_obj.analyte] += [qc_obj.conc_unit]
                    else:
                        all_units.append(qc_obj.conc_unit)
                    if qc_obj.valid_data:
                        qc_obj.validate_table_data()
                        analyte = qc_obj.analyte

                        if multiple_analyte and analyte is not None:
                            pass_failed_batches["qc_pass_batches"][analyte] += qc_obj.passed_batches
                            pass_failed_batches["qc_failed_batches"][analyte] += qc_obj.failed_batches
                            total_batches["qc_batches"][analyte] += qc_obj.total_batches
                            overall_result += alternate_cross_reference_batch_numbers(qc_obj, batch_performance_table,
                                                                                      analyte)
                            overall_result += alternate_cross_reference_assay_date(qc_obj, batch_performance_table,
                                                                                   analyte)
                        elif not multiple_analyte:
                            pass_failed_batches["qc_pass_batches"] += qc_obj.passed_batches
                            pass_failed_batches["qc_failed_batches"] += qc_obj.failed_batches
                            total_batches["qc_batches"] += qc_obj.total_batches
                            overall_result += alternate_cross_reference_batch_numbers(qc_obj, batch_performance_table)
                            overall_result += alternate_cross_reference_assay_date(qc_obj, batch_performance_table)
                except:
                    message = error_messages.get_message("data_error")
                    overall_result.append(red(message, None, qc_obj.table_title, qc_obj.table_type))

            overall_result += qc_obj.result
    else:
        overall_result += [red('Could not find table', None, 'Quality Control Samples', 'Quality Control Samples')]

    if qc_header_values:
        if multiple_analyte:
            qc_header_values = {key: sorted(list(set(val))) for key, val in qc_header_values.items()}
        else:
            qc_header_values = sorted(list(set(qc_header_values)))

    if 'sa' in report_type:
        if reassay_tables:
            for table in reassay_tables:
                reassay_obj = ReanalysisAndReassay(table, total_samples, sample_table, options, analytes,
                                                   multiple_analyte, multi_analyte_total_samples, template_type)
                reassay_obj.validate_table_format()
                if reassay_obj.valid_format:
                    try:
                        reassay_obj.process_table()
                        excel_writer, notation_worksheet, row_count = write_to_excel(excel_writer, notation_worksheet,
                                                                                     reassay_obj, color=color,
                                                                                     style=style, count=row_count)
                        if multiple_analyte:
                            all_units[reassay_obj.analyte] += [reassay_obj.conc_unit]
                        else:
                            all_units.append(reassay_obj.conc_unit)
                        if reassay_obj.valid_data:
                            reassay_obj.validate_table_data()
                    except:
                        message = error_messages.get_message("data_error")
                        overall_result.append(red(message, None, reassay_obj.table_title, reassay_obj.table_type))

                overall_result += reassay_obj.result
        else:
            overall_result += [red('Could not find table', None, 'Reanalysis and Reassay', 'Reanalysis and Reassay')]

    if 'mv' in report_type:
        if stability_tables:
            for table in stability_tables:
                stability_obj = Stability(table, analytes, template_type)
                stability_obj.validate_table_format()
                if stability_obj.valid_format:
                    stability_obj.process_table()
                    try:
                        excel_writer, notation_worksheet, row_count = write_to_excel(excel_writer, notation_worksheet,
                                                                                     stability_obj, color=color,
                                                                                     style=style, count=row_count)
                        if multiple_analyte:
                            all_units[stability_obj.analyte] += [stability_obj.conc_unit]
                        else:
                            all_units.append(stability_obj.conc_unit)
                        if stability_obj.valid_data:
                            stability_obj.validate_table_data()
                            stability_obj.validate_freeze_thaw_number(ft_cycle)
                    except:
                        message = error_messages.get_message("data_error")
                        overall_result.append(red(message, None, stability_obj.table_title, stability_obj.table_type))

                overall_result += stability_obj.result

        else:
            overall_result += [red('Could not find table', None, 'Stability', 'Stability')]

        if blood_stability_tables:
            for table in blood_stability_tables:
                blood_stability_obj = BloodStability(table, analytes, template_type)
                blood_stability_obj.validate_table_format()
                if blood_stability_obj.valid_format:
                    try:
                        blood_stability_obj.process_table()
                        excel_writer, notation_worksheet, row_count = write_to_excel(excel_writer, notation_worksheet,
                                                                                     blood_stability_obj, color=color,
                                                                                     style=style, count=row_count)
                        if multiple_analyte:
                            all_units[blood_stability_obj.analyte] += [blood_stability_obj.conc_unit]
                        else:
                            all_units.append(blood_stability_obj.conc_unit)
                        if blood_stability_obj.valid_data:
                            blood_stability_obj.validate_table_data()
                            blood_stability_obj.validate_table_title()
                    except:
                        message = error_messages.get_message("data_error")
                        overall_result.append(
                            red(message, None, blood_stability_obj.table_title, blood_stability_obj.table_type))

                overall_result += blood_stability_obj.result
        else:
            overall_result += [red('Could not find table', None, 'Blood Stability', 'Blood Stability')]

        if stock_stability_tables:
            for table in stock_stability_tables:
                stock_stability_obj = StockStability(table, analytes, template_type)
                stock_stability_obj.validate_table_format()
                if stock_stability_obj.valid_format:
                    try:
                        stock_stability_obj.process_table()
                        excel_writer, notation_worksheet, row_count = write_to_excel(excel_writer, notation_worksheet,
                                                                                     stock_stability_obj, color=color,
                                                                                     style=style, count=row_count)
                        if multiple_analyte:
                            all_units[stock_stability_obj.analyte] += [stock_stability_obj.conc_unit]
                        else:
                            all_units.append(stock_stability_obj.conc_unit)
                        if stock_stability_obj.valid_data:
                            stock_stability_obj.validate_table_data()
                    except:
                        message = error_messages.get_message("data_error")
                        overall_result.append(
                            red(message, None, stock_stability_obj.table_title, stock_stability_obj.table_type))

                overall_result += stock_stability_obj.result
        else:
            overall_result += [red('Could not find table', None, 'Stock Stability', 'Stock Stability')]

        if specificity_tables:
            for table in specificity_tables:
                specificity_obj = Specificity(table, multiple_analyte, analytes, template_type)
                specificity_obj.validate_table_format()
                if specificity_obj.valid_format:
                    try:
                        specificity_obj.process_table()
                        excel_writer, notation_worksheet, row_count = write_to_excel(excel_writer, notation_worksheet,
                                                                                     specificity_obj, color=color,
                                                                                     style=style, count=row_count)
                        if multiple_analyte:
                            all_units[specificity_obj.analyte] += [specificity_obj.conc_unit]
                        else:
                            all_units.append(specificity_obj.conc_unit)
                        if specificity_obj.valid_data:
                            specificity_obj.validate_table_data(LLOQ, ULOQ)
                    except:
                        message = error_messages.get_message("data_error")
                        overall_result.append(
                            red(message, None, specificity_obj.table_title, specificity_obj.table_type))

                overall_result += specificity_obj.result
        else:
            overall_result += [red('Could not find table', None, 'Specificity', 'Specificity')]

        if selectivity_tables:
            all_selectivity_info = list()
            for index, table in enumerate(selectivity_tables):
                if template_type == "ppd":
                    selectivity_obj = PPDSelectivity(table, LLOQ, analytes, multiple_analyte, template_type)
                else:
                    selectivity_obj = Selectivity(table, LLOQ, analytes, multiple_analyte, template_type)

                selectivity_obj.validate_table_format()
                if selectivity_obj.valid_format:
                    try:
                        selectivity_obj.process_table()
                        excel_writer, notation_worksheet, row_count = write_to_excel(excel_writer, notation_worksheet,
                                                                                     selectivity_obj, color=color,
                                                                                     style=style, count=row_count)
                        if not selectivity_obj.fortified_table.empty:
                            table["table_rows"] = selectivity_obj.fortified_table
                            selectivity_tables.append(table)

                        if multiple_analyte:
                            all_units[selectivity_obj.analyte] += [selectivity_obj.conc_unit]
                        else:
                            all_units.append(selectivity_obj.conc_unit)
                        if selectivity_obj.valid_data:
                            selectivity_obj.validate_table_data()
                            all_selectivity_info.append(selectivity_obj.all_selectivity_info)
                            overall_result += selectivity_obj.check_final_selectivity(all_selectivity_info,
                                                                                      len(selectivity_tables), index)
                    except:
                        message = error_messages.get_message("data_error")
                        overall_result.append(
                            red(message, None, selectivity_obj.table_title, selectivity_obj.table_type))

                overall_result += selectivity_obj.result
        else:
            overall_result += [red('Could not find table', None, 'Selectivity', 'Selectivity')]

        if dil_integrity_tables:
            for table in dil_integrity_tables:
                dil_integrity_obj = DilutionIntegrity(table, analytes, template_type)
                dil_integrity_obj.validate_table_format()
                if dil_integrity_obj.valid_format:
                    try:
                        dil_integrity_obj.process_table()
                        excel_writer, notation_worksheet, row_count = write_to_excel(excel_writer, notation_worksheet,
                                                                                     dil_integrity_obj, color=color,
                                                                                     style=style, count=row_count)
                        analyte = dil_integrity_obj.analyte
                        if multiple_analyte:
                            all_units[analyte] += [dil_integrity_obj.conc_unit]
                        else:
                            all_units.append(dil_integrity_obj.conc_unit)
                        if dil_integrity_obj.valid_data:
                            dil_integrity_obj.validate_table_data()
                            if multiple_analyte:
                                dil_integrity_dict[analyte]["unit"] = dil_integrity_obj.conc_unit
                                dil_integrity_dict[analyte]["nominal"] = dil_integrity_obj.nominal_conc
                            else:
                                dil_integrity_dict["unit"] = dil_integrity_obj.conc_unit
                                dil_integrity_dict["nominal"] = dil_integrity_obj.nominal_conc
                    except:
                        message = error_messages.get_message("data_error")
                        overall_result.append(
                            red(message, None, dil_integrity_obj.table_title, dil_integrity_obj.table_type))

                overall_result += dil_integrity_obj.result
        else:
            overall_result += [yellow('Could not find table', None, 'Dilution Integrity', 'Dilution Integrity')]

        if dil_linear_tables:
            for table in dil_linear_tables:
                dil_lin_obj = DilutionLinearity(table, analytes, template_type)
                dil_lin_obj.validate_table_format()
                if dil_lin_obj.valid_format:
                    try:
                        dil_lin_obj.process_table()
                        excel_writer, notation_worksheet, row_count = write_to_excel(excel_writer, notation_worksheet,
                                                                                     dil_lin_obj, color=color,
                                                                                     style=style, count=row_count)
                        if multiple_analyte:
                            all_units[dil_lin_obj.analyte] += [dil_lin_obj.conc_unit]
                        else:
                            all_units.append(dil_lin_obj.conc_unit)
                        if dil_lin_obj.valid_data:
                            dil_lin_obj.validate_table_data()
                    except:
                        message = error_messages.get_message("data_error")
                        overall_result.append(red(message, None, dil_lin_obj.table_title, dil_lin_obj.table_type))

                overall_result += dil_lin_obj.result
        else:
            if 'lm' in report_type:
                overall_result += [red('Could not find table', None, 'Dilution Linearity', 'Dilution Linearity')]

        if recovery_tables:
            for table in recovery_tables:
                recovery_obj = ExtractionRecovery(table, template_type)
                recovery_obj.validate_table_format()
                if recovery_obj.valid_format:
                    try:
                        recovery_obj.process_table()
                        excel_writer, notation_worksheet, row_count = write_to_excel(excel_writer, notation_worksheet,
                                                                                     recovery_obj, color=color,
                                                                                     style=style, count=row_count)
                        if recovery_obj.valid_data:
                            recovery_obj.validate_table_data()
                    except:
                        message = error_messages.get_message("data_error")
                        overall_result.append(red(message, None, recovery_obj.table_title, recovery_obj.table_type))

                overall_result += recovery_obj.result
        else:
            overall_result += [yellow('Could not find table', None, 'Recovery', 'Recovery')]

        if matrix_effect_tables:
            for table in matrix_effect_tables:
                matrix_obj = MatrixEffect(table, analytes, template_type)
                matrix_obj.validate_table_format()
                if matrix_obj.valid_format:
                    try:
                        matrix_obj.process_table()
                        excel_writer, notation_worksheet, row_count = write_to_excel(excel_writer, notation_worksheet,
                                                                                     matrix_obj, color=color,
                                                                                     style=style, count=row_count)
                        if multiple_analyte:
                            all_units[matrix_obj.analyte] += [matrix_obj.conc_unit]
                        else:
                            all_units.append(matrix_obj.conc_unit)
                        if matrix_obj.valid_data:
                            matrix_obj.validate_table_data()
                    except:
                        message = error_messages.get_message("data_error")
                        overall_result.append(red(message, None, matrix_obj.table_title, matrix_obj.table_type))

                overall_result += matrix_obj.result
        else:
            overall_result += [yellow('Could not find table', None, 'Matrix Effect', 'Matrix Effect')]

        dummy_df = pd.DataFrame()
        if val_summary_tables:
            if not multiple_analyte:
                for table in val_summary_tables:
                    table_df = table["table_rows"]
                    dummy_df = pd.concat([dummy_df, table_df]).reset_index(drop=True)

                val_summary_table = val_summary_tables[0]
                val_summary_table["table_rows"] = dummy_df
                val_summary_tables = [val_summary_table]

            for parsed_tables in val_summary_tables:
                val_summary_obj = ValidationSummary(parsed_tables, possible_anticoagulants, possible_sample_matrices,
                                                    possible_sample_preps, possible_sample_species,
                                                    possible_regressions, dil_integrity_dict, ULOQ, LLOQ,
                                                    qc_header_values, cs_header_values, ap_stats, units, options,
                                                    multiple_analyte, analytes, template_type)
                val_summary_obj.process_table()
                excel_writer, notation_worksheet, row_count = write_to_excel(excel_writer, notation_worksheet,
                                                                             val_summary_obj, color=color, style=style,
                                                                             count=row_count)
                val_summary_obj.validate_table_data()
                overall_result += val_summary_obj.result

                if val_summary_obj.found_sample_matrix:
                    all_found_summary_values["found_sample_matrix"] += val_summary_obj.found_sample_matrix
                    all_found_summary_values["found_sample_matrix"] = list(
                        set(all_found_summary_values["found_sample_matrix"]))

                if val_summary_obj.found_anticoagulant:
                    all_found_summary_values["found_anticoagulant"] += val_summary_obj.found_anticoagulant
                    all_found_summary_values["found_anticoagulant"] = list(
                        set(all_found_summary_values["found_anticoagulant"]))

                if val_summary_obj.found_sample_prep:
                    all_found_summary_values["found_sample_prep"] += val_summary_obj.found_sample_prep
                    all_found_summary_values["found_sample_prep"] = list(
                        set(all_found_summary_values["found_sample_prep"]))

                if val_summary_obj.found_species:
                    all_found_summary_values["found_species"] += val_summary_obj.found_species
                    all_found_summary_values["found_species"] = list(set(all_found_summary_values["found_species"]))

                if val_summary_obj.found_regression:
                    all_found_summary_values["found_regression_weighting"].append(val_summary_obj.found_regression)
                    all_found_summary_values["found_regression_weighting"] = list(
                        set(all_found_summary_values["found_regression_weighting"]))

                if all_found_summary_values["found_LLOQ"] == 0:
                    all_found_summary_values["found_LLOQ"] = val_summary_obj.found_LLOQ

                if all_found_summary_values["found_ULOQ"] == 0:
                    all_found_summary_values["found_ULOQ"] = val_summary_obj.found_ULOQ
                if all_found_summary_values["found_freeze_thaws"] == 0:
                    all_found_summary_values["found_freeze_thaws"] = val_summary_obj.found_freeze_thaws
                if all_found_summary_values["found_short_stability"] == 0:
                    all_found_summary_values["found_short_stability"] = val_summary_obj.found_short_stability
                if all_found_summary_values["found_long_stability"] == 0:
                    all_found_summary_values["found_long_stability"] = val_summary_obj.found_long_stability
                if all_found_summary_values["found_run_length"] == 0:
                    all_found_summary_values["found_run_length"] = val_summary_obj.found_run_length
                if val_summary_obj.found_inter_precisions:
                    all_found_summary_values["found_inter_precisions"] += val_summary_obj.found_inter_precisions
                if val_summary_obj.found_intra_precisions:
                    all_found_summary_values["found_intra_precisions"] += val_summary_obj.found_intra_precisions
                if val_summary_obj.found_inter_accuracies:
                    all_found_summary_values["found_inter_accuracies"] += val_summary_obj.found_inter_accuracies
                if val_summary_obj.found_intra_accuracies:
                    all_found_summary_values["found_intra_accuracies"] += val_summary_obj.found_intra_accuracies
                if val_summary_obj.found_cs_concentrations:
                    all_found_summary_values["found_cs_concentrations"] += val_summary_obj.found_cs_concentrations
                    all_found_summary_values["found_cs_concentrations"] = list(
                        set(all_found_summary_values["found_cs_concentrations"]))
                if val_summary_obj.found_qc_concentrations:
                    all_found_summary_values["found_qc_concentrations"] += val_summary_obj.found_qc_concentrations
                    all_found_summary_values["found_qc_concentrations"] = list(
                        set(all_found_summary_values["found_qc_concentrations"]))
        else:
            overall_result += [red('Could not find table', None, 'Validation Summary', 'Validation Summary')]

    if not batch_performance_table.empty:

        overall_result += batch_performance_obj.batch_summary_check_accept_reject_ratio(total_batches,
                                                                                        pass_failed_batches)
        overall_result += batch_performance_obj.batch_summary_check_accept_reject_sequential(pass_failed_batches)
        analysts = batch_performance_obj.get_analyst_batch()
        try:
            run_length = options["found_run_length"]
        except:
            run_length = 0

        if run_length == 0:
            run_length = all_found_summary_values["found_run_length"]

        if run_length > 0:
            overall_result += batch_performance_obj.check_validated_run_length(run_length)

    if analysts:
        for analyst in analysts:
            analyst_name = analyst['analyst']
            failed = analyst['failed']
            total = analyst['total']
            analytes_batch_str = analyst.get('analytes_batch_str', "")
            analytes_batch_str = f"({analytes_batch_str}) " if analytes_batch_str != "" else ""
            if total > 0:
                percentage = round(100 * (failed / total))
                if percentage > 50:
                    overall_result += [red(f"{percentage}% of the verifiable assays performed by {analyst_name} "
                                           f"{analytes_batch_str}failed to meet acceptance criteria", None, 'Analysts',
                                           'Analysts')]
                elif percentage > 25:
                    overall_result += [yellow(f"{percentage}% of the verifiable assays performed by {analyst_name} "
                                              f"{analytes_batch_str}failed to meet acceptance criteria", None,
                                              'Analysts', 'Analysts')]
                elif percentage >= 0:
                    overall_result += [green(f"{percentage}% of the verifiable assays performed by {analyst_name} "
                                             f"{analytes_batch_str}failed to meet acceptance criteria", None,
                                             'Analysts', 'Analysts')]

    if tables:
        overall_result += check_across_titles(tables, possible_anticoagulants, possible_sample_matrices,
                                              possible_sample_preps, possible_sample_species, options)

        if multiple_analyte:
            overall_result += check_multi_analyte_across_title(tables, analytes)

    if not overall_result:
        overall_result += [red('RedThread was unable to analyze this visualization.', None, 'Report', 'Report')]

    try:
        mv_cs = options['found_cs_concentrations']

        missing_values = []
        for val in cs_header_values:
            if val not in mv_cs:
                missing_values.append(val)

        if missing_values:
            missing_values = str(missing_values).replace('[', '').replace(']', '').replace("Decimal(", "").replace("[",
                                                                                                                   "").replace(
                "]", "").replace("'", "").replace(")", "")
            if units:
                missing_values = missing_values.replace(',', " " + units + ",")
                missing_values = missing_values + " " + units
        overall_result += [
            yellow('The following concentrations were not validated in the Method Validation: ' + missing_values, None,
                   'Concentration Samples', 'Calibration Curve Samples')]
    except:
        pass

    overall_result = remove_duplicate_results(overall_result)

    overall_result += cross_check_unit(units, all_units, multiple_analyte, analytes)

    return overall_result, all_found_summary_values, excel_writer, notation_worksheet, row_count
