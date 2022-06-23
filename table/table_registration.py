import collections
from collections import defaultdict

import funcy
import pandas as pd

from checks import utils
from checks.check_result import yellow, green
from checks.utils import *


class FlagProperties:
    def __init__(self, table_type):
        self.table_type = table_type
        self.file_path = "config file/table_checks.ini"
        self.msg_dict = self.get_table_dict()

    def get_table_dict(self):
        """
        Read the file passed as parameter as a properties file.
        """
        config = configparser.RawConfigParser()
        config.read(self.file_path)
        data = config[self.table_type]
        return data

    def get_message(self, key, **kwargs):
        """
        **kwargs
        Batch: batch number
        cycle: cycle number
        reported_value: table value
        calc_value: calculated value
        threshold: test threshold
        added_str: string for batch failed/passed
        column: table column
        """
        kwargs = [('{' + key + '}', kwargs[key]) for key in kwargs]
        message = self.msg_dict[key]
        try:
            for item in kwargs:
                message = message.replace(item[0], str(item[1]))
        except Exception as e:
            pass

        return message


class Flags:
    found_sd = False
    found_re = False
    found_te = False
    found_cv = False
    found_ar = False
    found_n = False
    found_mean = False
    valid_format = True
    valid_data = True
    error_flag = False


class Table(Flags):
    def __init__(self, parsed_table, template_type):
        self.missing_col = None
        self.template_type = template_type
        self.fortified_table = pd.DataFrame()
        self.static_df = pd.DataFrame()
        self.final_static_df = pd.DataFrame()
        self.data_df = pd.DataFrame()
        self.result = list()
        self.prev_col = list()
        self.error_messages = FlagProperties("Errors")
        self.threshold_values = FlagProperties("Thresholds")
        self.valid_difference_values = FlagProperties("ValidDifferences")
        self.table_title = parsed_table["table_title"]
        self.table_type = parsed_table["table_type"]
        self.analysis_type = parsed_table["analysis_type"]
        self.tb_title = f'{parsed_table["table_title"]} {parsed_table["tb_title"]}'
        self.table_subtype = parsed_table["table_subtype"].lower()
        self.table_df = parsed_table["table_rows"]
        self.analyte = parsed_table.get("analyte", None)
        self.template = parsed_table.get("template", "")
        if self.template_type == "ppd":
            self.preprocess_ppd_table()
        self.cleaning_table()
        self.conc_unit = utils.find_units(list(self.table_df.columns))

    def cleaning_table(self):
        self.table_df = utils.remove_header(self.table_df, self.template_type)
        if self.table_type == "Blood Stability" or self.table_type == "Interference Test":
            self.table_df = utils.remove_footer_sb_table(self.table_df)
        else:
            self.table_df = utils.remove_footer(self.table_df)
        if self.table_type != "Summary":
            self.table_df, self.missing_col = utils.format_table_data_frame(self.table_df, self.table_type)

    def preprocess_ppd_table(self):
        self.table_df, self.prev_col = utils.ppd_merge_split_header(self.table_df)
        self.remove_consecutive_duplicate_column()
        if self.table_type == "Accuracy and Precision":
            self.table_df.columns = [str(x).replace(y, "").strip() for x, y in
                                     zip(self.table_df.columns, self.prev_col)]
            columns = [x for x in self.table_df.columns if str(x).strip() != ""]
            self.table_df = self.table_df[columns]

        elif self.table_type == "Selectivity":
            self.table_df, self.fortified_table = self.split_fortified_table(self.table_df)

        self.table_df = utils.ppd_merge_concentration_row_to_header(self.table_df)

    def check_required_column(self, required_col):
        missing_col = required_col.intersection(self.missing_col)
        if missing_col:
            message = self.error_messages.get_message("missing_col", col_names=", ".join(missing_col))
            self.result.append(red(message, None, self.table_title, self.table_type))
            self.valid_format = False

    def check_unite(self):
        if self.conc_unit == "":
            message = self.error_messages.get_message("missing_conc_unit")
            self.result.append(red(message, None, self.table_title, self.table_type))
            self.valid_format = False

    @staticmethod
    def split_fortified_table(table):
        fortified_table = pd.DataFrame()
        index = None
        try:
            for col in table.columns:
                fortified_table = table[table[col].str.lower().str.contains(r"fortified")]
                if not fortified_table.empty:
                    index = fortified_table.index.tolist()[0]
                    break
            if index:
                non_fortified_table = table[:index].reset_index(drop=True)
                fortified_table = table[index:].reset_index(drop=True)
                fortified_table.columns = fortified_table.iloc[0]
                fortified_table = fortified_table[1:]
            else:
                non_fortified_table = table
        except AttributeError:
            non_fortified_table = table
        return non_fortified_table, fortified_table

    def remove_consecutive_duplicate_column(self):
        duplicate_columns = utils.get_duplicate_columns(self.table_df)
        self.table_df = self.table_df.loc[:, ~duplicate_columns]
        self.prev_col = [val for i, val in enumerate(self.prev_col) if not duplicate_columns[i]]
        self.prev_col = [x for x in self.prev_col if x is not None]


class Stability(Table, Flags):
    def __init__(self, parsed_table, analytes, template_type):
        Table.__init__(self, parsed_table, template_type)
        self.stability_type = None
        self.orientation = ""
        self.found_duplicate_run = False
        self.hemolysis_effect = dict()
        self.analyte = get_analyte(analytes, self.tb_title)
        self.sample_count = Decimal(3)

    def process_table(self):
        self.stability_type = utils.get_stability_type(self.tb_title)
        self.orientation = utils.check_stability_orientation(self.table_df)
        calc_static_run_df = pd.DataFrame()
        final_static_run_df = pd.DataFrame()
        final_static_df = pd.DataFrame()
        final_df = pd.DataFrame()
        static_df = pd.DataFrame()
        try:
            if self.orientation == "VR":
                tables, duplicate_id = utils.split_run_df(self.table_df)
                for i in range(len(tables)):
                    table = tables[i]["table"]
                    data, static_df = split_static_table(table)
                    static_df = utils.drop_blank_col(static_df)
                    if i == len(tables) - 1:
                        try:
                            static_df, final_static_df = split_static_table(static_df, "two")
                        except:
                            pass
                    table, calc_static, self.found_re, found_cycle = format_df_calc_stats(data)
                    table = utils.fill_val(table)
                    final_df = pd.concat([final_df, table]).reset_index(drop=True)
                    if not static_df.empty:
                        calc_static_run_df = pd.concat([calc_static_run_df, calc_static]).reset_index(drop=True)

                        # Static Dataframe
                        static_df = utils.process_static_df(static_df, self.table_type)
                        static_df.insert(loc=0, column="run_id", value=table["run_id"][0])
                        final_static_run_df = pd.concat([final_static_run_df, static_df]).reset_index(drop=True)

            elif self.orientation == "HZ":
                final_df = pd.DataFrame(columns=["cycle_id", "column", "nominal", "conc", "re"])
                data_df, static_df = utils.split_static_table(self.table_df)
                data_df = utils.fill_val(data_df)
                static_df, missing_col = utils.format_table_data_frame(static_df, self.table_type)
                static_df = utils.drop_blank_col(static_df)
                static_col = static_df.iloc[:, 0]
                column = list(data_df.columns)

                self.found_re = utils.find_re_column(data_df)
                cycles = sorted(list(utils.get_cycles(column)))
                count = 0
                for cycle in cycles:
                    col = list(filter(lambda x: cycle in str(x), column))
                    for col in col:
                        index = data_df.columns.get_loc(col)
                        data_df.columns.values[index + 1] = f'{count} RE {cycle}'
                        count += 1

                for cycle in cycles:
                    df = data_df.filter(regex=cycle)
                    column = list(df.columns)
                    c = 0
                    conc = utils.find_nominal_conc(column)
                    for i, key in enumerate(conc):
                        run_df = pd.DataFrame(columns=["cycle_id", "column", "nominal", "conc", "re"])
                        values = df[column[c]].apply(utils.parse_decimal)
                        nominal = conc[key]
                        run_df["conc"] = values
                        run_df["column"] = key
                        run_df["nominal"] = nominal
                        run_df["cycle_id"] = cycle[0]
                        if self.found_re:
                            run_df["re"] = df[column[c + 1]].apply(utils.parse_decimal)
                            c += 2
                        else:
                            c += 1
                        final_df = pd.concat([final_df, run_df]).reset_index(drop=True)

                        calc_static_run_df = pd.concat(
                            [calc_static_run_df, utils.build_static_df(values, nominal)]).reset_index(drop=True)

                    # static_df dataframe
                    sub_static_df = static_df.filter(regex=cycle)
                    sub_static_df.insert(loc=0, column="column", value=static_col)
                    sub_static_df = utils.process_static_df(sub_static_df, self.table_type)
                    sub_static_df.insert(loc=0, column="cycle_id", value=cycle[0])
                    final_static_run_df = pd.concat([final_static_run_df, sub_static_df]).reset_index(
                        drop=True)  # End static_df dataframe

                final_static_run_df.insert(loc=1, column="run_id", value=data_df["run_id"][0])
                final_df.insert(loc=1, column="run_id", value=data_df["run_id"][0])

            final_df["calc_re"] = utils.calculate_re(final_df["conc"], final_df["nominal"])
            if self.found_re:
                final_df["per_diff_re"] = utils.calculate_per_diff(final_df["re"], final_df["calc_re"])

            if not final_static_run_df.empty:
                self.static_df, self.found_sd = utils.concat_static_df(final_static_run_df, calc_static_run_df)
            self.data_df = final_df.fillna("")

            if not final_static_df.empty:
                final_static_df = utils.process_static_df(final_static_df, self.table_type)
                final_calc_static_df = pd.DataFrame()
                samples = final_df["column"].unique()
                for sample in samples:
                    values = final_df[final_df["column"] == sample]["conc"].to_list()
                    nominal = final_df[final_df["column"] == sample]["nominal"].to_list()[0]
                    final_calc_static_df = pd.concat(
                        [final_calc_static_df, utils.build_static_df(values, nominal)]).reset_index(drop=True)

                final_static_df, self.found_sd = utils.concat_static_df(final_static_df, final_calc_static_df)

                self.final_static_df = final_static_df

            self.found_n = find_count_column(static_df)
        except Exception as e:
            message = self.error_messages.get_message("data_error")
            self.result.append(red(message, None, self.table_title, self.table_type))
            self.valid_data = False

    def validate_table_data(self):
        validation_failed_flag_count = 0

        result = list()
        data_flag, stats_flag, final_stats_flag = self.get_error_dict()
        valid_count_diff = Decimal(self.valid_difference_values.get_message("difference_n"))
        valid_mean_diff = Decimal(self.valid_difference_values.get_message("difference_mean"))
        valid_re_cv_diff = Decimal(self.valid_difference_values.get_message("difference_re_cv"))
        failed_count = 0
        batch_failed = 0
        failed_sample = 0
        total_count = 0
        total_sample_count = 0
        failed_column = False
        first = False
        stability_type = self.stability_type
        if stability_type == "FT":
            added_str = "F/T"
        elif stability_type == "FR":
            added_str = "Frozen"
        elif stability_type == "LT":
            added_str = "Long-term"
        else:
            added_str = ""

        threshold = Decimal(self.threshold_values.get_message("re_cv_threshold"))
        if "lm" in self.analysis_type:
            test_threshold = threshold + Decimal(5)
        else:
            test_threshold = threshold

        data_df = self.data_df
        batches = data_df["run_id"]
        all_columns = data_df["column"]
        all_conc = data_df["conc"]
        all_calc_res = data_df["calc_re"]
        if self.found_re:
            reported_res = data_df["re"]
            per_diff_res = data_df["per_diff_re"]

        if self.orientation == "HZ":
            cycles_ids = data_df["cycle_id"]
        previous_batch = batches[0]
        previous_column = all_columns[0]
        self.hemolysis_effect = dict.fromkeys(list(batches))
        stability_failed_batch_info_dict = {batch: False for batch in set(batches)}

        for index, column in enumerate(all_columns):
            batch = batches[index]
            column = all_columns[index]
            conc = all_conc[index]
            calc_re = all_calc_res[index]

            # set parameter for flag
            if self.orientation == "VR":
                failed_batch_para = {"batch": previous_batch, "added_str": added_str}
                threshold_red_para = {"batch": batch, "column": column, "conc": conc, "threshold": test_threshold}
                threshold_yellow_para = {"batch": batch, "column": column, "conc": conc,
                                         "threshold": test_threshold - Decimal(5)}
            elif self.orientation == "HZ":
                cycle = cycles_ids[index]
                failed_batch_para = {"batch": previous_batch, "added_str": added_str, "cycle": cycle}
                threshold_red_para = {"batch": batch, "column": column, "conc": conc, "threshold": test_threshold,
                                      "cycle": cycle}
                threshold_yellow_para = {"batch": batch, "column": column, "conc": conc,
                                         "threshold": test_threshold - Decimal(5), "cycle": cycle}

            if self.found_re:
                reported_re = reported_res[index]
                per_diff_re = per_diff_res[index]

                if self.orientation == "VR":
                    re_para = {"batch": batch, "column": column, "reported_value": reported_re,
                               "calc_value": utils.format_value(calc_re), "conc": conc}

                elif self.orientation == "HZ":
                    re_para = {"batch": batch, "column": column, "reported_value": reported_re,
                               "calc_value": utils.format_value(calc_re), "cycle": cycle, "conc": conc}

                if per_diff_re > valid_re_cv_diff:
                    if abs(calc_re) <= Decimal(1) and abs(reported_re) <= Decimal(1):
                        message = data_flag.get_message("re_rounding", **re_para)

                        result.append(yellow(message, None, self.table_title, self.table_type))
                    else:
                        message = data_flag.get_message("re_error", **re_para)
                        result.append(red(message, None, self.table_title, self.table_type))
                        if self.stability_type == "HM":
                            self.hemolysis_effect[batch] = True

            if stability_type is not None and added_str != "":
                if previous_batch == batch:
                    total_count += 1
                    if previous_column == column:
                        total_sample_count += 1
                        first = True

                if previous_column != column:
                    if failed_sample > (total_sample_count - (total_sample_count // 2)):
                        if first:
                            message = data_flag.get_message("failed_batch", **failed_batch_para)
                            result.append(red(message, None, self.table_title, self.table_type))
                            first = False
                        failed_column = True
                        batch_failed += 1
                        stability_failed_batch_info_dict[previous_batch] = True
                    previous_column = column
                    total_sample_count = 0
                    failed_sample = 0

                elif previous_batch != batch:
                    if not failed_column:
                        if failed_count > (total_count - (total_count * 2 / 3)):
                            message = data_flag.get_message("failed_batch", **failed_batch_para)

                            result.append(red(message, None, self.table_title, self.table_type))
                            batch_failed += 1
                            stability_failed_batch_info_dict[previous_batch] = True

                    total_count = 0
                    failed_count = 0
                    previous_batch = batch

            if str(calc_re) != "" and calc_re is not None:
                if abs(calc_re) > test_threshold:
                    message = data_flag.get_message("re", **threshold_red_para)
                    result.append(red(message, None, self.table_title, self.table_type))
                    failed_count += 1
                    failed_sample += 1
                    validation_failed_flag_count += 1
                    if self.stability_type == "HM":
                        self.hemolysis_effect[batch] = True

                elif (test_threshold - Decimal(5)) <= abs(calc_re) <= test_threshold:
                    message = data_flag.get_message("re", **threshold_yellow_para)
                    result.append(yellow(message, None, self.table_title, self.table_type))

        sample_count = self.sample_count
        static_df = self.static_df
        if not static_df.empty:
            all_columns = static_df["column"]
            reported_mean = static_df["mean"]
            reported_sd = static_df["sd"]
            reported_cv = static_df["cv"]
            reported_re = static_df["re"]

            calculated_mean = static_df["calc_mean"]
            calculated_sd = static_df["calc_sd"]
            calculated_cv = static_df["calc_cv"]
            calculated_re = static_df["calc_re"]
            calculated_n = static_df["calc_n"]

            per_diff_mean = static_df["per_diff_mean"]
            per_diff_sd = static_df["per_diff_sd"]
            per_diff_cv = static_df["per_diff_cv"]
            per_diff_re = static_df["per_diff_re"]
            batches = static_df["run_id"]

            if self.found_n:
                reported_n = static_df["n"]
                per_diff_n = static_df["per_diff_n"]

            if self.orientation == "HZ":
                cycles_ids = static_df["cycle_id"]

            for index, column in enumerate(all_columns):
                batch = batches[index]
                overall_reported_mean = reported_mean[index]
                overall_reported_sd = reported_sd[index]
                overall_reported_cv = reported_cv[index]
                overall_reported_re = reported_re[index]

                overall_clc_mean = calculated_mean[index]
                overall_clc_sd = calculated_sd[index]
                overall_clc_cv = calculated_cv[index]
                overall_clc_re = calculated_re[index]
                overall_clc_n = calculated_n[index]

                overall_per_diff_mean = per_diff_mean[index]
                overall_per_diff_sd = per_diff_sd[index]
                overall_per_diff_cv = per_diff_cv[index]
                overall_per_diff_re = per_diff_re[index]

                if self.orientation == "HZ":
                    cycle = cycles_ids[index]
                    count_error_para = {"batch": batch, "column": column, "reported_value": overall_reported_n,
                                        "calc_value": overall_clc_n, "cycle": cycle}
                    mean_para = {"batch": batch, "column": column, "reported_value": overall_reported_mean,
                                 "calc_value": utils.format_value(overall_clc_mean), "cycle": cycle}
                    cv_para = {"batch": batch, "column": column, "reported_value": overall_reported_cv,
                               "calc_value": utils.format_value(overall_clc_cv), "cycle": cycle}
                    re_para = {"batch": batch, "column": column, "reported_value": overall_reported_re,
                               "calc_value": utils.format_value(overall_clc_re), "cycle": cycle}
                    sd_para = {"batch": batch, "column": column, "reported_value": overall_reported_sd,
                               "calc_value": utils.format_value(overall_clc_sd), "cycle": cycle}
                    threshold_red_para = {"batch": batch, "column": column, "threshold": test_threshold, "cycle": cycle}
                    threshold_yellow_para = {"batch": batch, "column": column, "threshold": test_threshold - Decimal(5),
                                             "cycle": cycle}
                    failed_batch_para = {"batch": batch, "added_str": added_str, "cycle": cycle}

                elif self.orientation == "VR":
                    mean_para = {"batch": batch, "column": column, "reported_value": overall_reported_mean,
                                 "calc_value": utils.format_value(overall_clc_mean)}
                    cv_para = {"batch": batch, "column": column, "reported_value": overall_reported_cv,
                               "calc_value": utils.format_value(overall_clc_cv)}
                    re_para = {"batch": batch, "column": column, "reported_value": overall_reported_re,
                               "calc_value": utils.format_value(overall_clc_re)}
                    sd_para = {"batch": batch, "column": column, "reported_value": overall_reported_sd,
                               "calc_value": utils.format_value(overall_clc_sd)}
                    threshold_red_para = {"batch": batch, "column": column, "threshold": test_threshold}
                    threshold_yellow_para = {"batch": batch, "column": column, "threshold": test_threshold - Decimal(5)}
                    failed_batch_para = {"batch": batch, "added_str": added_str}

                if self.found_n:
                    overall_reported_n = reported_n[index]
                    overall_per_diff_n = per_diff_n[index]
                    count_error_para = {"batch": batch, "column": column, "reported_value": overall_reported_n,
                                        "calc_value": overall_clc_n}

                    if overall_per_diff_n > valid_count_diff:
                        message = stats_flag.get_message("count_error", **count_error_para)
                        result.append(red(message, None, self.table_title, self.table_type))

                if overall_per_diff_mean > valid_mean_diff:
                    if abs(overall_clc_mean) <= Decimal(1) and abs(overall_reported_mean) <= Decimal(1):
                        message = stats_flag.get_message("mean_rounding", **mean_para)
                        result.append(yellow(message, None, self.table_title, self.table_type))

                    else:
                        message = stats_flag.get_message("mean_error", **mean_para)
                        result.append(red(message, None, self.table_title, self.table_type))

                if overall_per_diff_cv > valid_re_cv_diff:
                    if abs(overall_clc_cv) <= Decimal(1) and abs(overall_reported_cv) <= Decimal(1):
                        message = stats_flag.get_message("cv_rounding", **cv_para)
                        result.append(yellow(message, None, self.table_title, self.table_type))

                    else:
                        message = stats_flag.get_message("cv_error", **cv_para)
                        result.append(red(message, None, self.table_title, self.table_type))

                if overall_per_diff_re > valid_re_cv_diff:
                    if abs(overall_clc_re) <= Decimal(1) and abs(overall_reported_re) <= Decimal(1):
                        message = stats_flag.get_message("re_rounding", **re_para)
                        result.append(yellow(message, None, self.table_title, self.table_type))

                    else:
                        message = stats_flag.get_message("re_error", **re_para)
                        result.append(red(message, None, self.table_title, self.table_type))

                if overall_per_diff_sd > valid_re_cv_diff:
                    if abs(overall_clc_sd) <= Decimal(1) and abs(overall_reported_sd) <= Decimal(1):
                        message = stats_flag.get_message("sd_rounding", **sd_para)
                        result.append(yellow(message, None, self.table_title, self.table_type))

                    else:
                        message = stats_flag.get_message("sd_error", **sd_para)
                        result.append(red(message, None, self.table_title, self.table_type))

                if abs(overall_clc_re) > test_threshold:
                    message = stats_flag.get_message("re", **threshold_red_para)
                    result.append(red(message, None, self.table_title, self.table_type))
                    validation_failed_flag_count += 1
                    stability_failed_batch_info_dict[batch] = True

                    if self.stability_type == "HM":
                        self.hemolysis_effect[batch] = True

                    if stability_type is not None and added_str != "":
                        message = stats_flag.get_message("failed_batch", **failed_batch_para)
                        result.append(red(message, None, self.table_title, self.table_type))
                        batch_failed += 1

                elif (test_threshold - Decimal(5)) <= abs(overall_clc_re) <= test_threshold:
                    message = stats_flag.get_message("re_yellow", **threshold_yellow_para)
                    result.append(yellow(message, None, self.table_title, self.table_type))

                if abs(overall_clc_cv) > test_threshold:
                    message = stats_flag.get_message("cv", **threshold_red_para)
                    result.append(red(message, None, self.table_title, self.table_type))
                    validation_failed_flag_count += 1
                    stability_failed_batch_info_dict[batch] = True

                    if stability_type is not None and added_str != "":
                        message = stats_flag.get_message("failed_batch", **failed_batch_para)
                        result.append(red(message, None, self.table_title, self.table_type))
                        batch_failed += 1
                    if self.stability_type == "HM":
                        self.hemolysis_effect[batch] = True

                elif (test_threshold - Decimal(5)) <= abs(overall_clc_cv) <= test_threshold:
                    message = stats_flag.get_message("cv_yellow", **threshold_yellow_para)
                    result.append(yellow(message, None, self.table_title, self.table_type))

                # if self.stability_type == "HM":
                required_count_para = {"batch": batch, "column": column, "sample_count": sample_count}
                if overall_clc_n < sample_count:
                    message = stats_flag.get_message("required_count", **required_count_para)
                    result.append(red(message, None, self.table_title, self.table_type))
                    batch_failed += 1
                    validation_failed_flag_count += 1
                    stability_failed_batch_info_dict[batch] = True

                    if self.stability_type == "HM":
                        self.hemolysis_effect[batch] = True

        if not self.final_static_df.empty:
            final_static = self.final_static_df
            all_column = final_static["column"]
            reported_mean = final_static["mean"]
            reported_cv = final_static["cv"]
            reported_re = final_static["re"]
            reported_n = final_static["n"]

            calculated_mean = final_static["calc_mean"]
            calculated_cv = final_static["calc_cv"]
            calculated_re = final_static["calc_re"]
            calculated_n = final_static["calc_n"]

            per_diff_mean = final_static["per_diff_mean"]
            per_diff_cv = final_static["per_diff_cv"]
            per_diff_re = final_static["per_diff_re"]
            per_diff_n = final_static["per_diff_n"]

            if self.found_sd:
                reported_sd = final_static["sd"]
                calculated_sd = final_static["calc_sd"]
                per_diff_sd = final_static["per_diff_sd"]

            for index, column in enumerate(all_column):
                overall_reported_mean = reported_mean[index]
                overall_reported_cv = reported_cv[index]
                overall_reported_re = reported_re[index]
                overall_reported_n = reported_n[index]

                overall_clc_mean = calculated_mean[index]
                overall_clc_cv = calculated_cv[index]
                overall_clc_re = calculated_re[index]
                overall_clc_n = calculated_n[index]

                overall_per_diff_mean = per_diff_mean[index]
                overall_per_diff_cv = per_diff_cv[index]
                overall_per_diff_re = per_diff_re[index]
                overall_per_diff_n = per_diff_n[index]

                count_error_para = {"column": column, "reported_value": overall_reported_n,
                                    "calc_value": utils.format_value(overall_clc_n)}

                mean_para = {"column": column, "reported_value": overall_reported_mean,
                             "calc_value": utils.format_value(overall_clc_mean)}

                cv_para = {"column": column, "reported_value": overall_reported_cv,
                           "calc_value": utils.format_value(overall_clc_cv)}

                re_para = {"column": column, "reported_value": overall_reported_re,
                           "calc_value": utils.format_value(overall_clc_re)}

                threshold_red_para = {"column": column, "threshold": test_threshold}
                threshold_yellow_para = {"column": column, "threshold": test_threshold - Decimal(5)}

                if self.found_sd:
                    overall_reported_sd = reported_sd[index]
                    overall_clc_sd = calculated_sd[index]
                    overall_per_diff_sd = per_diff_sd[index]

                    sd_para = {"column": column, "reported_value": overall_reported_sd,
                               "calc_value": utils.format_value(overall_clc_sd)}

                    if overall_per_diff_sd > valid_re_cv_diff:
                        if abs(overall_clc_sd) <= Decimal(1) and abs(overall_reported_sd) <= Decimal(1):
                            message = final_stats_flag.get_message("sd_rounding", **sd_para)
                            result.append(yellow(message, None, self.table_title, self.table_type))

                        else:
                            message = final_stats_flag.get_message("sd_error", **sd_para)
                            result.append(red(message, None, self.table_title, self.table_type))

                if overall_per_diff_n > valid_count_diff:
                    message = final_stats_flag.get_message("count_error", **count_error_para)
                    result.append(red(message, None, self.table_title, self.table_type))

                if overall_per_diff_mean > valid_mean_diff:
                    message = final_stats_flag.get_message("mean_rounding", **mean_para)

                    if abs(overall_clc_mean) <= Decimal(1) and abs(overall_reported_mean) <= Decimal(1):
                        result.append(yellow(message, None, self.table_title, self.table_type))

                    else:
                        message = final_stats_flag.get_message("mean_error", **mean_para)
                        result.append(red(message, None, self.table_title, self.table_type))

                if overall_per_diff_cv > valid_re_cv_diff:
                    if abs(overall_clc_cv) <= Decimal(1) and abs(overall_reported_cv) <= Decimal(1):
                        message = final_stats_flag.get_message("cv_rounding", **cv_para)
                        result.append(yellow(message, None, self.table_title, self.table_type))

                    else:
                        message = final_stats_flag.get_message("cv_error", **cv_para)
                        result.append(red(message, None, self.table_title, self.table_type))

                if overall_per_diff_re > valid_re_cv_diff:
                    if abs(overall_clc_re) <= Decimal(1) and abs(overall_reported_re) <= Decimal(1):
                        message = final_stats_flag.get_message("re_rounding", **re_para)
                        result.append(yellow(message, None, self.table_title, self.table_type))

                    else:
                        message = final_stats_flag.get_message("re_error", **re_para)
                        result.append(red(message, None, self.table_title, self.table_type))

                if abs(overall_clc_re) > test_threshold:
                    message = final_stats_flag.get_message("re", **threshold_red_para)
                    result.append(red(message, None, self.table_title, self.table_type))
                    validation_failed_flag_count += 1

                elif (test_threshold - Decimal(5)) <= abs(overall_clc_re) <= test_threshold:
                    message = final_stats_flag.get_message("re_yellow", **threshold_yellow_para)
                    result.append(yellow(message, None, self.table_title, self.table_type))

                if abs(overall_clc_cv) > test_threshold:
                    message = final_stats_flag.get_message("cv", **threshold_red_para)
                    result.append(red(message, None, self.table_title, self.table_type))
                    validation_failed_flag_count += 1

                elif (test_threshold - Decimal(5)) <= abs(overall_clc_cv) <= test_threshold:
                    message = final_stats_flag.get_message("cv_yellow", **threshold_yellow_para)
                    result.append(yellow(message, None, self.table_title, self.table_type))

        # if self.stability_type == "HM": # according to token SD-136 this method applied for all tables
        required_column_level = list(filter(
            lambda x: "low" in str(x).lower() or "lqc" in str(x).lower() or "high" in str(x).lower() or "hqc" in str(
                x).lower(), self.table_df.columns))
        if len(required_column_level) < 2:
            message = stats_flag.get_message("required_conc_level", batch=batches[0])
            result.append(red(message, None, self.table_title, self.table_type))
            self.hemolysis_effect = {key: True for (key, value) in self.hemolysis_effect.items()}
            batch_failed += 1

        if stability_type is not None and added_str != "":
            if batch_failed == 0:
                pass_batch_para = {"added_str": added_str}
                message = final_stats_flag.get_message("all_batch_pass", **pass_batch_para)
                result.append(green(message, None, self.table_title, self.table_type))

        if validation_failed_flag_count == 0:
            result.append(green('All values within nominal range', None, self.table_title, self.table_type))

        if stability_type == "HP" or stability_type == "PS":
            for batch, has_failed in stability_failed_batch_info_dict.items():
                if not has_failed:
                    if stability_type == "HP":
                        message = stats_flag.get_message('not_lipidemic_effect', batch=batch)
                    elif stability_type == "PS":
                        message = stats_flag.get_message('stability_experiment_pass', batch=batch)

                    result.append(green(message, None, self.table_title, self.table_type))

        if self.stability_type == "HM":
            for batch, effect in self.hemolysis_effect.items():
                if not effect:
                    message = data_flag.get_message("hemolysis_not_observed", batch=batch)
                    result.append(green(message, None, self.table_title, self.table_type))
        self.result += result

    def get_error_dict(self):
        if self.stability_type == "HM":
            data_flag = FlagProperties("hm_stability_data")
            stats_flag = FlagProperties("hm_stability_batch_stats")
            final_stats_flag = FlagProperties("hm_stability_stats")
        elif self.orientation == "HZ":
            data_flag = FlagProperties("cycle_stability_data")
            stats_flag = FlagProperties("cycle_stability_batch_stats")
            final_stats_flag = FlagProperties("stability_stats")
        else:
            data_flag = FlagProperties("stability_data")
            stats_flag = FlagProperties("stability_batch_stats")
            final_stats_flag = FlagProperties("stability_stats")

        return data_flag, stats_flag, final_stats_flag

    def validate_freeze_thaw_number(self, ft_cycle):
        if self.stability_type == "FT":
            match = False
            try:
                cycles_ids = self.data_df["cycle_id"].unique()
                cycles_id = max([utils.parse_int(x) for x in cycles_ids if utils.parse_int(x) is not None])
                if ft_cycle is not None and cycles_id == int(ft_cycle):
                    pass
                else:
                    self.result.append(
                        red(f"Contains the Freeze/Thaw Cycle {cycles_ids[-1]} instead of {ft_cycle}", None,
                            self.table_title, self.table_type))
            except KeyError:
                title_array = self.tb_title.lower().split()
                title_array = [str(x).replace("(", "").replace(")", "") for x in title_array]
                try:
                    cycle = [x for x in title_array if "cycle" in x.lower()][0]
                    index = title_array.index(cycle)
                    cycle = utils.parse_int(title_array[(index + 1)])
                    if parse_int(ft_cycle) == cycle and cycle:
                        pass
                    else:
                        self.result.append(
                            red(f"Contains the Freeze/Thaw Cycle {cycle} instead of {ft_cycle}", None, self.table_title,
                                self.table_type))
                except (IndexError,):
                    self.result.append(red("Not found Cycle Id in the table", None, self.table_title, self.table_type))

    def validate_table_format(self):
        error_messages = self.error_messages
        required_col = {"Run ID"}
        missing_col = required_col.intersection(self.missing_col)
        if missing_col:
            message = error_messages.get_message("missing_col", col_names=", ".join(missing_col))
            self.result.append(red(message, None, self.table_title, self.table_type))
            self.valid_format = False

        conc_units = utils.find_units(self.table_df.columns)
        if conc_units == "":
            message = error_messages.get_message("missing_conc_unit")
            self.result.append(red(message, None, self.table_title, self.table_type))
            self.valid_format = False


class Calibration(Table, Flags):
    def __init__(self, parsed_table, analytes, template_type):
        Table.__init__(self, parsed_table, template_type)
        self.LLOQ = 0
        self.ULOQ = 0
        self.te_threshold = Decimal(self.threshold_values.get_message("te_threshold"))
        self.pass_batches = []
        self.failed_batches = []
        self.total_batches = []
        self.nominal_val = []
        self.found_duplicate_run = False
        self.found_stats_re = False
        self.analyte = get_analyte(analytes, self.tb_title)

    def validate_table_format(self):
        error_massages = self.error_messages
        required_column = {"Run ID"}
        missing_col = required_column.intersection(self.missing_col)
        if missing_col:
            message = error_massages.get_message("missing_col", col_names=' '.join(missing_col))
            self.result.append(red(message, None, self.table_title, self.table_type))
            self.valid_format = False

        if self.conc_unit == "":
            message = error_massages.get_message("missing_conc_unit")
            self.result.append(red(message, None, self.table_title, self.table_type))
            self.valid_format = False

    def process_table(self):
        table = self.table_df
        static_df = pd.DataFrame()
        final_df = pd.DataFrame()
        calc_static_df = pd.DataFrame()
        nominal_conc = utils.find_nominal_conc(table.columns)
        try:
            table, static_df = utils.split_static_table(table)
            table = utils.fill_val(table)
            all_run_tables, duplicate_id = utils.split_run_df(table)

            for index, run_table in enumerate(all_run_tables):
                data_table = run_table["table"]
                data_table, self.found_re = format_calibration_table(data_table, nominal_conc)
                final_df = pd.concat([final_df, data_table]).reset_index(drop=True)

            final_df["conc"] = final_df["conc"].apply(utils.parse_decimal)
            final_df["calc_re"] = utils.calculate_re(final_df["conc"], final_df["nominal"])
            if self.found_re:
                final_df["per_diff_re"] = utils.calculate_per_diff(final_df["re"], final_df["calc_re"])

            self.nominal_val = list(final_df["nominal"].unique())
            self.ULOQ = max(self.nominal_val)
            self.LLOQ = min(self.nominal_val)

            # remove outlier
            if "lm" in self.analysis_type:
                threshold_1 = Decimal(25)
                threshold_2 = Decimal(20)
                dummy_df = final_df[
                    ((final_df["nominal"].eq(self.LLOQ)) & (final_df["calc_re"].apply(abs) <= threshold_1)) | (
                            (final_df["nominal"].eq(self.ULOQ)) & (final_df["calc_re"].apply(abs) <= threshold_1)) | (
                            final_df["calc_re"].apply(abs) <= threshold_2)].reset_index(drop=True)
            else:
                threshold_1 = Decimal(20)
                threshold_2 = Decimal(15)
                dummy_df = final_df[
                    ((final_df["nominal"].eq(self.LLOQ)) & (final_df["calc_re"].apply(abs) <= threshold_1)) | (
                            final_df["calc_re"].apply(abs) <= threshold_2)].reset_index(drop=True)

            # *value excluded from overall calculation
            if not static_df.empty:
                samples = list(final_df["column"].unique())
                header_values = []
                for sample in samples:
                    try:
                        header_values.append(final_df[final_df["column"] == sample]["nominal"].to_list()[0])
                        values = dummy_df[dummy_df["column"] == sample]["conc"].to_list()
                        nominal = dummy_df[dummy_df["column"] == sample]["nominal"].to_list()
                        calc_static_df = pd.concat(
                            [calc_static_df, utils.build_static_df(values, nominal[0])]).reset_index(drop=True)
                    except IndexError:
                        calc_static_df = pd.concat([calc_static_df, utils.build_static_df([], [])]).reset_index(
                            drop=True)

                static_df = utils.process_static_df(static_df, self.table_type)
                static_df, self.found_sd = utils.concat_static_df(static_df, calc_static_df)
                self.found_te = utils.find_te_column(static_df)
                if self.found_te:
                    static_df["te"] = static_df["te"].apply(utils.parse_decimal)
                    static_df["calc_te"] = abs(static_df["calc_re"]) + abs(static_df["calc_cv"])
                    static_df["per_diff_te"] = abs((static_df["calc_te"] - static_df["te"]) / static_df["calc_te"])

            static_df = static_df.replace("", np.nan).dropna(how='all', axis=1)
            static_df = static_df.replace("NA", np.nan)
            self.found_cv = utils.find_cv_column(static_df)
            self.found_n = find_count_column(static_df)
            self.static_df = static_df.fillna("")
            self.data_df = final_df.fillna("")
            self.found_stats_re = find_re_column(static_df)
        except Exception as e:
            message = self.error_messages.get_message("data_error")
            self.result.append(red(message, None, self.table_title, self.table_type))
            self.valid_data = False

    def get_error_dict(self):
        data_flag = FlagProperties("calibration_data")
        stats_flag = FlagProperties("calibration_stats")
        return data_flag, stats_flag

    def validate_table_data(self):
        result = []
        total_sample = 0
        failed_sample = 0
        valid_count_diff = Decimal(self.valid_difference_values.get_message("difference_n"))
        valid_mean_diff = Decimal(self.valid_difference_values.get_message("difference_mean"))
        valid_re_cv_diff = Decimal(self.valid_difference_values.get_message("difference_re_cv"))
        valid_te_diff = Decimal(self.valid_difference_values.get_message("difference_te"))
        failed_batches = []
        table = self.data_df
        static_df = self.static_df

        data_flag, stats_flag = self.get_error_dict()
        threshold = Decimal(self.threshold_values.get_message("re_cv_threshold"))
        if "lm" in self.analysis_type:
            threshold = threshold + Decimal(5)

        all_columns = table["column"]
        batches = table["run_id"]
        nominal_concentrations = table["nominal"]
        all_conc = table["conc"]
        calc_res = table["calc_re"]
        total_batches = list(set(batches))
        if self.found_re:
            reported_res = table["re"]
            per_diff_res = table["per_diff_re"]
        previous_batch = batches[0]
        for index in range(len(all_columns)):
            failed_batch_para = {"batch": previous_batch}
            batch = batches[index]
            nominal = nominal_concentrations[index]
            column = all_columns[index]
            conc = all_conc[index]
            calc_re = calc_res[index]

            if previous_batch != batch:
                if (failed_sample / total_sample) > 0.25:
                    failed_batches.append(previous_batch)
                    message = data_flag.get_message("failed_batch", **failed_batch_para)
                    result.append(red(message, None, self.table_title, self.table_type))

                failed_sample = 0
                total_sample = 0
                previous_batch = batch

            total_sample += 1

            if self.found_re:
                reported_re = reported_res[index]
                per_diff_re = per_diff_res[index]
                re_error = {"batch": batch, "column": column, "reported_value": reported_re,
                            "calc_value": utils.format_value(calc_re), "conc": conc}

                if str(per_diff_re) != "":
                    if per_diff_re > valid_re_cv_diff:
                        if abs(calc_re) <= Decimal(1) and abs(reported_re) <= Decimal(1):
                            message = data_flag.get_message("re_rounding", **re_error)
                            result.append(yellow(message, None, self.table_title, self.table_type))
                        else:
                            message = data_flag.get_message("re_error", **re_error)
                            result.append(red(message, None, self.table_title, self.table_type))

            if str(calc_re) != "":
                if "lm" in self.analysis_type:
                    if str(self.LLOQ) == str(nominal) or str(self.ULOQ) == str(nominal):
                        test_threshold = threshold + Decimal(5)
                    else:
                        test_threshold = threshold
                else:
                    if str(self.LLOQ) == str(nominal):
                        test_threshold = threshold + Decimal(5)
                    else:
                        test_threshold = threshold
                re_red = {"batch": batch, "threshold": test_threshold, "conc": conc, "column": column}
                re_yellow = {"batch": batch, "threshold": test_threshold - Decimal(5), "conc": conc, "column": column}

                if abs(calc_re) > test_threshold:
                    message = data_flag.get_message("re", **re_red)
                    result.append(red(message, None, self.table_title, self.table_type))
                    failed_sample += 1

                    if str(self.LLOQ) == str(nominal) or str(self.ULOQ) == str(nominal):
                        truncate_range_para = {"batch": batch}
                        message = data_flag.get_message("truncate_range", **truncate_range_para)
                        result.append(yellow(message, None, self.table_title, self.table_type))

                elif (test_threshold - Decimal(5)) <= abs(calc_re) <= test_threshold:
                    message = data_flag.get_message("re", **re_yellow)
                    result.append(yellow(message, None, self.table_title, self.table_type))

            else:
                if str(self.LLOQ) == str(nominal) or str(self.ULOQ) == str(nominal):
                    truncate_range_para = {"batch": batch}
                    message = data_flag.get_message("truncate_range", **truncate_range_para)
                    result.append(yellow(message, None, self.table_title, self.table_type))
                failed_sample += 1

        if (failed_sample / total_sample) > 0.25:
            message = data_flag.get_message("failed_batch", **failed_batch_para)
            result.append(red(message, None, self.table_title, self.table_type))
            failed_batches.append(previous_batch)

        if not static_df.empty:
            all_column = static_df["column"]
            reported_mean = static_df["mean"]

            calculated_mean = static_df["calc_mean"]
            calculated_n = static_df["calc_n"]

            per_diff_mean = static_df["per_diff_mean"]

            nominal_value = static_df["nominal"]

            if self.found_n:
                reported_n = static_df["n"]
                per_diff_n = static_df["per_diff_n"]

            if self.found_cv:
                reported_cv = static_df["cv"]
                per_diff_cv = static_df["per_diff_cv"]
                calculated_cv = static_df["calc_cv"]
            else:
                message = stats_flag.get_message("cv_not_found")
                result.append(yellow(message, None, self.table_title, self.table_type))

            if self.found_sd:
                reported_sd = static_df["sd"]
                calculated_sd = static_df["calc_sd"]
                per_diff_sd = static_df["per_diff_sd"]

            if self.found_te:
                calculated_tes = static_df["calc_te"]
                per_diff_te = static_df["per_diff_te"]
                reported_te = static_df["te"]

            if self.found_stats_re:
                reported_re = static_df["re"]
                per_diff_re = static_df["per_diff_re"]
                calculated_re = static_df["calc_re"]

            for index in range(len(all_column)):
                te_threshold = self.te_threshold
                column = all_column[index]
                overall_reported_mean = reported_mean[index]
                overall_clc_mean = calculated_mean[index]
                overall_clc_n = calculated_n[index]
                nominal = nominal_value[index]
                overall_per_diff_mean = per_diff_mean[index]

                if "lm" in self.analysis_type:
                    if str(self.LLOQ) == str(nominal) or str(self.ULOQ) == str(nominal):
                        test_threshold = threshold + Decimal(5)
                        te_threshold += Decimal(10)
                    else:
                        test_threshold = threshold
                else:
                    if str(self.LLOQ) == str(nominal):
                        test_threshold = threshold + Decimal(5)
                        te_threshold += Decimal(10)
                    else:
                        test_threshold = threshold

                red_threshold = {"column": column, "threshold": test_threshold}
                yellow_threshold = {"column": column, "threshold": test_threshold - Decimal(5)}

                if self.found_cv:
                    overall_per_diff_cv = per_diff_cv[index]
                    overall_reported_cv = reported_cv[index]
                    overall_clc_cv = calculated_cv[index]

                    cv_error_para = {"column": column, "reported_value": overall_reported_cv,
                                     "calc_value": utils.format_value(overall_clc_cv)}

                    if overall_per_diff_cv > valid_re_cv_diff:
                        if abs(overall_reported_cv) <= Decimal(1) and abs(overall_clc_cv) <= Decimal(1):
                            message = stats_flag.get_message("cv_rounding", **cv_error_para)
                            result.append(yellow(message, None, self.table_title, self.table_type))

                        else:
                            message = stats_flag.get_message("cv_error", **cv_error_para)
                            result.append(red(message, None, self.table_title, self.table_type))

                    if str(overall_clc_cv).strip() != "":
                        if abs(overall_clc_cv) > test_threshold:
                            message = stats_flag.get_message("cv", **red_threshold)
                            result.append(red(message, None, self.table_title, self.table_type))
                        elif (test_threshold - Decimal(5)) <= abs(overall_clc_cv) <= test_threshold:
                            message = stats_flag.get_message("cv", **yellow_threshold)
                            result.append(yellow(message, None, self.table_title, self.table_type))

                if self.found_sd:
                    overall_reported_sd = reported_sd[index]
                    overall_clc_sd = calculated_sd[index]
                    overall_per_diff_sd = per_diff_sd[index]
                    sd_error_para = {"column": column, "reported_value": overall_reported_sd,
                                     "calc_value": utils.format_value(overall_clc_sd)}

                    if overall_per_diff_sd > valid_re_cv_diff:
                        if abs(overall_clc_sd) <= Decimal(1) and abs(overall_reported_sd) <= Decimal(1):
                            message = stats_flag.get_message("sd_rounding", **sd_error_para)
                            result.append(yellow(message, None, self.table_title, self.table_type))

                        else:
                            message = stats_flag.get_message("sd_error", **sd_error_para)
                            result.append(red(message, None, self.table_title, self.table_type))

                if self.found_te:
                    calculated_te = calculated_tes[index]
                    overall_per_diff_te = per_diff_te[index]
                    overall_reported_te = reported_te[index]
                    te_para = {"column": column, "te_threshold": te_threshold}
                    te_error_para = {"column": column, "reported_value": overall_reported_te,
                                     "calc_value": utils.format_value(calculated_te)}

                    if overall_per_diff_te > valid_te_diff:
                        message = stats_flag.get_message("te_error", **te_error_para)
                        result.append(red(message, None, self.table_title, self.table_type))

                    if abs(calculated_te) > te_threshold:
                        message = stats_flag.get_message("te", **te_para)
                        result.append(red(message, None, self.table_title, self.table_type))

                if self.found_n:
                    overall_reported_n = reported_n[index]
                    overall_per_diff_n = per_diff_n[index]
                    count_error_para = {"column": column, "reported_value": overall_reported_n,
                                        "calc_value": overall_clc_n}

                    if overall_per_diff_n > valid_count_diff:
                        message = stats_flag.get_message("count_error", **count_error_para)
                        result.append(red(message, None, self.table_title, self.table_type))

                mean_error_para = {"column": column, "reported_value": overall_reported_mean,
                                   "calc_value": utils.format_value(overall_clc_mean)}

                if overall_per_diff_mean > valid_mean_diff:
                    if abs(overall_clc_mean) <= Decimal(1) and abs(overall_reported_mean) <= Decimal(1):
                        message = stats_flag.get_message("mean_rounding", **mean_error_para)
                        result.append(yellow(message, None, self.table_title, self.table_type))

                    else:
                        message = stats_flag.get_message("mean_error", **mean_error_para)
                        result.append(red(message, None, self.table_title, self.table_type))

                if self.found_stats_re:
                    overall_per_diff_re = per_diff_re[index]
                    overall_reported_re = reported_re[index]
                    overall_clc_re = calculated_re[index]

                    re_error_para = {"column": column, "reported_value": overall_reported_re,
                                     "calc_value": utils.format_value(overall_clc_re)}

                    if overall_per_diff_re > valid_re_cv_diff:
                        if abs(overall_reported_re) <= Decimal(1) and abs(overall_clc_re) <= Decimal(1):
                            message = stats_flag.get_message("re_rounding", **re_error_para)
                            result.append(yellow(message, None, self.table_title, self.table_type))

                        else:
                            message = stats_flag.get_message("re_error", **re_error_para)
                            result.append(red(message, None, self.table_title, self.table_type))

                    red_threshold = {"column": column, "threshold": test_threshold}
                    yellow_threshold = {"column": column, "threshold": test_threshold - Decimal(5)}

                    if abs(overall_clc_re) > test_threshold:
                        message = stats_flag.get_message("re", **red_threshold)
                        result.append(red(message, None, self.table_title, self.table_type))

                    elif (test_threshold - Decimal(5)) <= abs(overall_clc_re) <= test_threshold:
                        message = stats_flag.get_message("re", **yellow_threshold)
                        result.append(yellow(message, None, self.table_title, self.table_type))

                else:
                    message = stats_flag.get_message("re_not_found")
                    result.append(yellow(message, None, self.table_title, self.table_type))

        else:
            result.append(yellow("Could not find overall statistics", None, self.table_title, self.table_type))
        if len(result) == 0:
            result.append(green('All values within nominal range', None, self.table_title, self.table_type))
        self.result += result
        passed_batches = [x for x in total_batches if x not in failed_batches]
        self.failed_batches = failed_batches
        self.pass_batches = passed_batches
        self.total_batches = total_batches


class ExtractionRecovery(Table, Flags):
    def __init__(self, parsed_table, template_type):
        Table.__init__(self, parsed_table, template_type)
        self.replicate = 0
        self.recovery = []
        self.overall_recovery = []

    def validate_table_format(self):
        error_messages = self.error_messages
        required_col = {"Run ID", "%Recovery", "Overall %Recovery", "QC Level"}
        missing_col = required_col.intersection(self.missing_col)
        if missing_col:
            message = error_messages.get_message("missing_col", col_names=", ".join(missing_col))
            self.result.append(red(message, None, self.table_title, self.table_type))
            self.valid_format = False

    def process_table(self):
        table_df = self.table_df
        static_df = pd.DataFrame()
        calc_static_df = pd.DataFrame()

        recovery = table_df["recovery"].apply(utils.parse_decimal).to_list()
        recovery = [x for x in recovery if x is not None]
        self.recovery = recovery

        overall_recovery = table_df["overall_recovery"].apply(utils.parse_decimal).to_list()
        overall_recovery = [x for x in overall_recovery if x is not None]
        self.overall_recovery = overall_recovery

        try:
            run_id = table_df["run_id"].to_list()
            run_id = [x for x in run_id if str(x).strip() != ""]
            if len(run_id) == 1:
                run_id = run_id[0]
            column = list(table_df.columns)
            column = list(
                filter(lambda x: "peak" in str(x).lower() and "area" in str(x).lower() or "ratio" in str(x).lower(),
                       column))
            qc_tables = utils.split_qc_df(table_df)
            for qc_table in qc_tables:
                sub_calc_static_df = pd.DataFrame()
                table = qc_table["table"]
                table, static = utils.split_static_table(table)
                for col in column:
                    values = table[col].apply(utils.parse_decimal).to_list()
                    sub_calc_static_df = pd.concat([sub_calc_static_df, utils.build_static_df(values, 0)])

                sub_calc_static_df["qc_level"] = qc_table["level"]
                calc_static_df = pd.concat([calc_static_df, sub_calc_static_df]).reset_index(drop=True)

                static = utils.drop_blank_col(static)
                static = static.T
                static = utils.remove_header(static)
                static = utils.remove_footer(static)
                static, missing_cols = utils.format_table_data_frame(static, self.table_type)
                static.insert(loc=0, column="run_id", value=[run_id] * len(column))
                static.insert(loc=1, column="column", value=column)
                static = static.reset_index(drop=True)
                static_df = pd.concat([static_df, static]).reset_index(drop=True)

            static_df, self.found_sd = utils.concat_static_df(static_df, calc_static_df)
            self.static_df = static_df.fillna("")
            self.found_n = utils.find_count_column(static_df)
            self.found_mean = utils.find_mean_column(static_df)
            self.found_cv = utils.find_cv_column(static_df)
            self.found_re = utils.find_stats_re_column(static_df)
        except Exception as e:
            message = self.error_messages.get_message("data_error")
            self.result.append(red(message, None, self.table_title, self.table_type))
            self.valid_data = False

    def validate_table_data(self):
        result = []
        valid_count_diff = Decimal(self.valid_difference_values.get_message("difference_n"))
        valid_mean_diff = Decimal(self.valid_difference_values.get_message("difference_mean"))
        valid_re_cv_diff = Decimal(self.valid_difference_values.get_message("difference_re_cv"))
        valid_te_diff = Decimal(self.valid_difference_values.get_message("difference_te"))

        static_df = self.static_df
        test_threshold = Decimal(self.threshold_values.get_message("re_cv_threshold"))
        if "lm" in self.analysis_type:
            test_threshold += Decimal(5)
        if not static_df.empty:
            batches = static_df["run_id"]
            qc_levels = static_df["qc_level"]
            all_column = static_df["column"]

            reported_mean = static_df["mean"]
            calculated_mean = static_df["calc_mean"]
            per_diff_mean = static_df["per_diff_mean"]

            if self.found_cv:
                reported_cv = static_df["cv"]
                per_diff_cv = static_df["per_diff_cv"]
            calculated_cv = static_df["calc_cv"]

            if self.found_re:
                reported_re = static_df["re"]
                per_diff_re = static_df["per_diff_re"]
            calculated_re = static_df["calc_re"]

            if self.found_n:
                reported_n = static_df["n"]
                per_diff_n = static_df["per_diff_n"]
            calculated_n = static_df["calc_n"]

            if self.found_sd:
                reported_sd = static_df["sd"]
                per_diff_sd = static_df["per_diff_sd"]
            calculated_sd = static_df["calc_sd"]
            stats_flag = self.get_error_dict()
            for index, column in enumerate(all_column):
                batch = batches[index]
                level = qc_levels[index]
                overall_reported_mean = reported_mean[index]
                overall_clc_mean = calculated_mean[index]
                overall_clc_cv = calculated_cv[index]
                overall_clc_re = calculated_re[index]
                overall_clc_n = calculated_n[index]
                overall_clc_sd = calculated_sd[index]
                overall_per_diff_mean = per_diff_mean[index]

                mean_error_para = {"level": level, "column": column, "batch": batch,
                                   "reported_value": overall_reported_mean,
                                   "calc_value": utils.format_value(overall_clc_mean)}

                red_threshold = {"batch": batch, "column": column, "level": level, "threshold": test_threshold}
                yellow_threshold = {"batch": batch, "column": column, "level": level,
                                    "threshold": test_threshold - Decimal(5)}

                if self.found_cv:
                    overall_reported_cv = reported_cv[index]
                    overall_per_diff_cv = per_diff_cv[index]
                    cv_error_para = {"level": level, "column": column, "batch": batch,
                                     "reported_value": overall_reported_cv,
                                     "calc_value": utils.format_value(overall_clc_cv)}

                    if overall_per_diff_cv > valid_re_cv_diff:
                        if abs(overall_reported_cv) <= Decimal(1) and abs(overall_clc_cv) <= Decimal(1):
                            message = stats_flag.get_message("cv_rounding", **cv_error_para)
                            result.append(yellow(message, None, self.table_title, self.table_type))

                        else:
                            message = stats_flag.get_message("cv_error", **cv_error_para)
                            result.append(red(message, None, self.table_title, self.table_type))

                if self.found_re:
                    overall_reported_re = reported_re[index]
                    overall_per_diff_re = per_diff_re[index]
                    re_error_para = {"level": level, "column": column, "batch": batch,
                                     "reported_value": overall_reported_re,
                                     "calc_value": utils.format_value(overall_clc_re)}

                    if overall_per_diff_re > valid_re_cv_diff:
                        if abs(overall_reported_re) <= Decimal(1) and abs(overall_clc_re) <= Decimal(1):
                            message = stats_flag.get_message("re_rounding", **re_error_para)
                            result.append(yellow(message, None, self.table_title, self.table_type))

                        else:
                            message = stats_flag.get_message("re_error", **re_error_para)
                            result.append(red(message, None, self.table_title, self.table_type))

                if self.found_n:
                    overall_reported_n = reported_n[index]
                    overall_per_diff_n = per_diff_n[index]
                    count_error_para = {"level": level, "column": column, "batch": batch,
                                        "reported_value": overall_reported_n, "calc_value": overall_clc_n}
                    if overall_per_diff_n > Decimal(0.01):
                        message = stats_flag.get_message("count_error", **count_error_para)
                        result.append(red(message, None, self.table_title, self.table_type))

                if self.found_sd:
                    overall_reported_sd = reported_sd[index]
                    overall_per_diff_sd = per_diff_sd[index]
                    sd_error_para = {"level": level, "column": column, "batch": batch,
                                     "reported_value": overall_reported_sd,
                                     "calc_value": utils.format_value(overall_clc_sd)}
                    if overall_per_diff_sd > valid_re_cv_diff:
                        if abs(overall_clc_sd) <= Decimal(1) and abs(overall_reported_sd) <= Decimal(1):
                            message = stats_flag.get_message("sd_rounding", **sd_error_para)
                            result.append(yellow(message, None, self.table_title, self.table_type))

                        else:
                            message = stats_flag.get_message("sd_error", **sd_error_para)
                            result.append(red(message, None, self.table_title, self.table_type))

                if self.replicate != 0:
                    replicate_para = {"column": column, "batch": batch, "level": level, "replicate": self.replicate}
                    if overall_clc_n < self.replicate:
                        message = stats_flag.get_message("replicate", **replicate_para)
                        result.append(red(message, None, self.table_title, self.table_type))

                if overall_per_diff_mean > valid_mean_diff:
                    if abs(overall_clc_mean) <= Decimal(1) and abs(overall_reported_mean) <= Decimal(1):
                        message = stats_flag.get_message("mean_rounding", **mean_error_para)
                        result.append(yellow(message, None, self.table_title, self.table_type))

                    else:
                        message = stats_flag.get_message("mean_error", **mean_error_para)
                        result.append(red(message, None, self.table_title, self.table_type))

                try:
                    if abs(overall_clc_re) > test_threshold:
                        message = stats_flag.get_message("re", **red_threshold)
                        result.append(red(message, None, self.table_title, self.table_type))

                    elif (test_threshold - Decimal(5)) <= abs(overall_clc_re) <= test_threshold:
                        message = stats_flag.get_message("re", **yellow_threshold)
                        result.append(yellow(message, None, self.table_title, self.table_type))
                except:
                    pass

                if abs(overall_clc_cv) > test_threshold:
                    message = stats_flag.get_message("cv", **red_threshold)
                    result.append(red(message, None, self.table_title, self.table_type))
                elif (test_threshold - Decimal(5)) <= abs(overall_clc_cv) <= test_threshold:
                    message = stats_flag.get_message("cv", **yellow_threshold)
                    result.append(yellow(message, None, self.table_title, self.table_type))

            if len(self.recovery) != 0:
                calculated_recoveries = []
                c = 0
                qc_levels = set(qc_levels)
                for index, level in enumerate(qc_levels):
                    f_mean = calculated_mean[c]
                    s_mean = calculated_mean[c + 1]
                    reported_recovery = self.recovery[index]
                    batch = batches[c]
                    c += 2
                    calc_recovery = (f_mean / s_mean) * 100
                    calculated_recoveries.append(calc_recovery)
                    recovery_para = {"batch": batch, "level": level, "reported_value": reported_recovery,
                                     "calc_value": utils.format_value(calc_recovery)}
                    per_diff_recovery = (calc_recovery - reported_recovery) / calc_recovery
                    if per_diff_recovery > Decimal(0.1):
                        message = stats_flag.get_message("recovery", **recovery_para)
                        result.append(red(message, None, self.table_title, self.table_type))

                try:
                    calc_overall_recovery = statistics.mean(calculated_recoveries)
                    per_diff_recovery = (calc_overall_recovery - self.overall_recovery[0]) / calc_overall_recovery
                    overall_recovery_para = {"batch": batches[0], "reported_value": self.overall_recovery[0],
                                             "calc_value": utils.format_value(calc_overall_recovery)}
                    if per_diff_recovery > Decimal(0.1):
                        message = stats_flag.get_message("overall_recovery_error", **overall_recovery_para)
                        result.append(red(message, None, self.table_title, self.table_type))
                except:
                    pass
            if len(result) == 0:
                result.append(green(f"All values within acceptable range", None, self.table_title, self.table_type))
        else:
            pass
        self.result += result

    def get_error_dict(self):
        stats_flag = FlagProperties("recovery_stats")
        return stats_flag


class DilutionLinearity(Table, Flags):
    def __init__(self, parsed_table, analytes, template_type):
        Table.__init__(self, parsed_table, template_type)
        self.orientation = ""
        self.analyte = get_analyte(analytes, self.tb_title)
        self.replicate = Decimal(5)
        self.ac_upper_threshold = Decimal(self.threshold_values.get_message("ac_upper_threshold"))
        self.ac_lower_threshold = Decimal(self.threshold_values.get_message("ac_lower_threshold"))
        self.found_accuracy = False

    def validate_table_format(self):
        error_messages = self.error_messages
        self.orientation = self.check_orientation()
        if self.orientation == "HZ":
            required_col = {"Run ID"}
            missing_col = required_col.intersection(self.missing_col)
            if missing_col:
                message = error_messages.get_message("missing_col", col_names=", ".join(missing_col))
                self.result.append(red(message, None, self.table_title, self.table_type))
                self.valid_format = False
        else:
            required_column = {"Run ID", "Samples", "Nominal", "%Accuracy", "Concentration", "%RE"}
            missing_col = list(set(required_column).intersection(set(self.missing_col)))
            missing_ac_re_col = {"%Accuracy", "%RE"}.intersection(missing_col)
            if len(missing_ac_re_col) == 2:
                message = error_messages.get_message("missing_col", col_names=' '.join(missing_col))
                self.result.append(red(message, None, self.table_title, self.table_type))
                self.valid_format = False
            elif len(missing_ac_re_col) == 1 and len(missing_col) > 1:
                missing_col = [x for x in missing_col if x not in missing_ac_re_col]
                message = error_messages.get_message("missing_col", col_names=' '.join(missing_col))
                self.result.append(red(message, None, self.table_title, self.table_type))
                self.valid_format = False

            elif len(missing_ac_re_col) == 0 and missing_col:
                message = error_messages.get_message("missing_col", col_names=' '.join(missing_col))
                self.result.append(red(message, None, self.table_title, self.table_type))
                self.valid_format = False

        if self.orientation == "HZ":
            conc_units = utils.find_units(self.table_df.columns)
            if conc_units == "":
                message = error_messages.get_message("missing_conc_unit")
                self.result.append(red(message, None, self.table_title, self.table_type))
                self.valid_format = False

    def process_table(self):
        table = self.table_df
        try:
            self.found_re = utils.find_re_column(table)
            if self.orientation == "VR":
                table = utils.remove_extra_char(table)
                self.found_accuracy = True if "ac" in list(table.columns) else False
                self.found_re = utils.find_re_column(table)
                self.found_cv = find_cv_column(table)
                table["calc_val"] = table["nominal"] / table["dilution"]
                if self.found_accuracy:
                    table["calc_ac"] = (table["conc"] / table["nominal"]) * 100
                    table["per_diff_ac"] = abs((table["calc_ac"] - table["ac"]) / table["calc_ac"])
                if self.found_re:
                    table["calc_re"] = ((table["conc"] - table["nominal"]) / table["nominal"]) * 100
                    table["per_diff_re"] = utils.calculate_per_diff(table["re"], table["calc_re"])

                self.data_df = table.fillna("").sort_values(["run_id", "dilution"]).reset_index(drop=True)

            elif self.orientation == "HZ":
                final_df = pd.DataFrame()
                static_df = pd.DataFrame()
                column = list(table.columns)
                nominal_conc = utils.find_nominal_conc(column)
                dilution_factor = utils.extract_dilution_factor(column)
                table_df, static_df = utils.split_static_table(table)
                calc_static_df = pd.DataFrame()

                while (len(nominal_conc) - 1) == len(dilution_factor):
                    dilution_factor.append(Decimal(1))
                count = 0
                for i, key in enumerate(nominal_conc):
                    run_df = pd.DataFrame(columns=["run_id", "column", "nominal", "dilution", "conc", "re"])
                    values = table_df[key].apply(utils.parse_decimal_2)
                    run_df["conc"] = values
                    run_df["run_id"] = table_df["run_id"]
                    run_df["column"] = key
                    nominal = nominal_conc[key]
                    run_df["nominal"] = nominal
                    run_df["dilution"] = dilution_factor[i]
                    if self.found_re:
                        try:
                            column_index = table_df.columns.get_loc(key)
                            table_df.columns.values[column_index + 1] = f"re_{count}"
                            run_df["re"] = table_df[f"re_{count}"].apply(utils.parse_decimal)
                            count += 1
                        except:
                            run_df["re"] = utils.parse_decimal("")
                    run_df = run_df.fillna("")
                    final_df = pd.concat([final_df, run_df]).reset_index(drop=True)
                    calc_static_df = pd.concat([calc_static_df, utils.build_static_df(values, nominal)]).reset_index(
                        drop=True)

                final_df = utils.fill_val(final_df)
                dummy_df = final_df.copy()
                dummy_df["conc"] = dummy_df["conc"].apply(utils.parse_decimal_2_1)
                final_df["calc_re"] = ((dummy_df["conc"] - dummy_df["nominal"]) / dummy_df["nominal"]) * 100
                if self.found_re:
                    final_df["per_diff_re"] = utils.calculate_per_diff(final_df["re"], final_df["calc_re"])

                # Static Calculation
                static_df = utils.process_static_df(static_df, self.table_type)
                static_df.insert(loc=0, column="run_id", value=final_df["run_id"][0])

                # Concat reported and calculated static dataframes
                final_static_df, self.found_sd = utils.concat_static_df(static_df, calc_static_df)
                self.data_df = final_df.fillna("")
                self.static_df = final_static_df.fillna("")
        except Exception as e:
            message = self.error_messages.get_message("data_error")
            self.result.append(red(message, None, self.table_title, self.table_type))
            self.valid_data = False

    def check_orientation(self):
        column = list(self.table_df.columns)
        concentration = utils.extract_nominal_conc(column)
        if concentration is None:
            return "VR"
        else:
            return "HZ"

    def validate_table_data(self):
        result = []
        valid_count_diff = Decimal(self.valid_difference_values.get_message("difference_n"))
        valid_mean_diff = Decimal(self.valid_difference_values.get_message("difference_mean"))
        valid_re_cv_diff = Decimal(self.valid_difference_values.get_message("difference_re_cv"))
        valid_te_diff = Decimal(self.valid_difference_values.get_message("difference_te"))
        threshold = Decimal(self.threshold_values.get_message("re_cv_threshold"))
        ac_upper_threshold = self.ac_upper_threshold
        ac_lower_threshold = self.ac_lower_threshold
        if "lm" in self.analysis_type:
            test_threshold = threshold + Decimal(5)
        else:
            test_threshold = threshold

        table = self.data_df
        if self.orientation == "VR":
            data_flag = self.get_error_dict()
            samples = table["samples"].to_list()
            dilutions = table["dilution"].to_list()
            batches = table["run_id"].to_list()
            if self.found_accuracy:
                all_accuracy = table["ac"]
                calculated_acs = table["calc_ac"]
                per_diff_acs = table["per_diff_ac"]
            elif self.found_re:
                reported_res = table["re"]
                calc_res = table["calc_re"]
                per_diff_res = table["per_diff_re"]
            count = 0
            batch = batches[0]
            sample = samples[0]
            pre_dilution = dilutions[0]
            for index, sample in enumerate(samples):
                dilution = dilutions[index]
                batch = batches[index]

                replicate_para = {"batch": batch, "sample": sample, "dilution": pre_dilution,
                                  "replicate": self.replicate}
                if pre_dilution != dilution:
                    if count < self.replicate:
                        message = data_flag.get_message("replicate", **replicate_para)
                        result.append(red(message, None, self.table_title, self.table_type))
                    count = 0
                    pre_dilution = dilution

                if self.found_accuracy:
                    accuracy = all_accuracy[index]
                    calc_ac = calculated_acs[index]
                    per_diff_ac = per_diff_acs[index]

                    ac_error_para = {"batch": batch, "sample": sample, "reported_value": accuracy,
                                     "calc_value": utils.format_value(calc_ac, decimal_point=1)}
                    ac_lower_para = {"batch": batch, "sample": sample, "threshold": ac_lower_threshold}
                    ac_upper_para = {"batch": batch, "sample": sample, "threshold": ac_upper_threshold}

                    if str(per_diff_ac).strip() != "":
                        if per_diff_ac > valid_re_cv_diff:
                            message = data_flag.get_message("ac_error", **ac_error_para)
                            result.append(red(message, None, self.table_title, self.table_type))

                    if str(calc_ac).strip() != "":
                        if abs(calc_ac) < ac_lower_threshold:
                            message = data_flag.get_message("ac_lower", **ac_lower_para)
                            result.append(red(message, None, self.table_title, self.table_type))

                        elif abs(calc_ac) > ac_upper_threshold:
                            message = data_flag.get_message("ac_upper", **ac_upper_para)
                            result.append(red(message, None, self.table_title, self.table_type))
                elif self.found_re:
                    reported_re = reported_res[index]
                    calc_re = calc_res[index]
                    per_diff_re = per_diff_res[index]
                    re_error_para = {"sample": sample, "batch": batch, "reported_value": reported_re,
                                     "calc_value": utils.format_value(calc_re)}
                    if str(per_diff_re) != "":
                        if per_diff_re > valid_re_cv_diff:
                            if abs(calc_re) <= Decimal(1) and abs(reported_re) <= Decimal(1):
                                message = data_flag.get_message("re_rounding", **re_error_para)
                                result.append(yellow(message, None, self.table_title, self.table_type))
                            else:
                                message = data_flag.get_message("re_error", **re_error_para)
                                result.append(red(message, None, self.table_title, self.table_type))

                    red_threshold_para = {"sample": sample, "batch": batch, "threshold": test_threshold}
                    yellow_threshold_para = {"sample": sample, "batch": batch, "threshold": test_threshold - Decimal(5)}
                    if str(calc_re).strip() != "":
                        if abs(calc_re) > test_threshold:
                            message = data_flag.get_message("re", **red_threshold_para)
                            result.append(red(message, None, self.table_title, self.table_type))

                        elif (test_threshold - Decimal(5)) <= abs(calc_re) <= test_threshold:
                            message = data_flag.get_message("re", **yellow_threshold_para)
                            result.append(yellow(message, None, self.table_title, self.table_type))

                count += 1

            if self.found_cv:
                pass
            else:
                message = data_flag.get_message("missing_cv")
                self.result.append(yellow(message, None, self.table_title, self.table_type))

            replicate_para = {"batch": batch, "sample": sample, "dilution": pre_dilution, "replicate": self.replicate}
            if count < self.replicate:
                message = data_flag.get_message("replicate", **replicate_para)
                result.append(red(message, None, self.table_title, self.table_type))

        elif self.orientation == "HZ":
            data_flag, stats_flag = self.get_error_dict()
            batches = table["run_id"]
            all_columns = table["column"]
            all_conc = table["conc"]
            all_calc_res = table["calc_re"]

            if self.found_re:
                reported_res = table["re"]
                per_diff_res = table["per_diff_re"]

            for index, column in enumerate(all_columns):
                batch = batches[index]
                conc = all_conc[index]
                column = all_columns[index]
                calc_re = all_calc_res[index]

                if self.found_re:
                    reported_re = reported_res[index]
                    per_diff_re = per_diff_res[index]
                    re_error_para = {"column": column, "batch": batch, "reported_value": reported_re,
                                     "calc_value": utils.format_value(calc_re), "conc": conc}
                    if str(per_diff_re) != "":
                        if per_diff_re > valid_re_cv_diff:
                            if abs(calc_re) <= Decimal(1) and abs(reported_re) <= Decimal(1):
                                message = data_flag.get_message("re_rounding", **re_error_para)
                                result.append(yellow(message, None, self.table_title, self.table_type))
                            else:
                                message = data_flag.get_message("re_error", **re_error_para)
                                result.append(red(message, None, self.table_title, self.table_type))

                red_threshold_para = {"column": column, "batch": batch, "threshold": test_threshold, "conc": conc}
                yellow_threshold_para = {"column": column, "batch": batch, "threshold": test_threshold - Decimal(5),
                                         "conc": conc}
                if str(calc_re).strip() != "":
                    if abs(calc_re) > test_threshold:
                        message = data_flag.get_message("re", **red_threshold_para)
                        result.append(red(message, None, self.table_title, self.table_type))

                    elif (test_threshold - Decimal(5)) <= abs(calc_re) <= test_threshold:
                        message = data_flag.get_message("re", **yellow_threshold_para)
                        result.append(yellow(message, None, self.table_title, self.table_type))

            static_df = self.static_df
            if not static_df.empty:
                all_column = static_df["column"]
                reported_mean = static_df["mean"]

                reported_cv = static_df["cv"]
                reported_re = static_df["re"]
                reported_n = static_df["n"]

                calculated_mean = static_df["calc_mean"]
                calculated_sd = static_df["calc_sd"]
                calculated_cv = static_df["calc_cv"]
                calculated_re = static_df["calc_re"]
                calculated_n = static_df["calc_n"]

                per_diff_mean = static_df["per_diff_mean"]

                per_diff_cv = static_df["per_diff_cv"]
                per_diff_re = static_df["per_diff_re"]
                per_diff_n = static_df["per_diff_n"]

                if self.found_sd:
                    reported_sd = static_df["sd"]
                    per_diff_sd = static_df["per_diff_sd"]

                for index, column in enumerate(all_column):

                    overall_reported_mean = reported_mean[index]
                    overall_reported_cv = reported_cv[index]
                    overall_reported_re = reported_re[index]
                    overall_reported_n = reported_n[index]

                    overall_clc_mean = calculated_mean[index]
                    overall_clc_cv = calculated_cv[index]
                    overall_clc_re = calculated_re[index]
                    overall_clc_n = calculated_n[index]

                    overall_per_diff_mean = per_diff_mean[index]
                    overall_per_diff_cv = per_diff_cv[index]
                    overall_per_diff_re = per_diff_re[index]
                    overall_per_diff_n = per_diff_n[index]

                    replicate_para = {"column": column, "replicate": self.replicate}
                    count_error_para = {"column": column, "reported_value": overall_reported_n,
                                        "calc_value": overall_clc_n}
                    mean_error_para = {"column": column, "reported_value": overall_reported_mean,
                                       "calc_value": utils.format_value(overall_clc_mean)}
                    re_error_para = {"column": column, "reported_value": overall_reported_re,
                                     "calc_value": utils.format_value(overall_clc_re)}
                    cv_error_para = {"column": column, "reported_value": overall_reported_cv,
                                     "calc_value": utils.format_value(overall_clc_cv)}
                    red_threshold_para = {"column": column, 'threshold': test_threshold}
                    yellow_threshold_para = {"column": column, 'threshold': test_threshold - Decimal(5)}

                    if self.found_sd:
                        overall_reported_sd = reported_sd[index]
                        overall_clc_sd = calculated_sd[index]
                        overall_per_diff_sd = per_diff_sd[index]

                        sd_error_para = {"column": column, "reported_value": overall_reported_sd,
                                         "calc_value": utils.format_value(overall_clc_sd)}
                        if overall_per_diff_sd > valid_re_cv_diff:
                            if abs(overall_clc_sd) <= Decimal(1) and abs(overall_reported_sd) <= Decimal(1):
                                message = stats_flag.get_message("sd_rounding", **sd_error_para)
                                result.append(yellow(message, None, self.table_title, self.table_type))

                            else:
                                message = stats_flag.get_message("sd_error", **sd_error_para)
                                result.append(red(message, None, self.table_title, self.table_type))

                    if overall_clc_n < self.replicate:
                        message = stats_flag.get_message("replicate", **replicate_para)
                        result.append(red(message, None, self.table_title, self.table_type))

                    if overall_per_diff_n > valid_count_diff:
                        message = stats_flag.get_message("count_error", **count_error_para)
                        result.append(red(message, None, self.table_title, self.table_type))

                    if str(overall_per_diff_mean).strip() != "":
                        if overall_per_diff_mean > valid_mean_diff:
                            if abs(overall_clc_mean) <= Decimal(1) and abs(overall_reported_mean <= Decimal(1)):
                                message = stats_flag.get_message("mean_rounding", **mean_error_para)
                                result.append(yellow(message, None, self.table_title, self.table_type))

                            else:
                                message = stats_flag.get_message("mean_error", **mean_error_para)
                                result.append(red(message, None, self.table_title, self.table_type))

                    if str(overall_per_diff_cv).strip() != "":
                        if overall_per_diff_cv > valid_re_cv_diff:
                            if abs(overall_reported_cv) <= Decimal(1) and abs(overall_clc_cv) <= Decimal(1):
                                message = stats_flag.get_message("cv_rounding", **cv_error_para)
                                result.append(yellow(message, None, self.table_title, self.table_type))

                            else:
                                message = stats_flag.get_message("cv_error", **cv_error_para)
                                result.append(red(message, None, self.table_title, self.table_type))

                    if str(overall_per_diff_re).strip() != "":
                        if overall_per_diff_re > valid_re_cv_diff:
                            if abs(overall_reported_re) <= Decimal(1) and abs(overall_clc_re) <= Decimal(1):
                                message = stats_flag.get_message("re_rounding", **re_error_para)
                                result.append(yellow(message, None, self.table_title, self.table_type))

                            else:
                                message = stats_flag.get_message("re_error", **re_error_para)
                                result.append(red(message, None, self.table_title, self.table_type))

                    if abs(overall_clc_re) > test_threshold:
                        message = stats_flag.get_message("re", **red_threshold_para)
                        result.append(red(message, None, self.table_title, self.table_type))

                    elif (test_threshold - Decimal(5)) <= abs(overall_clc_re) <= test_threshold:
                        message = stats_flag.get_message("re", **yellow_threshold_para)
                        result.append(yellow(message, None, self.table_title, self.table_type))

                    if abs(overall_clc_cv) > test_threshold:
                        message = stats_flag.get_message("cv", **red_threshold_para)
                        result.append(red(message, None, self.table_title, self.table_type))

                    elif (test_threshold - Decimal(5)) <= abs(overall_clc_cv) <= test_threshold:
                        message = stats_flag.get_message("cv", **yellow_threshold_para)
                        result.append(yellow(message, None, self.table_title, self.table_type))

        if len(result) == 0:
            result.append(green('All values within nominal range', None, self.table_title, self.table_type))

        self.result += result

    def get_error_dict(self):
        if self.orientation == "VR":
            data_flag = FlagProperties("vr_dilution_lin")
            return data_flag
        elif self.orientation == "HZ":
            data_flag = FlagProperties("hz_dilution_lin_data")
            stats_flag = FlagProperties("hz_dilution_lin_stats")
            return data_flag, stats_flag


class PPDAccuracyAndPrecision(Table, Flags):
    def __init__(self, parsed_table, analytes, template_type, LLOQ, ULOQ):
        Table.__init__(self, parsed_table, template_type)
        self.total_batches = []
        self.orientation = ""
        self.found_final_te = False
        self.LLOQ = LLOQ
        self.ULOQ = ULOQ
        self.te_threshold = Decimal(self.threshold_values.get_message("te_threshold"))
        self.analyte = get_analyte(analytes, self.tb_title)
        self.header_values = []
        self.failed_batch_numbers = []
        self.passed_batch_numbers = []
        self.intra_precisions = {}
        self.intra_accuracies = {}
        self.inter_precisions = {}
        self.inter_accuracies = {}
        self.inter_tes = {}
        self.intra_batch_info = {}
        self.inter_batch_info = {}
        if "sm" in self.analysis_type:
            self.required_levels = 4
            self.required_replicates = 5
            self.required_batches = 3
        elif "lm" in self.analysis_type:
            self.required_levels = 5
            self.required_replicates = 3
            self.required_batches = 6
        self.found_duplicate_run = False
        self.table_columns = [col for col in self.table_df.columns if "statistics" not in str(col).lower()]
        self.table_df = self.table_df[self.table_columns]
        self.rep_col = list(filter(lambda x: "rep" in str(x).lower(), self.table_columns))

    def validate_table_format(self):
        error_messages = self.error_messages
        required_col = {"Run ID", "Sample", "Mean", "%RE", "%CV", "SD", "N"}
        self.check_required_column(required_col)

        if self.valid_format:
            self.conc_unit = utils.find_units(self.table_df["sample"].tolist())
            if self.conc_unit == "":
                message = error_messages.get_message("missing_conc_unit")
                self.result.append(red(message, None, self.table_title, self.table_type))
                self.valid_format = False

    def flatten_table_data(self, table):
        final_df = pd.DataFrame(columns=["sample", "run_id", "column", "conc"])
        for _, row in table.iterrows():

            run_df = pd.DataFrame(columns=["sample", "run_id", "column", "conc"])
            conc_values = []
            column_val = []

            for col in self.rep_col:
                conc_values.append(row[col])
                column_val.append(col)

            run_df["conc"] = conc_values
            run_df["column"] = column_val
            run_df["sample"] = row["sample"]
            run_df["run_id"] = row["run_id"]
            final_df = pd.concat([final_df, run_df]).reset_index(drop=True)
        final_df["conc"] = final_df["conc"].apply(utils.parse_decimal)
        nominal_conc = utils.find_nominal_conc(final_df["sample"].unique())
        if nominal_conc:
            cols, vals = list(nominal_conc.keys()), list(nominal_conc.values())
            final_df = final_df.assign(nominal=np.where(final_df['sample'] == cols[0], vals[0], np.nan))
        return final_df

    def clean_table(self, table):
        column = self.table_columns
        rep_col_index = column.index(self.rep_col[0])
        rep_header_index = table.index[table[column[rep_col_index]] == column[rep_col_index]].tolist()
        rep_header_index = rep_header_index[0] if rep_header_index else None
        if rep_header_index is not None:
            table = table[rep_header_index + 1:].reset_index(drop=True)
        return table

    @staticmethod
    def fill_blank_samples(table):
        table["sample"] = table["sample"].replace("", np.nan)
        table.ffill(inplace=True)
        return table

    @staticmethod
    def recreate_sample(table):
        samples = table["sample"]
        sample_count = len(samples)
        sample = " ".join([x for x in samples if str(x).strip() != ""])
        table["sample"] = [sample] * sample_count
        return table

    def split_table_on_sample(self, table):
        prev_col = self.prev_col
        column = self.table_columns
        index = [prev_col.index(x) for x in prev_col if str(x).strip() != ""][0]
        next_table_indexes = table.index[table[column[index]] == prev_col[index]].tolist()
        sample_tables = []
        for index in next_table_indexes:
            first_table = table[:index].reset_index(drop=True)
            table = table[index + 1:].reset_index(drop=True)
            table = self.clean_table(table)
            sample_tables.append(first_table)

        sample_tables.append(table)
        return sample_tables

    def split_anova_stats(self, table):
        column = self.table_columns
        anova_stats = table[table[self.rep_col[0]].str.lower().str.contains(r"((?=.*inter)|(?=.*anova))(?=.*statistics)")]
        intra_stats = table[table[self.rep_col[0]].str.lower().str.contains(r"((?=.*intra)|(?=.*pooled))(?=.*statistics)")]
        table = table.drop(intra_stats.index.tolist() + anova_stats.index.tolist())
        anova_col = [x for x in column if x not in self.rep_col and str(x) != "run_id"]
        anova_stats = anova_stats[anova_col].reset_index(drop=True)
        return table, anova_stats

    def split_data_intra_stats(self, table):
        column = table.columns
        stats_col = list(filter(lambda x: "rep" not in str(x).lower(), column))
        table_col = ["sample", "run_id"] + self.rep_col
        data_table = table[table_col]
        stats_table = table[stats_col]
        return data_table, stats_table

    @staticmethod
    def calculate_inter_statics(table):
        stats_df = pd.DataFrame()
        samples = table["sample"].unique()
        nominal = table["nominal"].unique()[0]
        for sample in samples:
            values = table[table["sample"] == sample]["conc"]
            stats_df = pd.concat([stats_df, utils.build_static_df(values, nominal)]).reset_index(drop=True)
        return stats_df

    @staticmethod
    def calculate_intra_statics(table):
        stats_df = pd.DataFrame()
        run_ids = table["run_id"].unique()
        nominal = table["nominal"].unique()[0]
        for run_id in run_ids:
            values = table[table["run_id"] == run_id]["conc"]
            stats_df = pd.concat([stats_df, utils.build_static_df(values, nominal)]).reset_index(drop=True)
        return stats_df

    def process_table(self):
        table_df = self.table_df
        final_data_table = pd.DataFrame()
        calc_intra_static = pd.DataFrame()
        calc_inter_static = pd.DataFrame()
        final_inter_static_table = pd.DataFrame()
        final_intra_static_table = pd.DataFrame()
        try:
            sample_tables = self.split_table_on_sample(table_df)
            for sample_table in sample_tables:
                sample_table.replace("", np.nan, inplace=True)
                sample_table = sample_table.dropna(how='all')
                sample_table = sample_table.fillna("")

                intra_table, anova_stats = self.split_anova_stats(sample_table)
                intra_table = self.recreate_sample(intra_table)

                data_table, intra_stats = self.split_data_intra_stats(intra_table)
                anova_stats["sample"] = data_table["sample"][0] * anova_stats.shape[0]

                flatten_data_table = self.flatten_table_data(data_table)

                final_data_table = pd.concat([final_data_table, flatten_data_table]).reset_index(drop=True)
                calc_inter_static = pd.concat(
                    [calc_inter_static, self.calculate_inter_statics(flatten_data_table)]).reset_index(drop=True)

                calc_intra_static = pd.concat(
                    [calc_intra_static, self.calculate_intra_statics(flatten_data_table)]).reset_index(drop=True)

                final_intra_static_table = pd.concat([final_intra_static_table, intra_stats]).reset_index(drop=True)
                final_inter_static_table = pd.concat([final_inter_static_table, anova_stats]).reset_index(drop=True)

            final_intra_static_table, self.found_sd = utils.concat_static_df(final_intra_static_table,
                                                                             calc_intra_static)
            final_inter_static_table, self.found_sd = utils.concat_static_df(final_inter_static_table,
                                                                             calc_inter_static)

            final_data_table["calc_re"] = utils.calculate_re(final_data_table["conc"], final_data_table["nominal"])

            final_inter_static_table["calc_te"] = abs(final_inter_static_table["calc_re"]) + abs(
                final_inter_static_table["calc_cv"])

            final_intra_static_table["calc_te"] = abs(final_intra_static_table["calc_re"]) + abs(
                final_intra_static_table["calc_cv"])

            self.data_df = final_data_table.fillna("")
            self.static_df = final_intra_static_table.fillna("")
            self.final_static_df = final_inter_static_table.fillna("")
        except Exception as e:
            message = self.error_messages.get_message("data_error")
            self.result.append(red(message, None, self.table_title, self.table_type))
            self.valid_data = False

    def validate_table_data(self):
        data_flag, batch_stats_flag, stats_flag = self.get_error_dict()
        result = []

        valid_count_diff = Decimal(self.valid_difference_values.get_message("difference_n"))
        valid_mean_diff = Decimal(self.valid_difference_values.get_message("difference_mean"))
        valid_re_cv_diff = Decimal(self.valid_difference_values.get_message("difference_re_cv"))
        valid_te_diff = Decimal(self.valid_difference_values.get_message("difference_te"))

        lloq_intra_cvs = []
        uloq_intra_cvs = []
        lmh_intra_cvs = []

        lloq_intra_biases = []
        uloq_intra_biases = []
        lmh_intra_biases = []

        lmh_inter_cvs = []
        uloq_inter_cv = None
        lloq_inter_cv = None

        lloq_inter_bias = None
        uloq_inter_bias = None
        lmh_inter_biases = []

        lmh_inter_tes = []
        lloq_inter_te = None
        uloq_inter_te = None

        threshold = Decimal(self.threshold_values.get_message("re_cv_threshold"))
        if "lm" in self.analysis_type:
            threshold += Decimal(5)

        table = self.data_df
        if not table.empty:
            red_count = 0
            self.header_values = table["nominal"].unique()
            batches = table["run_id"].to_list()
            samples = table["sample"]
            header_names = list(table["sample"].unique())
            all_conc = table["conc"]
            all_calc_res = table["calc_re"]
            nominal_concentrations = table["nominal"]
            if self.found_re:
                reported_res = table["re"]
                per_diff_res = table["per_diff_re"]

            self.total_batches = list(set(batches))
            total_column = len(header_names)
            dict_val = [0]*len(header_names)
            yellow_columns = dict(zip(header_names, dict_val))
            red_columns = dict(zip(header_names, dict_val))
            previous_batch = batches[0]
            previous_sample = samples[0]
            failed_sample = 0
            total_sample = 0
            failed_column = 0

            failed_batch_para = {"batch": previous_batch, "col_level": self.get_column_level(previous_sample)}
            levels_para = {"batch": previous_batch, "replicate": self.required_levels}
            for index, sample in enumerate(samples):
                batch = batches[index]
                conc = all_conc[index]
                calc_re = all_calc_res[index]
                nominal = nominal_concentrations[index]

                failed_batch_para = {"batch": previous_batch, "col_level": self.get_column_level(previous_sample)}
                levels_para = {"batch": previous_batch, "required_levels": self.required_levels}

                if previous_sample != sample:
                    col_values = table[(table["run_id"] == previous_batch) & (table["sample"] == previous_sample)][
                        "conc"].replace("", np.nan).dropna()
                    if (failed_sample / total_sample) > 0.5:
                        failed_column += 1
                    if (total_sample - failed_sample) < self.required_replicates:
                        if self.table_subtype.lower() != "inter":
                            message = data_flag.get_message("failed_batch", **failed_batch_para)
                            result.append(red(message, None, self.table_title, self.table_type))
                        red_count += 1
                    elif col_values.empty:
                        red_count += 1

                    previous_sample = sample
                    failed_sample = 0
                    total_sample = 0

                if previous_batch != batch:
                    if total_column < self.required_levels:
                        if "lm" in self.analysis_type:
                            message = data_flag.get_message("lm_required_level", **levels_para)
                        else:
                            message = data_flag.get_message("sm_required_level", **levels_para)
                        result.append(red(message, None, self.table_title, self.table_type))
                        self.failed_batch_numbers.append(previous_batch)
                        red_count += 1

                    if (total_column - failed_column) < self.required_levels:
                        message = data_flag.get_message("required_pass_level", **levels_para)
                        result.append(red(message, None, self.table_title, self.table_type))
                        self.failed_batch_numbers.append(previous_batch)
                        red_count += 1

                    self.intra_batch_info[previous_batch] = red_count
                    red_count = 0
                    previous_batch = batch
                    failed_column = 0

                total_sample += 1

                if self.found_re:
                    reported_re = reported_res[index]
                    per_diff_re = per_diff_res[index]
                    re_error_para = {"batch": batch, "sample": sample, "reported_value": reported_re,
                                     "calc_value": utils.format_value(calc_re), "conc": conc}

                    if per_diff_re > valid_re_cv_diff:
                        if abs(calc_re) <= Decimal(1) and abs(reported_re) <= Decimal(1):
                            message = data_flag.get_message("re_rounding", **re_error_para)
                            result.append(yellow(message, None, self.table_title, self.table_type))
                        else:
                            message = data_flag.get_message("re_error", **re_error_para)
                            result.append(red(message, None, self.table_title, self.table_type))

                if str(calc_re) != "" and calc_re is not None:
                    if "lloq" in str(sample).lower() or "uloq" in str(
                            sample).lower() or self.LLOQ == nominal or self.ULOQ == nominal:
                        test_threshold = threshold + Decimal(5)
                    else:
                        test_threshold = threshold

                    re_red = {"batch": batch, "sample": sample, "conc": conc, "threshold": test_threshold}
                    re_yellow = {"batch": batch, "sample": sample, "conc": conc,
                                 "threshold": test_threshold - Decimal(5)}
                    if abs(calc_re) > test_threshold:
                        message = data_flag.get_message("re", **re_red)
                        result.append(red(message, None, self.table_title, self.table_type))
                        failed_sample += 1
                        red_columns[sample] += 1

                    elif (test_threshold - Decimal(5)) <= abs(calc_re) <= test_threshold:
                        message = data_flag.get_message("re", **re_yellow)
                        result.append(yellow(message, None, self.table_title, self.table_type))
                        yellow_columns[sample] += 1

            if (failed_sample / total_sample) > 0.5:
                failed_column += 1

            if (total_sample - failed_sample) < self.required_replicates:
                if self.table_subtype.lower() != "inter":
                    message = data_flag.get_message("failed_batch", **failed_batch_para)
                    result.append(red(message, None, self.table_title, self.table_type))
                red_count += 1

            if total_column < self.required_levels:
                if "lm" in self.analysis_type:
                    message = data_flag.get_message("lm_required_level", **levels_para)
                else:
                    message = data_flag.get_message("sm_required_level", **levels_para)
                result.append(red(message, None, self.table_title, self.table_type))
                self.failed_batch_numbers.append(previous_batch)
                red_count += 1

            if (total_column - failed_column) < self.required_levels:
                message = data_flag.get_message("required_pass_level", **levels_para)
                result.append(red(message, None, self.table_title, self.table_type))
                self.failed_batch_numbers.append(previous_batch)
                red_count += 1

            for header_name in header_names:
                qc_fail_para = {"sample": header_name}
                if yellow_columns[header_name] > 1:
                    message = data_flag.get_message("qc_near_fail", **qc_fail_para)
                    result.append(yellow(message, None, self.table_title, self.table_type))

                if red_columns[header_name] > 1:
                    message = data_flag.get_message("qc_fail", **qc_fail_para)
                    result.append(yellow(message, None, self.table_title, self.table_type))
            self.intra_batch_info[previous_batch] = red_count

        static_df = self.static_df
        if not static_df.empty:
            red_count = 0
            samples = static_df["sample"]
            reported_mean = static_df["mean"]
            reported_sd = static_df["sd"]
            reported_cv = static_df["cv"]
            reported_re = static_df["re"]
            reported_n = static_df["n"]

            calculated_mean = static_df["calc_mean"]
            calculated_sd = static_df["calc_sd"]
            calculated_cv = static_df["calc_cv"]
            calculated_re = static_df["calc_re"]
            calculated_n = static_df["calc_n"]

            per_diff_mean = static_df["per_diff_mean"]
            per_diff_sd = static_df["per_diff_sd"]
            per_diff_cv = static_df["per_diff_cv"]
            per_diff_re = static_df["per_diff_re"]
            per_diff_n = static_df["per_diff_n"]
            nominal_concentrations = static_df["nominal"]

            batches = static_df["run_id"]
            previous_batch = batches[0]

            calculated_tes = static_df["calc_te"]

            if self.found_te:
                per_diff_te = static_df["per_diff_te"]
                reported_te = static_df["te"]

            self.found_te = check_te_column(static_df, previous_batch)
            for index, sample in enumerate(samples):
                te_threshold = self.te_threshold
                batch = batches[index]
                overall_reported_mean = reported_mean[index]
                overall_reported_sd = reported_sd[index]
                overall_reported_cv = reported_cv[index]
                overall_reported_re = reported_re[index]
                overall_reported_n = reported_n[index]

                overall_clc_mean = calculated_mean[index]
                overall_clc_sd = calculated_sd[index]
                overall_clc_cv = calculated_cv[index]
                overall_clc_re = calculated_re[index]
                overall_clc_n = calculated_n[index]

                overall_per_diff_mean = per_diff_mean[index]
                overall_per_diff_sd = per_diff_sd[index]
                overall_per_diff_cv = per_diff_cv[index]
                overall_per_diff_re = per_diff_re[index]
                overall_per_diff_n = per_diff_n[index]

                nominal = nominal_concentrations[index]

                calculated_te = calculated_tes[index]

                if previous_batch != batch:
                    self.intra_batch_info[previous_batch] += red_count
                    red_count = 0
                    previous_batch = batch
                    self.found_te = check_te_column(static_df, previous_batch)

                if "lloq" in str(sample).lower() or self.LLOQ == nominal:
                    lloq_intra_cvs.append(overall_clc_cv)
                    lloq_intra_biases.append(overall_clc_re)

                    test_threshold = threshold + Decimal(5)
                    te_threshold += Decimal(10)

                elif "uloq" in str(sample).lower() or self.ULOQ == nominal:
                    uloq_intra_cvs.append(overall_clc_cv)
                    uloq_intra_biases.append(overall_clc_re)

                    test_threshold = threshold + Decimal(5)
                    te_threshold += Decimal(10)

                else:
                    lmh_intra_cvs.append(overall_clc_cv)
                    lmh_intra_biases.append(overall_clc_re)

                    test_threshold = threshold

                if self.found_te:
                    overall_per_diff_te = per_diff_te[index]
                    overall_reported_te = reported_te[index]

                    te_error_para = {"batch": batch, "sample": sample, "reported_value": overall_reported_te,
                                     "calc_value": utils.format_value(calculated_te)}

                    if str(overall_per_diff_te) != "":
                        if overall_per_diff_te > valid_te_diff:
                            message = batch_stats_flag.get_message("te_error", **te_error_para)
                            result.append(red(message, None, self.table_title, self.table_type))

                else:
                    pass

                te_red = {"batch": batch, "sample": sample, "te_threshold": te_threshold}
                te_yellow = {"batch": batch, "sample": sample, "te_threshold": te_threshold - Decimal(5)}

                if abs(calculated_te) > te_threshold:
                    message = batch_stats_flag.get_message("te", **te_red)
                    result.append(red(message, None, self.table_title, self.table_type))

                elif (te_threshold - Decimal(5)) <= abs(calculated_te) <= te_threshold:
                    message = batch_stats_flag.get_message("te", **te_yellow)
                    result.append(yellow(message, None, self.table_title, self.table_type))

                count_error_para = {"batch": batch, "sample": sample, "reported_value": overall_reported_n,
                                    "calc_value": overall_clc_n}
                mean_para = {"batch": batch, "sample": sample, "reported_value": overall_reported_mean,
                             "calc_value": utils.format_value(overall_clc_mean)}
                cv_para = {"batch": batch, "sample": sample, "reported_value": overall_reported_cv,
                           "calc_value": utils.format_value(overall_clc_cv)}
                re_para = {"batch": batch, "sample": sample, "reported_value": overall_reported_re,
                           "calc_value": utils.format_value(overall_clc_re)}
                sd_para = {"batch": batch, "sample": sample, "reported_value": overall_reported_sd,
                           "calc_value": utils.format_value(overall_clc_sd)}
                replicate_para = {"batch": batch, "sample": sample, "replicate": self.required_replicates}

                if overall_per_diff_n > valid_count_diff:
                    message = batch_stats_flag.get_message("count_error", **count_error_para)
                    result.append(red(message, None, self.table_title, self.table_type))

                if overall_clc_n < Decimal(self.required_replicates):
                    message = batch_stats_flag.get_message("replicate", **replicate_para)
                    result.append(red(message, None, self.table_title, self.table_type))
                    self.failed_batch_numbers.append(batch)
                    red_count += 1

                if overall_per_diff_mean > valid_mean_diff:
                    if abs(overall_clc_mean) <= Decimal(1) and abs(overall_reported_mean) <= Decimal(1):
                        message = batch_stats_flag.get_message("mean_rounding", **mean_para)
                        result.append(yellow(message, None, self.table_title, self.table_type))

                    else:
                        message = batch_stats_flag.get_message("mean_error", **mean_para)
                        result.append(red(message, None, self.table_title, self.table_type))

                if overall_per_diff_cv > valid_re_cv_diff:
                    if abs(overall_clc_cv) <= Decimal(1) and abs(overall_reported_cv) <= Decimal(1):
                        message = batch_stats_flag.get_message("cv_rounding", **cv_para)
                        result.append(yellow(message, None, self.table_title, self.table_type))

                    else:
                        # red_count += 1
                        message = batch_stats_flag.get_message("cv_error", **cv_para)
                        result.append(red(message, None, self.table_title, self.table_type))

                if overall_per_diff_re > valid_re_cv_diff:
                    if abs(overall_clc_re) <= Decimal(1) and abs(overall_reported_re) <= Decimal(1):
                        message = batch_stats_flag.get_message("re_rounding", **re_para)
                        result.append(yellow(message, None, self.table_title, self.table_type))

                    else:
                        # red_count += 1
                        message = batch_stats_flag.get_message("re_error", **re_para)
                        result.append(red(message, None, self.table_title, self.table_type))

                if overall_per_diff_sd > valid_re_cv_diff:
                    if abs(overall_clc_sd) <= Decimal(1) and abs(overall_reported_sd) <= Decimal(1):
                        message = batch_stats_flag.get_message("sd_rounding", **sd_para)
                        result.append(yellow(message, None, self.table_title, self.table_type))

                    else:
                        message = batch_stats_flag.get_message("sd_error", **sd_para)
                        result.append(red(message, None, self.table_title, self.table_type))

                threshold_red_para = {"batch": batch, "sample": sample, "threshold": test_threshold}
                threshold_yellow_para = {"batch": batch, "sample": sample, "threshold": test_threshold - Decimal(5)}
                failed_batch_para = {"batch": batch, "col_level": self.get_column_level(sample)}

                if abs(overall_clc_re) > test_threshold:
                    message = batch_stats_flag.get_message("re", **threshold_red_para)
                    result.append(red(message, None, self.table_title, self.table_type))
                    batch_failed_message = batch_stats_flag.get_message("failed_batch", **failed_batch_para)
                    result.append(red(batch_failed_message, None, self.table_title, self.table_type))
                    self.failed_batch_numbers.append(batch)
                    red_count += 1

                elif (test_threshold - Decimal(5)) <= abs(overall_clc_re) <= test_threshold:
                    message = batch_stats_flag.get_message("re", **threshold_yellow_para)
                    result.append(yellow(message, None, self.table_title, self.table_type))

                if str(overall_clc_cv) != "":
                    if abs(overall_clc_cv) > test_threshold:
                        message = batch_stats_flag.get_message("cv", **threshold_red_para)
                        result.append(red(message, None, self.table_title, self.table_type))
                        batch_failed_message = batch_stats_flag.get_message("failed_batch", **failed_batch_para)
                        result.append(red(batch_failed_message, None, self.table_title, self.table_type))
                        self.failed_batch_numbers.append(batch)
                        red_count += 1

                    elif (test_threshold - Decimal(5)) <= abs(overall_clc_cv) <= test_threshold:
                        message = batch_stats_flag.get_message("cv", **threshold_yellow_para)
                        result.append(yellow(message, None, self.table_title, self.table_type))

            self.intra_batch_info[previous_batch] += red_count
            self.failed_batch_numbers = list(set(self.failed_batch_numbers))
            self.passed_batch_numbers = [x for x in self.total_batches if x not in self.failed_batch_numbers]

            if len(self.passed_batch_numbers) < self.required_batches:
                message = batch_stats_flag.get_message("required_batches", required_batches=self.required_batches)
                result.append(red(message, None, self.table_title, self.table_type))

            message = batch_stats_flag.get_message("no_pooled_verification")
            result.append(red(message, None, self.table_title, self.table_type))
            message = batch_stats_flag.get_message("no_ancillary_verification")
            result.append(red(message, None, self.table_title, self.table_type))

        final_static_df = self.final_static_df
        inter_red_count = 0
        if not final_static_df.empty:
            samples = final_static_df["sample"]
            reported_mean = final_static_df["mean"]
            reported_cv = final_static_df["cv"]
            reported_re = final_static_df["re"]
            reported_n = final_static_df["n"]

            calculated_mean = final_static_df["calc_mean"]
            calculated_cv = final_static_df["calc_cv"]
            calculated_re = final_static_df["calc_re"]
            calculated_n = final_static_df["calc_n"]

            per_diff_mean = final_static_df["per_diff_mean"]
            per_diff_cv = final_static_df["per_diff_cv"]
            per_diff_re = final_static_df["per_diff_re"]
            per_diff_n = final_static_df["per_diff_n"]
            nominal_concentrations = final_static_df["nominal"]

            if self.found_sd:
                reported_sd = final_static_df["sd"]
                calculated_sd = final_static_df["calc_sd"]
                per_diff_sd = final_static_df["per_diff_sd"]

            calculated_tes = final_static_df["calc_te"]
            if self.found_final_te:
                per_diff_te = final_static_df["per_diff_te"]
                reported_te = final_static_df["te"]

            for index, sample in enumerate(samples):
                te_threshold = self.te_threshold
                overall_reported_mean = reported_mean[index]
                overall_reported_cv = reported_cv[index]
                overall_reported_re = reported_re[index]
                overall_reported_n = reported_n[index]

                overall_clc_mean = calculated_mean[index]
                overall_clc_cv = calculated_cv[index]
                overall_clc_re = calculated_re[index]
                overall_clc_n = calculated_n[index]

                overall_per_diff_mean = per_diff_mean[index]
                overall_per_diff_cv = per_diff_cv[index]
                overall_per_diff_re = per_diff_re[index]
                overall_per_diff_n = per_diff_n[index]
                nominal = nominal_concentrations[index]

                calculated_te = calculated_tes[index]
                if self.found_final_te:
                    overall_per_diff_te = per_diff_te[index]
                    overall_reported_te = reported_te[index]

                if "lloq" in str(sample).lower() or self.LLOQ == nominal:
                    lloq_inter_cv = overall_clc_cv
                    lloq_inter_bias = overall_clc_re
                    if self.found_final_te:
                        lloq_inter_te = calculated_te

                    test_threshold = threshold + Decimal(5)
                    te_threshold += Decimal(10)

                elif "uloq" in str(sample).lower() or self.ULOQ == nominal:
                    uloq_inter_cv = overall_clc_cv
                    uloq_inter_bias = overall_clc_re
                    if self.found_final_te:
                        uloq_inter_te = calculated_te

                    test_threshold = threshold + Decimal(5)
                    te_threshold += Decimal(10)

                else:
                    lmh_inter_cvs.append(overall_clc_cv)
                    lmh_inter_biases.append(overall_clc_re)
                    if self.found_final_te:
                        lmh_inter_tes.append(calculated_te)
                    test_threshold = threshold

                count_error_para = {"sample": sample, "reported_value": overall_reported_n, "calc_value": overall_clc_n}

                mean_para = {"sample": sample, "reported_value": overall_reported_mean,
                             "calc_value": utils.format_value(overall_clc_mean)}

                cv_para = {"sample": sample, "reported_value": overall_reported_cv,
                           "calc_value": utils.format_value(overall_clc_cv)}

                re_para = {"sample": sample, "reported_value": overall_reported_re,
                           "calc_value": utils.format_value(overall_clc_re)}

                threshold_red_para = {"sample": sample, "threshold": test_threshold}
                threshold_yellow_para = {"sample": sample, "threshold": test_threshold - Decimal(5)}

                if self.found_final_te:
                    te_error_para = {"sample": sample, "reported_value": overall_reported_te,
                                     "calc_value": utils.format_value(calculated_te)}

                    if str(overall_per_diff_te) != "":
                        if overall_per_diff_te > valid_te_diff:
                            message = stats_flag.get_message("te_error", **te_error_para)
                            result.append(red(message, None, self.table_title, self.table_type))

                else:
                    pass

                if abs(calculated_te) > te_threshold:
                    message = stats_flag.get_message("te", sample=sample, te_threshold=te_threshold)
                    result.append(red(message, None, self.table_title, self.table_type))

                elif (te_threshold - Decimal(5)) <= abs(calculated_te) <= te_threshold:
                    message = stats_flag.get_message("te", sample=sample, te_threshold=te_threshold - Decimal(5))
                    result.append(yellow(message, None, self.table_title, self.table_type))

                if self.found_sd:
                    overall_reported_sd = reported_sd[index]
                    overall_clc_sd = calculated_sd[index]
                    overall_per_diff_sd = per_diff_sd[index]

                    sd_para = {"sample": sample, "reported_value": overall_reported_sd,
                               "calc_value": utils.format_value(overall_clc_sd)}
                    if overall_per_diff_sd > valid_re_cv_diff:
                        message = stats_flag.get_message("sd_rounding", **sd_para)
                        if abs(overall_clc_sd) <= Decimal(1) and abs(overall_reported_sd) <= Decimal(1):
                            result.append(yellow(message, None, self.table_title, self.table_type))

                        else:
                            message = stats_flag.get_message("sd_error", **sd_para)
                            result.append(red(message, None, self.table_title, self.table_type))

                if overall_per_diff_n > valid_count_diff:
                    message = stats_flag.get_message("count_error", **count_error_para)
                    result.append(red(message, None, self.table_title, self.table_type))

                if overall_per_diff_mean > valid_mean_diff:
                    if abs(overall_clc_mean) <= Decimal(1) and abs(overall_reported_mean) <= Decimal(1):
                        message = stats_flag.get_message("mean_rounding", **mean_para)
                        result.append(yellow(message, None, self.table_title, self.table_type))

                    else:
                        message = stats_flag.get_message("mean_error", **mean_para)
                        result.append(red(message, None, self.table_title, self.table_type))

                if overall_per_diff_cv > valid_re_cv_diff:
                    if abs(overall_clc_cv) <= Decimal(1) and abs(overall_reported_cv) <= Decimal(1):
                        message = stats_flag.get_message("cv_rounding", **cv_para)
                        result.append(yellow(message, None, self.table_title, self.table_type))

                    else:
                        # inter_red_count += 1
                        message = stats_flag.get_message("cv_error", **cv_para)
                        result.append(red(message, None, self.table_title, self.table_type))

                if overall_per_diff_re > valid_re_cv_diff:
                    if abs(overall_clc_re) <= Decimal(1) and abs(overall_reported_re) <= Decimal(1):
                        message = stats_flag.get_message("re_rounding", **re_para)
                        result.append(yellow(message, None, self.table_title, self.table_type))

                    else:
                        # inter_red_count += 1
                        message = stats_flag.get_message("re_error", **re_para)
                        result.append(red(message, None, self.table_title, self.table_type))
                failed_batch_message = stats_flag.get_message("failed_batch", sample=sample)

                if abs(overall_clc_re) > test_threshold:
                    message = stats_flag.get_message("re", **threshold_red_para)
                    result.append(red(message, None, self.table_title, self.table_type))
                    result.append(red(failed_batch_message, None, self.table_title, self.table_type))
                    inter_red_count += 1
                elif (test_threshold - Decimal(5)) <= abs(overall_clc_re) <= test_threshold:
                    message = stats_flag.get_message("re", **threshold_yellow_para)
                    result.append(yellow(message, None, self.table_title, self.table_type))

                if abs(overall_clc_cv) > test_threshold:
                    message = stats_flag.get_message("cv", **threshold_red_para)
                    result.append(red(message, None, self.table_title, self.table_type))
                    result.append(red(failed_batch_message, None, self.table_title, self.table_type))
                    inter_red_count += 1

                elif (test_threshold - Decimal(5)) <= abs(overall_clc_cv) <= test_threshold:
                    message = stats_flag.get_message("cv", **threshold_yellow_para)
                    result.append(yellow(message, None, self.table_title, self.table_type))

            message = stats_flag.get_message("no_inter_te_verification")
            result.append(red(message, None, self.table_title, self.table_type))

        if self.table_subtype != "inter":
            for batch, val in self.intra_batch_info.items():
                if val == 0:
                    message = batch_stats_flag.get_message("ap_demonstrated", batch=batch)
                    self.result.append(green(message, None, self.table_title, self.table_type))

        if self.table_subtype != "intra" and self.table_subtype != "inter":
            if sum(self.intra_batch_info.values()) == 0 and inter_red_count == 0:
                message = stats_flag.get_message('ap_demonstrated', batch_count=len(self.intra_batch_info))
                self.result.append(green(message, None, self.table_title, self.table_type))

        if self.table_subtype == "intra":
            message = batch_stats_flag.get_message("no_statistics")
            self.result.append(yellow(message, None, self.table_title, self.table_type))

        elif self.table_subtype == "inter":
            for batch in set(batches):
                message = stats_flag.get_message("no_statistics", batch=batch)
                self.result.append(yellow(message, None, self.table_title, self.table_type))
        else:
            pass
        if len(result) == 0:
            result.append(green('All values within nominal range', None, self.table_title, self.table_type))

            # intra accuracy and precision
        if len(lloq_intra_cvs) > 0:
            lloq_intra_cvs = [x for x in lloq_intra_cvs if str(x).strip() != ""]
            lloq_intra_high_precision = round(max(lloq_intra_cvs), 2)
            lloq_intra_low_precision = round(min(lloq_intra_cvs), 2)
            self.intra_precisions["lloq_precision"] = [lloq_intra_low_precision, lloq_intra_high_precision]

        if len(uloq_intra_cvs) > 0:
            uloq_intra_cvs = [x for x in uloq_intra_cvs if str(x).strip() != ""]
            uloq_intra_high_precision = round(max(uloq_intra_cvs), 2)
            uloq_intra_low_precision = round(min(uloq_intra_cvs), 2)
            self.intra_precisions["uloq_precision"] = [uloq_intra_low_precision, uloq_intra_high_precision]

        if len(lmh_intra_cvs) > 0:
            lmh_intra_cvs = [x for x in lmh_intra_cvs if str(x).strip() != ""]
            lmh_intra_high_precision = round(max(lmh_intra_cvs), 2)
            lmh_intra_low_precision = round(min(lmh_intra_cvs), 2)
            self.intra_precisions["lmh_precision"] = [lmh_intra_low_precision, lmh_intra_high_precision]

        if len(lloq_intra_biases) > 0:
            lloq_intra_biases = [x for x in lloq_intra_biases if str(x).strip() != ""]
            lloq_intra_high_accuracy = round(max(lloq_intra_biases), 2)
            lloq_intra_low_accuracy = round(min(lloq_intra_biases), 2)
            self.intra_accuracies["lloq_accuracy"] = [lloq_intra_low_accuracy, lloq_intra_high_accuracy]

        if len(uloq_intra_biases) > 0:
            uloq_intra_biases = [x for x in uloq_intra_biases if str(x).strip() != ""]
            uloq_intra_high_accuracy = round(max(uloq_intra_biases), 2)
            uloq_intra_low_accuracy = round(min(uloq_intra_biases), 2)
            self.intra_accuracies["uloq_accuracy"] = [uloq_intra_low_accuracy, uloq_intra_high_accuracy]

        if len(lmh_intra_biases) > 0:
            lmh_intra_biases = [x for x in lmh_intra_biases if str(x).strip() != ""]
            lmh_intra_high_accuracy = round(max(lmh_intra_biases), 2)
            lmh_intra_low_accuracy = round(min(lmh_intra_biases), 2)
            self.intra_accuracies["lmh_accuracy"] = [lmh_intra_low_accuracy, lmh_intra_high_accuracy]

        # inter accuracy and precision
        if lloq_inter_cv is not None and str(lloq_inter_cv).strip() != "":
            self.inter_precisions["lloq_precision"] = round(lloq_inter_cv, 2)

        if uloq_inter_cv is not None and str(uloq_inter_cv).strip() != "":
            self.inter_precisions["uloq_precision"] = round(uloq_inter_cv, 2)

        if len(lmh_inter_cvs) > 0:
            lmh_inter_cvs = [x for x in lmh_inter_cvs if str(x).strip() != ""]
            lmh_inter_high_precision = round(max(lmh_inter_cvs), 2)
            lmh_inter_low_precision = round(min(lmh_inter_cvs), 2)
            self.inter_precisions["lmh_precision"] = [lmh_inter_low_precision, lmh_inter_high_precision]

        if lloq_inter_bias is not None and str(lloq_inter_bias).strip() != "":
            self.inter_accuracies["lloq_accuracy"] = round(lloq_inter_bias, 2)

        if uloq_inter_bias is not None and str(uloq_inter_bias).strip() != "":
            self.inter_accuracies["uloq_accuracy"] = round(uloq_inter_bias, 2)

        if lloq_inter_te is not None and str(lloq_inter_te).strip() != "":
            self.inter_tes["lloq_te"] = round(lloq_inter_te, 2)

        if uloq_inter_te is not None and str(uloq_inter_te).strip() != "":
            self.inter_tes["uloq_te"] = round(uloq_inter_te, 2)

        if len(lmh_inter_tes) > 0:
            lmh_inter_tes = [x for x in lmh_inter_tes if str(x).strip() != ""]
            lmh_inter_high_te = round(max(lmh_inter_tes), 2)
            lmh_inter_low_te = round(min(lmh_inter_tes), 2)
            self.inter_tes["lmh_te"] = [lmh_inter_low_te, lmh_inter_high_te]

        if len(lmh_inter_biases) > 0:
            lmh_inter_biases = [x for x in lmh_inter_biases if str(x).strip() != ""]
            lmh_inter_high_accuracy = round(max(lmh_inter_biases), 2)
            lmh_inter_low_accuracy = round(min(lmh_inter_biases), 2)
            self.inter_accuracies["lmh_accuracy"] = [lmh_inter_low_accuracy, lmh_inter_high_accuracy]

        self.result += result

    @staticmethod
    def get_column_level(column):
        val = re.search(r"^qc\s*[0-9]*", column, re.I)
        if val:
            return val.group(0)
        else:
            return ""

    @staticmethod
    def get_error_dict():
        data_flag = FlagProperties("ppd_ap_data")
        stats_flag = FlagProperties("ppd_ap_batch_stats")
        final_stats_flag = FlagProperties("ppd_ap_stats")
        return data_flag, stats_flag, final_stats_flag


class AccuracyAndPrecision(Table, Flags):
    def __init__(self, parsed_table, analytes, template_type, LLOQ, ULOQ):
        Table.__init__(self, parsed_table, template_type)
        self.total_batches = []
        self.orientation = ""
        self.conc_unit = ""
        self.found_final_te = False
        self.LLOQ = LLOQ
        self.ULOQ = ULOQ
        self.te_threshold = Decimal(self.threshold_values.get_message("te_threshold"))
        self.analyte = get_analyte(analytes, self.tb_title)
        self.header_values = []
        self.failed_batch_numbers = []
        self.passed_batch_numbers = []
        self.intra_precisions = {}
        self.intra_accuracies = {}
        self.inter_precisions = {}
        self.inter_accuracies = {}
        self.inter_tes = {}
        self.intra_batch_info = {}
        self.inter_batch_info = {}
        if "sm" in self.analysis_type:
            self.required_levels = 4
            self.required_replicates = 5
            self.required_batches = 3
        elif "lm" in self.analysis_type:
            self.required_levels = 5
            self.required_replicates = 3
            self.required_batches = 6
        self.found_duplicate_run = False

    def validate_table_format(self):
        error_messages = self.error_messages
        required_col = {"Run ID"}
        missing_col = required_col.intersection(self.missing_col)
        if missing_col:
            message = error_messages.get_message("missing_col", col_names=", ".join(missing_col))
            self.result.append(red(message, None, self.table_title, self.table_type))
            self.valid_format = False

        conc_units = utils.find_units(self.table_df.columns)
        if conc_units == "":
            message = error_messages.get_message("missing_conc_unit")
            self.result.append(red(message, None, self.table_title, self.table_type))
            self.valid_format = False
        else:
            self.conc_unit = conc_units

    def process_table(self):
        table = self.table_df
        self.orientation = utils.check_orientation(table)
        final_df = pd.DataFrame()
        calc_static_run_df = pd.DataFrame()
        final_static_run_df = pd.DataFrame()
        static_df = pd.DataFrame()
        self.found_re = utils.find_re_column(table)
        final_static_df = pd.DataFrame()
        try:
            if self.orientation == "HZ":
                batch_tables, duplicate_id = utils.split_run_df(table)
                if len(duplicate_id) > 0:
                    self.result.append(
                        red(f"Unable to process table due to replicated batch run ID ({' '.join(duplicate_id)}) "
                            f"found for more than one set of results", None, self.table_title, self.table_type))
                    self.found_duplicate_run = True
                else:
                    for i in range(len(batch_tables)):
                        table_df = batch_tables[i]["table"]
                        data, static = utils.split_static_table(table_df)
                        static = utils.drop_blank_col(static)
                        try:
                            if i == len(batch_tables) - 1:
                                static, final_static_df = utils.split_static_table(static, "two")
                        except Exception as e:
                            pass
                        if self.table_subtype == "inter" and not static.empty:
                            final_static_df = static
                            static = pd.DataFrame()
                        table_df, calc_static, found_re, found_cycle = utils.format_df_calc_stats(data)
                        final_df = pd.concat([final_df, table_df]).reset_index(drop=True)

                        table_df = table_df.reset_index(drop=True)

                        if not static.empty:
                            static_df = utils.process_static_df(static, self.table_type)
                            static_df.insert(loc=0, column="run_id", value=table_df["run_id"][0])
                            calc_static_run_df = pd.concat([calc_static_run_df, calc_static]).reset_index(drop=True)
                            final_static_run_df = pd.concat([final_static_run_df, static_df]).reset_index(drop=True)

                    final_df["calc_re"] = utils.calculate_re(final_df["conc"], final_df["nominal"])
                    if self.found_re:
                        final_df["per_diff_re"] = utils.calculate_per_diff(final_df["re"], final_df["calc_re"])

                    if not final_static_run_df.empty:
                        static_df, self.found_sd = utils.concat_static_df(final_static_run_df, calc_static_run_df)

                    self.found_te = utils.find_te_column(static_df)
                    if self.found_te:
                        static_df["te"] = static_df["te"].apply(utils.parse_decimal)
                        static_df["calc_te"] = abs(static_df["calc_re"]) + abs(static_df["calc_cv"])
                        static_df["per_diff_te"] = abs((static_df["calc_te"] - static_df["te"]) / static_df["calc_te"])

                    if not final_static_df.empty:
                        final_static_df = utils.process_static_df(final_static_df, self.table_type)
                        final_calc_static_df = pd.DataFrame()
                        samples = final_df["column"].unique()
                        for sample in samples:
                            values = final_df[final_df["column"] == sample]["conc"].to_list()
                            nominal = final_df[final_df["column"] == sample]["nominal"].to_list()[0]
                            final_calc_static_df = pd.concat(
                                [final_calc_static_df, utils.build_static_df(values, nominal)]).reset_index(drop=True)
                        final_static_df, self.found_sd = utils.concat_static_df(final_static_df, final_calc_static_df)
                        self.found_final_te = utils.find_te_column(final_static_df)
                        if self.found_final_te:
                            final_static_df["te"] = final_static_df["te"].apply(utils.parse_decimal)
                            final_static_df["calc_te"] = abs(final_static_df["calc_re"]) + abs(
                                final_static_df["calc_cv"])
                            final_static_df["per_diff_te"] = utils.calculate_per_diff(final_static_df["te"],
                                                                                      final_static_df["calc_te"])
            final_df = final_df.fillna("")
            final_df = utils.fill_val(final_df)
            self.data_df = final_df.fillna("")
            self.static_df = static_df.fillna("")
            self.final_static_df = final_static_df.fillna("")
        except Exception as e:
            message = self.error_messages.get_message("data_error")
            self.result.append(red(message, None, self.table_title, self.table_type))
            self.valid_data = False

    def validate_table_data(self):
        data_flag, batch_stats_flag, stats_flag = self.get_error_dict()
        result = []

        valid_count_diff = Decimal(self.valid_difference_values.get_message("difference_n"))
        valid_mean_diff = Decimal(self.valid_difference_values.get_message("difference_mean"))
        valid_re_cv_diff = Decimal(self.valid_difference_values.get_message("difference_re_cv"))
        valid_te_diff = Decimal(self.valid_difference_values.get_message("difference_te"))

        lloq_intra_cvs = []
        uloq_intra_cvs = []
        lmh_intra_cvs = []

        lloq_intra_biases = []
        uloq_intra_biases = []
        lmh_intra_biases = []

        lmh_inter_cvs = []
        uloq_inter_cv = None
        lloq_inter_cv = None

        lloq_inter_bias = None
        uloq_inter_bias = None
        lmh_inter_biases = []

        lmh_inter_tes = []
        lloq_inter_te = None
        uloq_inter_te = None

        threshold = Decimal(self.threshold_values.get_message("re_cv_threshold"))
        if "lm" in self.analysis_type:
            threshold += Decimal(5)

        table = self.data_df
        if not table.empty:
            red_count = 0
            self.header_values = table["nominal"].unique()
            batches = table["run_id"].to_list()
            all_columns = table["column"]
            header_names = list(table["column"].unique())
            all_conc = table["conc"]
            all_calc_res = table["calc_re"]
            nominal_concentrations = table["nominal"]
            if self.found_re:
                reported_res = table["re"]
                per_diff_res = table["per_diff_re"]

            self.total_batches = list(set(batches))
            total_column = len(header_names)
            dict_val = [0]*len(header_names)
            yellow_columns = dict(zip(header_names, dict_val))
            red_columns = dict(zip(header_names, dict_val))
            previous_batch = batches[0]
            previous_column = all_columns[0]
            failed_sample = 0
            total_sample = 0
            failed_column = 0

            failed_batch_para = {"batch": previous_batch, "col_level": self.get_column_level(previous_column)}
            levels_para = {"batch": previous_batch, "replicate": self.required_levels}

            for index, column in enumerate(all_columns):
                if self.orientation == "HZ":
                    batch = batches[index]
                    conc = all_conc[index]
                    calc_re = all_calc_res[index]
                    nominal = nominal_concentrations[index]

                    failed_batch_para = {"batch": previous_batch, "col_level": self.get_column_level(previous_column)}
                    levels_para = {"batch": previous_batch, "required_levels": self.required_levels}

                    if previous_column != column:
                        col_values = table[(table["run_id"] == previous_batch) & (table["column"] == previous_column)][
                            "conc"].replace("", np.nan).dropna()
                        if (failed_sample / total_sample) > 0.5:
                            failed_column += 1
                        if (total_sample - failed_sample) < self.required_replicates:
                            if self.table_subtype.lower() != "inter":
                                message = data_flag.get_message("failed_batch", **failed_batch_para)
                                result.append(red(message, None, self.table_title, self.table_type))
                            red_count += 1
                        elif col_values.empty:
                            red_count += 1

                        previous_column = column
                        failed_sample = 0
                        total_sample = 0

                    if previous_batch != batch:
                        if total_column < self.required_levels:
                            if "lm" in self.analysis_type:
                                message = data_flag.get_message("lm_required_level", **levels_para)
                            else:
                                message = data_flag.get_message("sm_required_level", **levels_para)
                            result.append(red(message, None, self.table_title, self.table_type))
                            self.failed_batch_numbers.append(previous_batch)
                            red_count += 1

                        if (total_column - failed_column) < self.required_levels:
                            message = data_flag.get_message("required_pass_level", **levels_para)
                            result.append(red(message, None, self.table_title, self.table_type))
                            self.failed_batch_numbers.append(previous_batch)
                            red_count += 1

                        self.intra_batch_info[previous_batch] = red_count
                        red_count = 0
                        previous_batch = batch
                        failed_column = 0

                    total_sample += 1

                    if self.found_re:
                        reported_re = reported_res[index]
                        per_diff_re = per_diff_res[index]
                        re_error_para = {"batch": batch, "column": column, "reported_value": reported_re,
                                         "calc_value": utils.format_value(calc_re), "conc": conc}

                        if per_diff_re > valid_re_cv_diff:
                            if abs(calc_re) <= Decimal(1) and abs(reported_re) <= Decimal(1):
                                message = data_flag.get_message("re_rounding", **re_error_para)
                                result.append(yellow(message, None, self.table_title, self.table_type))
                            else:
                                message = data_flag.get_message("re_error", **re_error_para)
                                result.append(red(message, None, self.table_title, self.table_type))

                    if str(calc_re) != "" and calc_re is not None:
                        if "lloq" in str(column).lower() or "uloq" in str(
                                column).lower() or self.LLOQ == nominal or self.ULOQ == nominal:
                            test_threshold = threshold + Decimal(5)
                        else:
                            test_threshold = threshold

                        re_red = {"batch": batch, "column": column, "conc": conc, "threshold": test_threshold}
                        re_yellow = {"batch": batch, "column": column, "conc": conc,
                                     "threshold": test_threshold - Decimal(5)}
                        if abs(calc_re) > test_threshold:
                            message = data_flag.get_message("re", **re_red)
                            result.append(red(message, None, self.table_title, self.table_type))
                            failed_sample += 1
                            red_columns[column] += 1

                        elif (test_threshold - Decimal(5)) <= abs(calc_re) <= test_threshold:
                            message = data_flag.get_message("re", **re_yellow)
                            result.append(yellow(message, None, self.table_title, self.table_type))
                            yellow_columns[column] += 1

            if (failed_sample / total_sample) > 0.5:
                failed_column += 1

            if (total_sample - failed_sample) < self.required_replicates:
                if self.table_subtype.lower() != "inter":
                    message = data_flag.get_message("failed_batch", **failed_batch_para)
                    result.append(red(message, None, self.table_title, self.table_type))
                red_count += 1

            if total_column < self.required_levels:
                if "lm" in self.analysis_type:
                    message = data_flag.get_message("lm_required_level", **levels_para)
                else:
                    message = data_flag.get_message("sm_required_level", **levels_para)
                result.append(red(message, None, self.table_title, self.table_type))
                self.failed_batch_numbers.append(previous_batch)
                red_count += 1

            if (total_column - failed_column) < self.required_levels:
                message = data_flag.get_message("required_pass_level", **levels_para)
                result.append(red(message, None, self.table_title, self.table_type))
                self.failed_batch_numbers.append(previous_batch)
                red_count += 1

            for header_name in header_names:
                qc_fail_para = {"column": header_name}
                if yellow_columns[header_name] > 1:
                    message = data_flag.get_message("qc_near_fail", **qc_fail_para)
                    result.append(yellow(message, None, self.table_title, self.table_type))

                if red_columns[header_name] > 1:
                    message = data_flag.get_message("qc_fail", **qc_fail_para)
                    result.append(yellow(message, None, self.table_title, self.table_type))
            self.intra_batch_info[previous_batch] = red_count

        static_df = self.static_df
        if not static_df.empty:
            red_count = 0
            all_column = static_df["column"]
            reported_mean = static_df["mean"]
            reported_sd = static_df["sd"]
            reported_cv = static_df["cv"]
            reported_re = static_df["re"]
            reported_n = static_df["n"]

            calculated_mean = static_df["calc_mean"]
            calculated_sd = static_df["calc_sd"]
            calculated_cv = static_df["calc_cv"]
            calculated_re = static_df["calc_re"]
            calculated_n = static_df["calc_n"]

            per_diff_mean = static_df["per_diff_mean"]
            per_diff_sd = static_df["per_diff_sd"]
            per_diff_cv = static_df["per_diff_cv"]
            per_diff_re = static_df["per_diff_re"]
            per_diff_n = static_df["per_diff_n"]

            nominal_concentrations = static_df["nominal"]

            batches = static_df["run_id"]
            previous_batch = batches[0]

            if self.found_te:
                calculated_tes = static_df["calc_te"]
                per_diff_te = static_df["per_diff_te"]
                reported_te = static_df["te"]
            self.found_te = check_te_column(static_df, previous_batch)
            for index, column in enumerate(all_column):
                te_threshold = self.te_threshold
                batch = batches[index]
                overall_reported_mean = reported_mean[index]
                overall_reported_sd = reported_sd[index]
                overall_reported_cv = reported_cv[index]
                overall_reported_re = reported_re[index]
                overall_reported_n = reported_n[index]

                overall_clc_mean = calculated_mean[index]
                overall_clc_sd = calculated_sd[index]
                overall_clc_cv = calculated_cv[index]
                overall_clc_re = calculated_re[index]
                overall_clc_n = calculated_n[index]

                overall_per_diff_mean = per_diff_mean[index]
                overall_per_diff_sd = per_diff_sd[index]
                overall_per_diff_cv = per_diff_cv[index]
                overall_per_diff_re = per_diff_re[index]
                overall_per_diff_n = per_diff_n[index]

                nominal = nominal_concentrations[index]

                if previous_batch != batch:
                    self.intra_batch_info[previous_batch] += red_count
                    red_count = 0
                    previous_batch = batch
                    self.found_te = check_te_column(static_df, previous_batch)

                if "lloq" in str(column).lower() or self.LLOQ == nominal:
                    lloq_intra_cvs.append(overall_clc_cv)
                    lloq_intra_biases.append(overall_clc_re)

                    test_threshold = threshold + Decimal(5)
                    te_threshold += Decimal(10)

                elif "uloq" in str(column).lower() or self.ULOQ == nominal:
                    uloq_intra_cvs.append(overall_clc_cv)
                    uloq_intra_biases.append(overall_clc_re)

                    test_threshold = threshold + Decimal(5)
                    te_threshold += Decimal(10)

                else:
                    lmh_intra_cvs.append(overall_clc_cv)
                    lmh_intra_biases.append(overall_clc_re)

                    test_threshold = threshold

                if self.found_te:
                    calculated_te = calculated_tes[index]
                    overall_per_diff_te = per_diff_te[index]
                    overall_reported_te = reported_te[index]

                    te_error_para = {"batch": batch, "column": column, "reported_value": overall_reported_te,
                                     "calc_value": utils.format_value(calculated_te)}
                    te_red = {"batch": batch, "column": column, "te_threshold": te_threshold}
                    te_yellow = {"batch": batch, "column": column, "te_threshold": te_threshold - Decimal(5)}

                    if str(overall_per_diff_te) != "":
                        if overall_per_diff_te > valid_te_diff:
                            message = batch_stats_flag.get_message("te_error", **te_error_para)
                            result.append(red(message, None, self.table_title, self.table_type))

                        if abs(calculated_te) > te_threshold:
                            message = batch_stats_flag.get_message("te", **te_red)
                            result.append(red(message, None, self.table_title, self.table_type))

                        elif (te_threshold - Decimal(5)) <= abs(calculated_te) <= te_threshold:
                            message = batch_stats_flag.get_message("te", **te_yellow)
                            result.append(yellow(message, None, self.table_title, self.table_type))
                else:
                    if "lm" in self.analysis_type:
                        message = batch_stats_flag.get_message("te_n_found", batch=batch)
                        result.append(yellow(message, None, self.table_title, self.table_type))

                if self.orientation == "HZ":

                    count_error_para = {"batch": batch, "column": column, "reported_value": overall_reported_n,
                                        "calc_value": overall_clc_n}
                    mean_para = {"batch": batch, "column": column, "reported_value": overall_reported_mean,
                                 "calc_value": utils.format_value(overall_clc_mean)}
                    cv_para = {"batch": batch, "column": column, "reported_value": overall_reported_cv,
                               "calc_value": utils.format_value(overall_clc_cv)}
                    re_para = {"batch": batch, "column": column, "reported_value": overall_reported_re,
                               "calc_value": utils.format_value(overall_clc_re)}
                    sd_para = {"batch": batch, "column": column, "reported_value": overall_reported_sd,
                               "calc_value": utils.format_value(overall_clc_sd)}
                    replicate_para = {"batch": batch, "column": column, "replicate": self.required_replicates}

                    if overall_per_diff_n > valid_count_diff:
                        message = batch_stats_flag.get_message("count_error", **count_error_para)
                        result.append(red(message, None, self.table_title, self.table_type))

                    if overall_clc_n < Decimal(self.required_replicates):
                        message = batch_stats_flag.get_message("replicate", **replicate_para)
                        result.append(red(message, None, self.table_title, self.table_type))
                        self.failed_batch_numbers.append(batch)
                        red_count += 1

                    if overall_per_diff_mean > valid_mean_diff:
                        if abs(overall_clc_mean) <= Decimal(1) and abs(overall_reported_mean) <= Decimal(1):
                            message = batch_stats_flag.get_message("mean_rounding", **mean_para)
                            result.append(yellow(message, None, self.table_title, self.table_type))

                        else:
                            message = batch_stats_flag.get_message("mean_error", **mean_para)
                            result.append(red(message, None, self.table_title, self.table_type))

                    if overall_per_diff_cv > valid_re_cv_diff:
                        if abs(overall_clc_cv) <= Decimal(1) and abs(overall_reported_cv) <= Decimal(1):
                            message = batch_stats_flag.get_message("cv_rounding", **cv_para)
                            result.append(yellow(message, None, self.table_title, self.table_type))

                        else:
                            # red_count += 1
                            message = batch_stats_flag.get_message("cv_error", **cv_para)
                            result.append(red(message, None, self.table_title, self.table_type))

                    if overall_per_diff_re > valid_re_cv_diff:
                        if abs(overall_clc_re) <= Decimal(1) and abs(overall_reported_re) <= Decimal(1):
                            message = batch_stats_flag.get_message("re_rounding", **re_para)
                            result.append(yellow(message, None, self.table_title, self.table_type))

                        else:
                            # red_count += 1
                            message = batch_stats_flag.get_message("re_error", **re_para)
                            result.append(red(message, None, self.table_title, self.table_type))

                    if overall_per_diff_sd > valid_re_cv_diff:
                        if abs(overall_clc_sd) <= Decimal(1) and abs(overall_reported_sd) <= Decimal(1):
                            message = batch_stats_flag.get_message("sd_rounding", **sd_para)
                            result.append(yellow(message, None, self.table_title, self.table_type))

                        else:
                            message = batch_stats_flag.get_message("sd_error", **sd_para)
                            result.append(red(message, None, self.table_title, self.table_type))

                    threshold_red_para = {"batch": batch, "column": column, "threshold": test_threshold}
                    threshold_yellow_para = {"batch": batch, "column": column, "threshold": test_threshold - Decimal(5)}
                    failed_batch_para = {"batch": batch, "col_level": self.get_column_level(column)}

                    if abs(overall_clc_re) > test_threshold:
                        message = batch_stats_flag.get_message("re", **threshold_red_para)
                        result.append(red(message, None, self.table_title, self.table_type))
                        batch_failed_message = batch_stats_flag.get_message("failed_batch", **failed_batch_para)
                        result.append(red(batch_failed_message, None, self.table_title, self.table_type))
                        self.failed_batch_numbers.append(batch)
                        red_count += 1

                    elif (test_threshold - Decimal(5)) <= abs(overall_clc_re) <= test_threshold:
                        message = batch_stats_flag.get_message("re", **threshold_yellow_para)
                        result.append(yellow(message, None, self.table_title, self.table_type))

                    if str(overall_clc_cv) != "":
                        if abs(overall_clc_cv) > test_threshold:
                            message = batch_stats_flag.get_message("cv", **threshold_red_para)
                            result.append(red(message, None, self.table_title, self.table_type))
                            batch_failed_message = batch_stats_flag.get_message("failed_batch", **failed_batch_para)
                            result.append(red(batch_failed_message, None, self.table_title, self.table_type))
                            self.failed_batch_numbers.append(batch)
                            red_count += 1

                        elif (test_threshold - Decimal(5)) <= abs(overall_clc_cv) <= test_threshold:
                            message = batch_stats_flag.get_message("cv", **threshold_yellow_para)
                            result.append(yellow(message, None, self.table_title, self.table_type))

            self.intra_batch_info[previous_batch] += red_count
            self.failed_batch_numbers = list(set(self.failed_batch_numbers))
            self.passed_batch_numbers = [x for x in self.total_batches if x not in self.failed_batch_numbers]

            if len(self.passed_batch_numbers) < self.required_batches:
                message = batch_stats_flag.get_message("required_batches", required_batches=self.required_batches)
                result.append(red(message, None, self.table_title, self.table_type))

        final_static_df = self.final_static_df
        inter_red_count = 0
        if not final_static_df.empty:
            all_column = final_static_df["column"]
            reported_mean = final_static_df["mean"]
            reported_cv = final_static_df["cv"]
            reported_re = final_static_df["re"]
            reported_n = final_static_df["n"]

            calculated_mean = final_static_df["calc_mean"]
            calculated_cv = final_static_df["calc_cv"]
            calculated_re = final_static_df["calc_re"]
            calculated_n = final_static_df["calc_n"]

            per_diff_mean = final_static_df["per_diff_mean"]
            per_diff_cv = final_static_df["per_diff_cv"]
            per_diff_re = final_static_df["per_diff_re"]
            per_diff_n = final_static_df["per_diff_n"]

            nominal_concentrations = final_static_df["nominal"]

            if self.found_sd:
                reported_sd = final_static_df["sd"]
                calculated_sd = final_static_df["calc_sd"]
                per_diff_sd = final_static_df["per_diff_sd"]

            if self.found_final_te:
                calculated_tes = final_static_df["calc_te"]
                per_diff_te = final_static_df["per_diff_te"]
                reported_te = final_static_df["te"]

            for index in range(len(all_column)):
                te_threshold = self.te_threshold
                column = all_column[index]

                overall_reported_mean = reported_mean[index]
                overall_reported_cv = reported_cv[index]
                overall_reported_re = reported_re[index]
                overall_reported_n = reported_n[index]

                overall_clc_mean = calculated_mean[index]
                overall_clc_cv = calculated_cv[index]
                overall_clc_re = calculated_re[index]
                overall_clc_n = calculated_n[index]

                overall_per_diff_mean = per_diff_mean[index]
                overall_per_diff_cv = per_diff_cv[index]
                overall_per_diff_re = per_diff_re[index]
                overall_per_diff_n = per_diff_n[index]

                nominal = nominal_concentrations[index]

                if self.found_final_te:
                    calculated_te = calculated_tes[index]
                    overall_per_diff_te = per_diff_te[index]
                    overall_reported_te = reported_te[index]

                if "lloq" in str(column).lower() or self.LLOQ == nominal:
                    lloq_inter_cv = overall_clc_cv
                    lloq_inter_bias = overall_clc_re
                    if self.found_final_te:
                        lloq_inter_te = calculated_te

                    test_threshold = threshold + Decimal(5)
                    te_threshold += Decimal(10)

                elif "uloq" in str(column).lower() or self.ULOQ == nominal:
                    uloq_inter_cv = overall_clc_cv
                    uloq_inter_bias = overall_clc_re
                    if self.found_final_te:
                        uloq_inter_te = calculated_te

                    test_threshold = threshold + Decimal(5)
                    te_threshold += Decimal(10)

                else:
                    lmh_inter_cvs.append(overall_clc_cv)
                    lmh_inter_biases.append(overall_clc_re)
                    if self.found_final_te:
                        lmh_inter_tes.append(calculated_te)
                    test_threshold = threshold

                count_error_para = {"column": column, "reported_value": overall_reported_n, "calc_value": overall_clc_n}

                mean_para = {"column": column, "reported_value": overall_reported_mean,
                             "calc_value": utils.format_value(overall_clc_mean)}

                cv_para = {"column": column, "reported_value": overall_reported_cv,
                           "calc_value": utils.format_value(overall_clc_cv)}

                re_para = {"column": column, "reported_value": overall_reported_re,
                           "calc_value": utils.format_value(overall_clc_re)}

                threshold_red_para = {"column": column, "threshold": test_threshold}
                threshold_yellow_para = {"column": column, "threshold": test_threshold - Decimal(5)}

                if self.found_final_te:
                    te_error_para = {"column": column, "reported_value": overall_reported_te,
                                     "calc_value": utils.format_value(calculated_te)}

                    if str(overall_per_diff_te) != "":
                        if overall_per_diff_te > valid_te_diff:
                            message = stats_flag.get_message("te_error", **te_error_para)
                            result.append(red(message, None, self.table_title, self.table_type))

                    if abs(calculated_te) > te_threshold:
                        message = stats_flag.get_message("te", column=column, te_threshold=te_threshold)
                        result.append(red(message, None, self.table_title, self.table_type))

                    elif (te_threshold - Decimal(5)) <= abs(calculated_te) <= te_threshold:
                        message = stats_flag.get_message("te", column=column, te_threshold=te_threshold - Decimal(5))
                        result.append(yellow(message, None, self.table_title, self.table_type))

                if self.found_sd:
                    overall_reported_sd = reported_sd[index]
                    overall_clc_sd = calculated_sd[index]
                    overall_per_diff_sd = per_diff_sd[index]

                    sd_para = {"column": column, "reported_value": overall_reported_sd,
                               "calc_value": utils.format_value(overall_clc_sd)}
                    if overall_per_diff_sd > valid_re_cv_diff:
                        message = stats_flag.get_message("sd_rounding", **sd_para)
                        if abs(overall_clc_sd) <= Decimal(1) and abs(overall_reported_sd) <= Decimal(1):
                            result.append(yellow(message, None, self.table_title, self.table_type))

                        else:
                            message = stats_flag.get_message("sd_error", **sd_para)
                            result.append(red(message, None, self.table_title, self.table_type))

                if overall_per_diff_n > valid_count_diff:
                    message = stats_flag.get_message("count_error", **count_error_para)
                    result.append(red(message, None, self.table_title, self.table_type))

                if overall_per_diff_mean > valid_mean_diff:
                    if abs(overall_clc_mean) <= Decimal(1) and abs(overall_reported_mean) <= Decimal(1):
                        message = stats_flag.get_message("mean_rounding", **mean_para)
                        result.append(yellow(message, None, self.table_title, self.table_type))

                    else:
                        message = stats_flag.get_message("mean_error", **mean_para)
                        result.append(red(message, None, self.table_title, self.table_type))

                if overall_per_diff_cv > valid_re_cv_diff:
                    if abs(overall_clc_cv) <= Decimal(1) and abs(overall_reported_cv) <= Decimal(1):
                        message = stats_flag.get_message("cv_rounding", **cv_para)
                        result.append(yellow(message, None, self.table_title, self.table_type))

                    else:
                        # inter_red_count += 1
                        message = stats_flag.get_message("cv_error", **cv_para)
                        result.append(red(message, None, self.table_title, self.table_type))

                if overall_per_diff_re > valid_re_cv_diff:
                    if abs(overall_clc_re) <= Decimal(1) and abs(overall_reported_re) <= Decimal(1):
                        message = stats_flag.get_message("re_rounding", **re_para)
                        result.append(yellow(message, None, self.table_title, self.table_type))

                    else:
                        # inter_red_count += 1
                        message = stats_flag.get_message("re_error", **re_para)
                        result.append(red(message, None, self.table_title, self.table_type))
                failed_batch_message = stats_flag.get_message("failed_batch", column=column)

                if abs(overall_clc_re) > test_threshold:
                    message = stats_flag.get_message("re", **threshold_red_para)
                    result.append(red(message, None, self.table_title, self.table_type))
                    result.append(red(failed_batch_message, None, self.table_title, self.table_type))
                    inter_red_count += 1
                elif (test_threshold - Decimal(5)) <= abs(overall_clc_re) <= test_threshold:
                    message = stats_flag.get_message("re", **threshold_yellow_para)
                    result.append(yellow(message, None, self.table_title, self.table_type))

                if abs(overall_clc_cv) > test_threshold:
                    message = stats_flag.get_message("cv", **threshold_red_para)
                    result.append(red(message, None, self.table_title, self.table_type))
                    result.append(red(failed_batch_message, None, self.table_title, self.table_type))
                    inter_red_count += 1

                elif (test_threshold - Decimal(5)) <= abs(overall_clc_cv) <= test_threshold:
                    message = stats_flag.get_message("cv", **threshold_yellow_para)
                    result.append(yellow(message, None, self.table_title, self.table_type))

        if self.table_subtype != "inter":
            for batch, val in self.intra_batch_info.items():
                if val == 0:
                    message = batch_stats_flag.get_message("ap_demonstrated", batch=batch)
                    self.result.append(green(message, None, self.table_title, self.table_type))

        if self.table_subtype != "intra" and self.table_subtype != "inter":
            if sum(self.intra_batch_info.values()) == 0 and inter_red_count == 0:
                message = stats_flag.get_message('ap_demonstrated', batch_count=len(self.intra_batch_info))
                self.result.append(green(message, None, self.table_title, self.table_type))

        if self.table_subtype == "intra":
            message = batch_stats_flag.get_message("no_statistics")
            self.result.append(yellow(message, None, self.table_title, self.table_type))

        elif self.table_subtype == "inter":
            for batch in set(batches):
                message = stats_flag.get_message("no_statistics", batch=batch)
                self.result.append(yellow(message, None, self.table_title, self.table_type))
        else:
            pass
        if len(result) == 0:
            result.append(green('All values within nominal range', None, self.table_title, self.table_type))

            # intra accuracy and precision
        if len(lloq_intra_cvs) > 0:
            lloq_intra_cvs = [x for x in lloq_intra_cvs if str(x).strip() != ""]
            lloq_intra_high_precision = round(max(lloq_intra_cvs), 2)
            lloq_intra_low_precision = round(min(lloq_intra_cvs), 2)
            self.intra_precisions["lloq_precision"] = [lloq_intra_low_precision, lloq_intra_high_precision]

        if len(uloq_intra_cvs) > 0:
            uloq_intra_cvs = [x for x in uloq_intra_cvs if str(x).strip() != ""]
            uloq_intra_high_precision = round(max(uloq_intra_cvs), 2)
            uloq_intra_low_precision = round(min(uloq_intra_cvs), 2)
            self.intra_precisions["uloq_precision"] = [uloq_intra_low_precision, uloq_intra_high_precision]

        if len(lmh_intra_cvs) > 0:
            lmh_intra_cvs = [x for x in lmh_intra_cvs if str(x).strip() != ""]
            lmh_intra_high_precision = round(max(lmh_intra_cvs), 2)
            lmh_intra_low_precision = round(min(lmh_intra_cvs), 2)
            self.intra_precisions["lmh_precision"] = [lmh_intra_low_precision, lmh_intra_high_precision]

        if len(lloq_intra_biases) > 0:
            lloq_intra_biases = [x for x in lloq_intra_biases if str(x).strip() != ""]
            lloq_intra_high_accuracy = round(max(lloq_intra_biases), 2)
            lloq_intra_low_accuracy = round(min(lloq_intra_biases), 2)
            self.intra_accuracies["lloq_accuracy"] = [lloq_intra_low_accuracy, lloq_intra_high_accuracy]

        if len(uloq_intra_biases) > 0:
            uloq_intra_biases = [x for x in uloq_intra_biases if str(x).strip() != ""]
            uloq_intra_high_accuracy = round(max(uloq_intra_biases), 2)
            uloq_intra_low_accuracy = round(min(uloq_intra_biases), 2)
            self.intra_accuracies["uloq_accuracy"] = [uloq_intra_low_accuracy, uloq_intra_high_accuracy]

        if len(lmh_intra_biases) > 0:
            lmh_intra_biases = [x for x in lmh_intra_biases if str(x).strip() != ""]
            lmh_intra_high_accuracy = round(max(lmh_intra_biases), 2)
            lmh_intra_low_accuracy = round(min(lmh_intra_biases), 2)
            self.intra_accuracies["lmh_accuracy"] = [lmh_intra_low_accuracy, lmh_intra_high_accuracy]

        # inter accuracy and precision
        if lloq_inter_cv is not None and str(lloq_inter_cv).strip() != "":
            self.inter_precisions["lloq_precision"] = round(lloq_inter_cv, 2)

        if uloq_inter_cv is not None and str(uloq_inter_cv).strip() != "":
            self.inter_precisions["uloq_precision"] = round(uloq_inter_cv, 2)

        if len(lmh_inter_cvs) > 0:
            lmh_inter_cvs = [x for x in lmh_inter_cvs if str(x).strip() != ""]
            lmh_inter_high_precision = round(max(lmh_inter_cvs), 2)
            lmh_inter_low_precision = round(min(lmh_inter_cvs), 2)
            self.inter_precisions["lmh_precision"] = [lmh_inter_low_precision, lmh_inter_high_precision]

        if lloq_inter_bias is not None and str(lloq_inter_bias).strip() != "":
            self.inter_accuracies["lloq_accuracy"] = round(lloq_inter_bias, 2)

        if uloq_inter_bias is not None and str(uloq_inter_bias).strip() != "":
            self.inter_accuracies["uloq_accuracy"] = round(uloq_inter_bias, 2)

        if lloq_inter_te is not None and str(lloq_inter_te).strip() != "":
            self.inter_tes["lloq_te"] = round(lloq_inter_te, 2)

        if uloq_inter_te is not None and str(uloq_inter_te).strip() != "":
            self.inter_tes["uloq_te"] = round(uloq_inter_te, 2)

        if len(lmh_inter_tes) > 0:
            lmh_inter_tes = [x for x in lmh_inter_tes if str(x).strip() != ""]
            lmh_inter_high_te = round(max(lmh_inter_tes), 2)
            lmh_inter_low_te = round(min(lmh_inter_tes), 2)
            self.inter_tes["lmh_te"] = [lmh_inter_low_te, lmh_inter_high_te]

        if len(lmh_inter_biases) > 0:
            lmh_inter_biases = [x for x in lmh_inter_biases if str(x).strip() != ""]
            lmh_inter_high_accuracy = round(max(lmh_inter_biases), 2)
            lmh_inter_low_accuracy = round(min(lmh_inter_biases), 2)
            self.inter_accuracies["lmh_accuracy"] = [lmh_inter_low_accuracy, lmh_inter_high_accuracy]

        self.result += result

    @staticmethod
    def get_column_level(column):
        if "lloq" in str(column).lower():
            return "LLOQ QC"
        elif "hqc" in str(column).lower() or "high" in str(column).lower():
            return "High QC"
        elif "lqc" in str(column).lower() or "low" in str(column).lower():
            return "LOW QC"
        elif "mqc" in str(column).lower() or "mid" in str(column).lower():
            return "MID QC"
        elif "uloq" in str(column).lower():
            return "ULOQ QC"
        else:
            return ""

    def get_error_dict(self):
        if self.orientation == "HZ":
            data_flag = FlagProperties("ap_data")
            stats_flag = FlagProperties("ap_batch_stats")
            final_stats_flag = FlagProperties("ap_stats")
            return data_flag, stats_flag, final_stats_flag


class QualityControl(Table, Flags):
    def __init__(self, parsed_table, analytes, template_type):
        Table.__init__(self, parsed_table, template_type)
        self.failed_batches = []
        self.passed_batches = []
        self.total_batches = []
        self.header_values = []
        self.ar_upper_limit = Decimal(self.threshold_values.get_message("ar_upper_threshold"))
        self.ar_lower_limit = Decimal(self.threshold_values.get_message("ar_lower_threshold"))
        self.analyte = get_analyte(analytes, self.tb_title)
        self.found_duplicate_run = False

    def validate_table_format(self):
        error_messages = self.error_messages
        required_col = {"Run ID"}
        missing_col = required_col.intersection(self.missing_col)
        if missing_col:
            message = error_messages.get_message("missing_col", col_names=", ".join(missing_col))
            self.result.append(red(message, None, self.table_title, self.table_type))
            self.valid_format = False

        conc_units = utils.find_units(self.table_df.columns)
        if conc_units == "":
            message = error_messages.get_message("missing_conc_unit")
            self.result.append(red(message, None, self.table_title, self.table_type))
            self.valid_format = False
        else:
            self.conc_unit = conc_units

    def process_table(self):
        table = self.table_df
        static_df = pd.DataFrame()
        final_df = pd.DataFrame()
        calc_static_df = pd.DataFrame()
        self.found_re = utils.find_re_column(table)
        nominal_conc = utils.find_nominal_conc(table.columns)
        try:
            try:
                table, static_df = utils.split_static_table(table)
            except:
                static_df = pd.DataFrame()

            all_run_tables, duplicate_id = utils.split_run_df(table)
            for run_table in all_run_tables:
                table = run_table["table"]
                table = utils.fill_val(table)
                table, found_re = format_calibration_table(table, nominal_conc)

                final_df = pd.concat([final_df, table]).reset_index(drop=True)

            dummy_df = final_df.copy()
            dummy_df["conc"] = dummy_df["conc"].apply(utils.parse_decimal)
            final_df["calc_re"] = ((dummy_df["conc"] - dummy_df["nominal"]) / dummy_df["nominal"]) * 100
            if self.found_re:
                final_df = final_df.fillna("")
                final_df["per_diff_re"] = utils.calculate_per_diff(final_df["re"], final_df["calc_re"])

            samples = list(dummy_df["column"].unique())
            for sample in samples:
                values = dummy_df[dummy_df["column"] == sample]["conc"].to_list()
                nominal = dummy_df[dummy_df["column"] == sample]["nominal"].to_list()
                calc_static_df = pd.concat([calc_static_df, utils.build_static_df(values, nominal[0])]).reset_index(
                    drop=True)

            if not static_df.empty:
                static_df = utils.process_static_df(static_df, self.table_type)
                static_df, self.found_sd = utils.concat_static_df(static_df, calc_static_df)
                self.found_ar = utils.find_ar_column(static_df)
                if self.found_ar:
                    static_df["ar"] = static_df["ar"].apply(utils.parse_decimal)
                    static_df["calc_ar"] = (static_df["calc_mean"] / static_df["nominal"]) * 100
                    static_df["per_diff_ar"] = abs((static_df["calc_ar"] - static_df["ar"]) / static_df["calc_ar"])

            self.data_df = final_df.fillna("")
            self.static_df = static_df
            self.found_cv = utils.find_cv_column(static_df)
            self.header_values = list(final_df["nominal"].unique())
        except Exception as e:
            message = self.error_messages.get_message("data_error")
            self.result.append(red(message, None, self.table_title, self.table_type))
            self.valid_data = False

    @staticmethod
    def get_error_dict():
        data_flag = FlagProperties("qc_data")
        stats_flag = FlagProperties("qc_stats")
        return data_flag, stats_flag

    def validate_table_data(self):
        validation_failed_flag_count = 0
        data_flag, stats_flag = self.get_error_dict()
        result = []

        valid_count_diff = Decimal(self.valid_difference_values.get_message("difference_n"))
        valid_mean_diff = Decimal(self.valid_difference_values.get_message("difference_mean"))
        valid_re_cv_diff = Decimal(self.valid_difference_values.get_message("difference_re_cv"))
        valid_te_diff = Decimal(self.valid_difference_values.get_message("difference_te"))
        valid_ar_diff = Decimal(self.valid_difference_values.get_message("difference_ar"))

        table = self.data_df
        threshold = Decimal(self.threshold_values.get_message("re_cv_threshold"))
        if "lm" in self.analysis_type:
            threshold += Decimal(5)

        all_columns = table["column"]
        batches = table["run_id"].to_list()
        all_conc = table["conc"]
        calc_res = table["calc_re"]
        if self.found_re:
            reported_res = table["re"]
            per_diff_res = table["per_diff_re"]

        self.total_batches = list(table["run_id"].unique())
        previous_batch = batches[0]
        previous_column = all_columns[0]

        total_run = 0
        failed_run = 0
        column_sample = 0
        column_sample_failed = 0
        sample = 0

        first_time = False
        column_failed = False

        if "dilution" in previous_column.lower():
            required_level = 3
            dilution = True
        else:
            required_level = 2
            dilution = False

        for index, column in enumerate(all_columns):
            batch = batches[index]
            conc = all_conc[index]
            calc_re = calc_res[index]

            if previous_column != column:
                try:
                    column_pct = (column_sample_failed / column_sample) * 100
                except ZeroDivisionError:
                    column_pct = 0
                dilution_para = {"batch": previous_batch, "pct": round(column_pct, 2)}
                level_para = {"batch": previous_batch, "column": previous_column, "required_level": required_level}

                if column_pct > 50:
                    if "dilution" in previous_column.lower():
                        message = data_flag.get_message("dilution_failed", **dilution_para)
                        result.append(yellow(message, None, self.table_title, self.table_type))

                    elif not first_time:
                        str1 = self.find_qc_level(previous_column)
                        if str(str1) == "":
                            str2 = "at a single concentration level"
                        else:
                            str2 = ""
                        column_failed_para = {"batch": previous_batch, "pct": round(column_pct, 2), "str1": str1,
                                              "str2": str2}
                        message = data_flag.get_message("failed_column", **column_failed_para)
                        result.append(red(message, None, self.table_title, self.table_type))

                        column_failed = True
                        first_time = True
                        self.failed_batches.append(previous_batch)
                if dilution:
                    if 1 < sample < required_level:
                        message = data_flag.get_message("required_level", **level_para)
                        result.append(red(message, None, self.table_title, self.table_type))
                    elif sample == 0:
                        message = data_flag.get_message("required_level", **level_para)
                        result.append(yellow(message, None, self.table_title, self.table_type))

                elif not dilution and "lloq" not in str(previous_column).lower() and "uloq" not in str(
                        previous_column).lower():
                    if sample < required_level:
                        message = data_flag.get_message("required_level", **level_para)
                        result.append(red(message, None, self.table_title, self.table_type))

                column_sample = 0
                column_sample_failed = 0
                previous_column = column
                sample = 0

                if "dilution" in previous_column.lower():
                    required_level = 3
                    dilution = True
                else:
                    required_level = 2
                    dilution = False

            if previous_batch != batch:
                if not column_failed:
                    try:
                        check_per = (failed_run / total_run) * 100
                    except ZeroDivisionError:
                        check_per = 0
                    if failed_run > (total_run - (total_run * 2 / 3)):
                        message = data_flag.get_message("failed_batch", batch=previous_batch)
                        result.append(red(message, None, self.table_title, self.table_type))
                        self.failed_batches.append(previous_batch)

                    elif int(check_per) == 33:
                        message = data_flag.get_message("accept_batch", batch=previous_batch)
                        result.append(yellow(message, None, self.table_title, self.table_type))

                total_run = 0
                failed_run = 0
                previous_batch = batch
                first_time = False
                column_failed = False

            if str(conc).strip() != "NA" and str(conc).strip() != "":
                column_sample += 1
            elif "dilution" in str(column).lower():
                column_sample += 1

            if not str(conc).isalpha() and str(conc) != "":
                sample += 1

            if not dilution and str(conc).strip() != "NA":
                total_run += 1

            if self.found_re:
                reported_re = reported_res[index]
                per_diff_re = per_diff_res[index]

                re_para = {"batch": batch, "column": column, "reported_value": reported_re,
                           "calc_value": utils.format_value(calc_re), "conc": conc}
                if str(per_diff_re) != "":
                    if per_diff_re > valid_re_cv_diff:
                        if abs(calc_re) <= Decimal(1) and abs(reported_re) <= Decimal(1):
                            message = data_flag.get_message("re_rounding", **re_para)
                            result.append(yellow(message, None, self.table_title, self.table_type))
                        else:
                            message = data_flag.get_message("re_error", **re_para)
                            result.append(red(message, None, self.table_title, self.table_type))

            if "lloq" in str(column).lower() or "uloq" in str(column).lower():
                test_threshold = threshold + Decimal(5)
            else:
                test_threshold = threshold
            re_red_para = {"batch": batch, "column": column, "conc": conc, "threshold": test_threshold}
            re_yellow_para = {"batch": batch, "column": column, "conc": conc, "threshold": test_threshold - Decimal(5)}
            if str(calc_re) != "":
                if abs(calc_re) > test_threshold:
                    message = data_flag.get_message("re", **re_red_para)
                    result.append(red(message, None, self.table_title, self.table_type))
                    column_sample_failed += 1
                    validation_failed_flag_count += 1

                    if not dilution:
                        failed_run += 1

                elif (test_threshold - Decimal(5)) <= abs(calc_re) <= test_threshold:
                    message = data_flag.get_message("re", **re_yellow_para)
                    result.append(yellow(message, None, self.table_title, self.table_type))

            if "aql" in str(conc).lower() or "bql" in str(conc).lower():
                message = data_flag.get_message("re", **re_red_para)
                result.append(red(message, None, self.table_title, self.table_type))
                failed_run += 1
                column_sample_failed += 1
                validation_failed_flag_count += 1

        try:
            column_pct = (column_sample_failed / column_sample) * 100
        except ZeroDivisionError:
            column_pct = 0
        if column_pct > 50:
            if "dilution" in previous_column.lower():
                message = data_flag.get_message("dilution_failed", batch=previous_batch, pct=round(column_pct, 2))
                result.append(yellow(message, None, self.table_title, self.table_type))
            elif not first_time:
                str1 = self.find_qc_level(previous_column)
                if str(str1) == "":
                    str2 = "at a single concentration level"
                else:
                    str2 = ""
                message = data_flag.get_message("failed_column", batch=previous_batch, pct=round(column_pct, 2),
                                                str1=str1, str2=str2)
                result.append(red(message, None, self.table_title, self.table_type))
                column_failed = True
                self.failed_batches.append(previous_batch)

            if sample < required_level:
                message = data_flag.get_message("required_level", batch=previous_batch, column=previous_column,
                                                required_level=required_level)
                result.append(red(message, None, self.table_title, self.table_type))

        if not column_failed:
            try:
                check_per = (failed_run / total_run) * 100
            except ZeroDivisionError:
                check_per = 0
            if failed_run > (total_run - (total_run * 2 / 3)):
                message = data_flag.get_message("failed_batch", batch=previous_batch)
                result.append(red(message, None, self.table_title, self.table_type))
                self.failed_batches.append(previous_batch)
            elif int(check_per) == 33:
                message = data_flag.get_message("accept_batch", batch=previous_batch)
                result.append(yellow(message, None, self.table_title, self.table_type))

        static_df = self.static_df
        ar_upper_limit = self.ar_upper_limit
        ar_lower_limit = self.ar_lower_limit
        if not static_df.empty:
            all_column = static_df["column"]
            reported_mean = static_df["mean"]
            reported_re = static_df["re"]
            reported_n = static_df["n"]

            calculated_mean = static_df["calc_mean"]
            calculated_re = static_df["calc_re"]
            calculated_n = static_df["calc_n"]

            per_diff_mean = static_df["per_diff_mean"]
            per_diff_re = static_df["per_diff_re"]
            per_diff_n = static_df["per_diff_n"]

            if self.found_cv:
                reported_cv = static_df["cv"]
                calculated_cv = static_df["calc_cv"]
                per_diff_cv = static_df["per_diff_cv"]

            if self.found_sd:
                reported_sd = static_df["sd"]
                calculated_sd = static_df["calc_sd"]
                per_diff_sd = static_df["per_diff_sd"]

            if self.found_ar:
                calculated_ars = static_df["calc_ar"]
                per_diff_ar = static_df["per_diff_ar"]
                reported_ar = static_df["ar"]

            for index, column in enumerate(all_column):

                overall_reported_mean = reported_mean[index]
                overall_reported_re = reported_re[index]
                overall_reported_n = reported_n[index]

                overall_clc_mean = calculated_mean[index]

                overall_clc_re = calculated_re[index]
                overall_clc_n = calculated_n[index]

                overall_per_diff_mean = per_diff_mean[index]
                overall_per_diff_re = per_diff_re[index]
                overall_per_diff_n = per_diff_n[index]

                mean_para = {"column": column, "reported_value": overall_reported_mean,
                             "calc_value": utils.format_value(overall_clc_mean)}

                re_para = {"column": column, "reported_value": overall_reported_re,
                           "calc_value": utils.format_value(overall_clc_re)}

                if self.found_sd:
                    overall_reported_sd = reported_sd[index]
                    overall_clc_sd = calculated_sd[index]
                    overall_per_diff_sd = per_diff_sd[index]

                    sd_para = {"column": column, "reported_value": overall_reported_sd,
                               "calc_value": utils.format_value(overall_clc_sd)}

                    if overall_per_diff_sd > valid_re_cv_diff:
                        if abs(overall_clc_sd) <= Decimal(1) and abs(overall_reported_sd) <= Decimal(1):
                            message = stats_flag.get_message("sd_rounding", **sd_para)
                            result.append(yellow(message, None, self.table_title, self.table_type))

                        else:
                            message = stats_flag.get_message("sd_error", **sd_para)
                            result.append(red(message, None, self.table_title, self.table_type))

                if "lloq" in str(column).lower() or "uloq" in str(column).lower():
                    test_threshold = threshold + Decimal(5)
                    upper_limit = ar_upper_limit + Decimal(5)
                    lower_limit = ar_lower_limit - Decimal(5)
                else:
                    test_threshold = threshold
                    upper_limit = ar_upper_limit
                    lower_limit = ar_lower_limit

                threshold_red = {"column": column, "threshold": test_threshold}
                threshold_yellow = {"column": column, "threshold": test_threshold - Decimal(5)}

                if self.found_ar:
                    calculated_ar = calculated_ars[index]
                    overall_per_diff_ar = per_diff_ar[index]
                    overall_reported_ar = reported_ar[index]

                    if overall_per_diff_ar > valid_ar_diff:
                        message = stats_flag.get_message("ar_error", column=column, reported_value=overall_reported_ar,
                                                         calc_value=utils.format_value(calculated_ar))
                        result.append(red(message, None, self.table_title, self.table_type))

                    if abs(calculated_ar) > upper_limit:
                        message = stats_flag.get_message("ar_upper", column=column, limit=upper_limit)
                        result.append(red(message, None, self.table_title, self.table_type))
                        validation_failed_flag_count += 1

                    elif (upper_limit - Decimal(5)) <= abs(calculated_ar) <= upper_limit:
                        message = stats_flag.get_message("ar_upper", column=column, limit=upper_limit - Decimal(5))
                        result.append(yellow(message, None, self.table_title, self.table_type))

                    elif abs(calculated_ar) < lower_limit:
                        message = stats_flag.get_message("ar_lower", column=column, limit=lower_limit)
                        result.append(red(message, None, self.table_title, self.table_type))
                        validation_failed_flag_count += 1

                if overall_per_diff_n > valid_count_diff:
                    message = stats_flag.get_message("count_error", column=column, reported_value=overall_reported_n,
                                                     calc_value=overall_clc_n)
                    result.append(red(message, None, self.table_title, self.table_type))

                if overall_per_diff_mean > valid_mean_diff:
                    if abs(overall_clc_mean) <= Decimal(1) and abs(overall_reported_mean) <= Decimal(1):
                        message = stats_flag.get_message("mean_rounding", **mean_para)
                        result.append(yellow(message, None, self.table_title, self.table_type))

                    else:
                        message = stats_flag.get_message("mean_error", **mean_para)
                        result.append(red(message, None, self.table_title, self.table_type))

                if self.found_cv:
                    overall_reported_cv = reported_cv[index]
                    overall_clc_cv = calculated_cv[index]
                    overall_per_diff_cv = per_diff_cv[index]
                    cv_para = {"column": column, "reported_value": overall_reported_cv,
                               "calc_value": utils.format_value(overall_clc_cv)}

                    if overall_per_diff_cv > valid_re_cv_diff:
                        if abs(overall_reported_cv) <= Decimal(1) and abs(overall_clc_cv) <= Decimal(1):
                            message = stats_flag.get_message("cv_rounding", **cv_para)
                            result.append(yellow(message, None, self.table_title, self.table_type))

                        else:
                            message = stats_flag.get_message("cv_error", **cv_para)
                            result.append(red(message, None, self.table_title, self.table_type))

                    if abs(overall_clc_cv) > test_threshold:
                        message = stats_flag.get_message("cv", **threshold_red)
                        result.append(red(message, None, self.table_title, self.table_type))
                        validation_failed_flag_count += 1

                    elif (test_threshold - Decimal(5)) <= abs(overall_clc_cv) <= test_threshold:
                        message = stats_flag.get_message("cv", **threshold_yellow)
                        result.append(yellow(message, None, self.table_title, self.table_type))

                if overall_per_diff_re > valid_re_cv_diff:
                    if abs(overall_reported_re) <= Decimal(1) and abs(overall_clc_re) <= Decimal(1):
                        message = stats_flag.get_message("re_rounding", **re_para)
                        result.append(yellow(message, None, self.table_title, self.table_type))

                    else:
                        message = stats_flag.get_message("re_error", **re_para)
                        result.append(red(message, None, self.table_title, self.table_type))

                if abs(overall_clc_re) > test_threshold:
                    message = stats_flag.get_message("re", **threshold_red)
                    result.append(red(message, None, self.table_title, self.table_type))
                    validation_failed_flag_count += 1

                elif (test_threshold - Decimal(5)) <= abs(overall_clc_re) <= test_threshold:
                    message = stats_flag.get_message("re", **threshold_yellow)
                    result.append(yellow(message, None, self.table_title, self.table_type))

        else:
            result.append(yellow(f"Could not find overall statistic", None, self.table_title, self.table_type))
        if validation_failed_flag_count == 0:
            result.append(green('All values within nominal range', None, self.table_title, self.table_type))

        self.passed_batches = [x for x in self.total_batches if x not in self.failed_batches]
        self.result += result

    @staticmethod
    def find_qc_level(column):
        column = str(column).lower()
        column_token = column.split(" ")
        str_value = list(filter(lambda x: "mid" in str(x) or "low" in str(x) or "high" in str(x), column_token))
        if len(str_value) > 0:
            return str(str_value[0]).capitalize()
        else:
            return ""


class PPDSelectivity(Table, Flags):
    def __init__(self, parsed_table, LLOQ, analytes, multiple_analyte, template_type):
        Table.__init__(self, parsed_table, template_type)
        self.LLOQ = LLOQ
        self.multiple_analyte = multiple_analyte
        self.analyte = get_analyte(analytes, self.tb_title)
        self.blank_selectivity = False
        self.is_threshold = Decimal(self.threshold_values.get_message("is_threshold"))
        self.analyte_threshold = Decimal(self.threshold_values.get_message("analyte_threshold"))
        self.analyte_is_selectivity = False
        self.analyte_selectivity = False
        self.is_selectivity = False
        self.required_sample = 6
        if self.multiple_analyte:
            self.LLOQ = LLOQ.get(self.analyte)
        if "lm" in self.analysis_type:
            self.required_sample = 10
        self.valid_required_sample = True
        self.found_duplicate_run = False
        self.valid_required_sample = True
        self.all_selectivity_info = defaultdict(list)
        self.table_df = self.table_df.loc[:, ~self.table_df.columns.duplicated()]

    def validate_table_format(self):
        required_col = {"Run ID"}
        self.check_unite()
        self.check_required_column(required_col)

    def process_table(self):
        final_data_df = pd.DataFrame()
        table_df = self.table_df
        try:
            table, stats = utils.split_ppd_selectivity(table_df)
            columns = utils.filter_conc_col(table.columns, self.conc_unit)
            samples = []

            [samples.append(x) for x in self.prev_col if str(x).strip() != "" and x not in samples]

            conc_column = [str(x).replace(y, "").strip() for x, y in zip(columns, samples)]
            if utils.parse_decimal(conc_column[0]):
                nominal_conc = utils.find_nominal_conc(conc_column)
            else:
                nominal_conc = {conc_column[0]: Decimal(0)}

            for index, sample in enumerate(samples):
                for conc_col, conc_val in nominal_conc.items():
                    data_df = pd.DataFrame(columns=["run_id", "sample", "column", "nominal", "conc"])
                    data_df["run_id"] = table["run_id"]
                    data_df["sample"] = sample
                    data_df["column"] = conc_col
                    data_df["nominal"] = conc_val
                    data_df["conc"] = table[columns[index]]
                    final_data_df = pd.concat([final_data_df, data_df]).reset_index(drop=True)
            stats_df = stats.T.reset_index(drop=True)
            stats_df.columns = stats_df.iloc[0]
            stats_df = stats_df[1:].reset_index(drop=True)
            final_data_df = pd.concat([final_data_df, stats_df], axis=1, join="inner")
            final_data_df, missing_col = utils.format_table_data_frame(final_data_df, self.table_type)
            self.found_re = utils.find_re_column(final_data_df)
            final_data_df["calc_re"] = utils.calculate_re(final_data_df["conc"].apply(utils.parse_decimal),
                                                          final_data_df["nominal"])
            if self.found_re:
                final_data_df["re"] = final_data_df["re"].apply(utils.parse_decimal)
                final_data_df["per_diff_re"] = utils.calculate_per_diff(final_data_df["re"], final_data_df["calc_re"])
            self.data_df = final_data_df
        except:
            message = self.error_messages.get_message("data_error")
            self.result.append(red(message, None, self.table_title, self.table_type))
            self.valid_data = False

    def get_selectivity_type(self, column=None):
        table_title = str(self.tb_title).lower()
        column = str(column).lower()
        if ("lloq" in table_title and "selectivity" in table_title) or "lloq" in column:
            selectivity = "LLOQ QC"
        elif ("hqc" in table_title and "selectivity" in table_title) or "hqc" in column:
            selectivity = "High QC"
        else:
            selectivity = None
        return selectivity

    def check_selectivity_type(self, column=None):
        table_title = str(self.tb_title).lower()
        column = str(column).lower()
        if (
                "lloq" in table_title and "selectivity" in table_title and 'unfortified' not in table_title and 'high' not in table_title) or "lloq" in column:
            selectivity = "lloq"
        elif ("hqc" in table_title and "selectivity" in table_title) or "hqc" in column:
            selectivity = "high"
        elif 'unfortified' in table_title or "blank" in table_title or "blank" in str(
                column).lower() or "unspiked" in str(column).lower():
            selectivity = "blank"
        else:
            selectivity = None
        return selectivity

    @staticmethod
    def bql(value):
        value = str(value).lower()
        pattern = re.compile(r"bql|blq|(?=.*<)(?=.*lloq)")
        return True if pattern.search(value) else False

    def validate_table_data(self):
        data_flag, stats_flag = self.get_error_dict()
        self.data_flag = data_flag
        table = self.data_df
        failed_sample = 0
        total_bql = 0
        not_bql = 0
        result = []
        valid_count_diff = Decimal(self.valid_difference_values.get_message("difference_n"))
        valid_mean_diff = Decimal(self.valid_difference_values.get_message("difference_mean"))
        valid_re_cv_diff = Decimal(self.valid_difference_values.get_message("difference_re_cv"))
        sm_selectivity_info = defaultdict(lambda: defaultdict(dict))
        selectivity_type = self.check_selectivity_type()
        # check sample repetition
        batch_tables, duplicate_id = utils.split_run_df(table)
        for batch_table in batch_tables:
            b_table = batch_table["table"]
            duplicate_sample = b_table[b_table["sample"].duplicated()]["sample"]
            for sample in duplicate_sample:
                message = data_flag.get_message("sample_count", sample=sample, required_sample=self.required_sample)
                result.append(red(message, None, self.table_title, self.table_type))

        try:
            unique_batches = table["run_id"].unique()
        except AttributeError:
            message = self.error_messages.get_message("data_error")
            self.result.append(red(message, None, self.table_title, self.table_type))
            return
        for u_batch in unique_batches:
            unique_column = table["column"].unique()

            if 3 <= len(table[(table["column"] == unique_column[0]) & (table["run_id"] == u_batch)][
                            "sample"].unique()) < self.required_sample:
                message = data_flag.get_message("required_sample", required_sample=self.required_sample, batch=u_batch)
                result.append(yellow(message, None, self.table_title, self.table_type))
                self.valid_required_sample = False
        threshold = Decimal(self.threshold_values.get_message("re_cv_threshold"))
        if "lm" in self.analysis_type:
            threshold += Decimal(5)

        batches = table["run_id"]
        all_samples = table["sample"].to_list()
        all_columns = table["column"]
        calc_res = table["calc_re"]
        concentrations = table["conc"]
        nominal_concentrations = table["nominal"]
        if self.found_re:
            per_diff_res = table["per_diff_re"]
            reported_res = table["re"]
        if selectivity_type == 'blank':
            selectivity_results = table["Result"]
        else:
            selectivity_results = [""] * len(concentrations)
        previous_col = all_columns[0]
        previous_batch = batches[0]
        previous_blank_batch = batches[0]
        total_samples = 0
        pass_batches = dict.fromkeys(set(batches), 0)
        for index, sample in enumerate(all_samples):
            column = all_columns[index]
            calc_re = calc_res[index]
            batch = batches[index]
            conc = concentrations[index]
            selectivity_result = selectivity_results[index]
            nominal_conc = nominal_concentrations[index]

            if selectivity_type == 'blank':
                if previous_blank_batch != batch:
                    if total_bql != 0 and (total_bql - not_bql) / total_bql < 0.8:
                        message = data_flag.get_message("blank_selectivity_failed", batch=previous_blank_batch)
                        result.append(red(message, None, self.table_title, self.table_type))
                    elif self.valid_required_sample:
                        if "lm" in self.analysis_type:
                            message = data_flag.get_message("blank_selectivity_pass", batch=previous_blank_batch)
                            result.append(green(message, None, self.table_title, self.table_type))
                            pass_batches[previous_blank_batch] += 1
                    previous_blank_batch = batch
                    total_bql = 0
                    not_bql = 0
                total_bql += 1
                if not self.bql(conc) and not self.bql(selectivity_result):
                    decimal_conc = parse_decimal(conc)
                    if decimal_conc is not None and decimal_conc > self.LLOQ:
                        not_bql += 1
                        message = data_flag.get_message("not_bql", batch=batch, sample=sample, value=conc,
                                                        unit=self.conc_unit, column=column)
                        result.append(red(message, None, self.table_title, self.table_type))
                        self.error_flag = True
            else:
                total_samples += 1
                if ("lloq" in str(column).lower() or nominal_conc == self.LLOQ) or "uloq" in str(column).lower():
                    test_threshold = threshold + Decimal(5)
                else:
                    test_threshold = threshold
                if previous_col != column:
                    if "lm" in self.analysis_type and self.valid_required_sample:
                        if (total_samples - failed_sample) / total_samples < 0.8:
                            message = data_flag.get_message("selectivity_failed", batch=previous_batch,
                                                            selectivity_type=self.get_selectivity_type(previous_col))
                            result.append(red(message, None, self.table_title, self.table_type))
                        else:
                            message = data_flag.get_message("selectivity_pass", batch=previous_batch,
                                                            selectivity_type=self.get_selectivity_type(previous_col))
                            result.append(green(message, None, self.table_title, self.table_type))
                            pass_batches[previous_batch] += 1
                    elif "sm" in self.analysis_type and self.valid_required_sample:
                        if (total_samples - failed_sample) < self.required_sample:
                            sm_selectivity_info[previous_batch][previous_col][1] = False
                        else:
                            sm_selectivity_info[previous_batch][previous_col][1] = True
                    total_samples = 0
                    failed_sample = 0
                    previous_batch = batch
                    previous_col = column
                if self.found_re:
                    per_diff_re = per_diff_res[index]
                    reported_re = reported_res[index]
                    re_para = {"batch": batch, "column": column, "sample": sample, "reported_value": reported_re,
                               "calc_value": utils.format_value(calc_re)}

                    if per_diff_re > valid_re_cv_diff:
                        if abs(calc_re) <= Decimal(1) and abs(reported_re) <= Decimal(1):
                            message = data_flag.get_message("re_rounding", **re_para)
                            result.append(yellow(message, None, self.table_title, self.table_type))

                        else:
                            message = data_flag.get_message("re_error", **re_para)
                            result.append(red(message, None, self.table_title, self.table_type))
                            self.error_flag = True

                if str(calc_re).strip() != "":
                    if abs(calc_re) > test_threshold:
                        failed_sample += 1
                        message = data_flag.get_message("re", batch=batch, sample=sample, threshold=test_threshold,
                                                        column=column, conc=conc)
                        result.append(red(message, None, self.table_title, self.table_type))
                        self.error_flag = True
                    elif (test_threshold - Decimal(5)) <= abs(calc_re) <= test_threshold:
                        message = data_flag.get_message("re", batch=batch, sample=sample,
                                                        threshold=test_threshold - Decimal(5), column=column, conc=conc)
                        result.append(yellow(message, None, self.table_title, self.table_type))

        if "lm" in self.analysis_type and self.valid_required_sample:
            if selectivity_type and selectivity_type != "blank":
                if (total_samples - failed_sample) / total_samples < 0.8:
                    message = data_flag.get_message("selectivity_failed", batch=previous_batch,
                                                    selectivity_type=self.get_selectivity_type(previous_col))
                    result.append(red(message, None, self.table_title, self.table_type))
                elif self.valid_required_sample:
                    message = data_flag.get_message("selectivity_pass", batch=previous_batch,
                                                    selectivity_type=self.get_selectivity_type(previous_col))
                    result.append(green(message, None, self.table_title, self.table_type))
                    pass_batches[previous_batch] += 1

            if total_bql != 0 and (total_bql - not_bql) / total_bql < 0.8:
                message = data_flag.get_message("blank_selectivity_failed", batch=previous_blank_batch)
                result.append(red(message, None, self.table_title, self.table_type))
            elif total_bql != 0 and self.valid_required_sample:
                message = data_flag.get_message("blank_selectivity_pass", batch=previous_blank_batch)
                result.append(green(message, None, self.table_title, self.table_type))
                pass_batches[previous_blank_batch] += 1

            num_col = len(set(all_columns))
            for batch, val in pass_batches.items():
                if val == num_col:
                    self.all_selectivity_info[batch].append(True)
                else:
                    self.all_selectivity_info[batch].append(False)

        elif "sm" in self.analysis_type:
            if (total_samples - failed_sample) < 6 or not self.valid_required_sample:
                sm_selectivity_info[previous_batch][previous_col][1] = False
            else:
                sm_selectivity_info[previous_batch][previous_col][1] = True

        static_df = self.static_df
        if not static_df.empty:
            all_column = static_df["column"]
            batches = static_df["run_id"]
            reported_mean = static_df["mean"]
            reported_cv = static_df["cv"]
            reported_re = static_df["re"]
            reported_n = static_df["n"]

            calculated_mean = static_df["calc_mean"]
            calculated_cv = static_df["calc_cv"]
            calculated_re = static_df["calc_re"]
            calculated_n = static_df["calc_n"]

            per_diff_mean = static_df["per_diff_mean"]
            per_diff_cv = static_df["per_diff_cv"]
            per_diff_re = static_df["per_diff_re"]
            per_diff_n = static_df["per_diff_n"]

            if self.found_sd:
                reported_sd = static_df["sd"]
                calculated_sd = static_df["calc_sd"]
                per_diff_sd = static_df["per_diff_sd"]
            sm_selectivity_info[batches[0]][all_column[0]][2] = True
            for index, column in enumerate(all_column):
                batch = batches[index]
                overall_reported_mean = reported_mean[index]
                overall_reported_cv = reported_cv[index]
                overall_reported_re = reported_re[index]
                overall_reported_n = reported_n[index]

                overall_clc_mean = calculated_mean[index]
                overall_clc_cv = calculated_cv[index]
                overall_clc_re = calculated_re[index]
                overall_clc_n = calculated_n[index]

                overall_per_diff_mean = per_diff_mean[index]
                overall_per_diff_cv = per_diff_cv[index]
                overall_per_diff_re = per_diff_re[index]
                overall_per_diff_n = per_diff_n[index]

                mean_para = {"column": column, "reported_value": overall_reported_mean,
                             "calc_value": utils.format_value(overall_clc_mean)}
                cv_para = {"column": column, "reported_value": overall_reported_cv,
                           "calc_value": utils.format_value(overall_clc_cv)}
                re_para = {"column": column, "reported_value": overall_reported_re,
                           "calc_value": utils.format_value(overall_clc_re)}

                if self.found_sd:
                    overall_reported_sd = reported_sd[index]
                    overall_clc_sd = calculated_sd[index]
                    overall_per_diff_sd = per_diff_sd[index]

                    sd_para = {"column": column, "reported_value": overall_reported_sd,
                               "calc_value": utils.format_value(overall_clc_sd)}

                    if overall_per_diff_sd > valid_re_cv_diff:
                        if abs(overall_clc_sd) <= Decimal(1) and abs(overall_reported_sd) <= Decimal(1):
                            message = stats_flag.get_message("sd_rounding", **sd_para)
                            result.append(yellow(message, None, self.table_title, self.table_type))

                        else:
                            message = stats_flag.get_message("sd_error", **sd_para)
                            result.append(red(message, None, self.table_title, self.table_type))
                            self.error_flag = True
                if "lloq" in str(column).lower() or "uloq" in str(column).lower():
                    test_threshold = threshold + Decimal(5)
                else:
                    test_threshold = threshold

                threshold_red = {"column": column, "threshold": test_threshold}
                threshold_yellow = {"column": column, "threshold": test_threshold - Decimal(5)}

                if overall_per_diff_n > valid_count_diff:
                    message = stats_flag.get_message("count_error", column=column, reported_value=overall_reported_n,
                                                     calc_value=overall_clc_n)
                    result.append(red(message, None, self.table_title, self.table_type))
                    self.error_flag = True
                if overall_per_diff_mean > valid_mean_diff:
                    if abs(overall_clc_mean) <= Decimal(1) and abs(overall_reported_mean) <= Decimal(1):
                        message = stats_flag.get_message("mean_rounding", **mean_para)
                        result.append(yellow(message, None, self.table_title, self.table_type))

                    else:
                        message = stats_flag.get_message("mean_error", **mean_para)
                        result.append(red(message, None, self.table_title, self.table_type))
                        self.error_flag = True
                if overall_per_diff_cv > valid_re_cv_diff:
                    if abs(overall_reported_cv) <= Decimal(1) and abs(overall_clc_cv) <= Decimal(1):
                        message = stats_flag.get_message("cv_rounding", **cv_para)
                        result.append(yellow(message, None, self.table_title, self.table_type))

                    else:
                        message = stats_flag.get_message("cv_error", **cv_para)
                        result.append(red(message, None, self.table_title, self.table_type))
                        self.error_flag = True

                if overall_per_diff_re > valid_re_cv_diff:
                    if abs(overall_reported_re) <= Decimal(1) and abs(overall_clc_re) <= Decimal(1):
                        message = stats_flag.get_message("re_rounding", **re_para)
                        result.append(yellow(message, None, self.table_title, self.table_type))

                    else:
                        message = stats_flag.get_message("re_error", **re_para)
                        result.append(red(message, None, self.table_title, self.table_type))
                        self.error_flag = True

                if abs(overall_clc_re) > test_threshold:
                    message = stats_flag.get_message("re", **threshold_red)
                    result.append(red(message, None, self.table_title, self.table_type))
                    sm_selectivity_info[batch][column][2] = False
                    self.error_flag = True

                elif (test_threshold - Decimal(5)) <= abs(overall_clc_re) <= test_threshold:
                    message = stats_flag.get_message("re", **threshold_yellow)
                    result.append(yellow(message, None, self.table_title, self.table_type))

                if abs(overall_clc_cv) > test_threshold:
                    message = stats_flag.get_message("cv", **threshold_red)
                    result.append(red(message, None, self.table_title, self.table_type))
                    sm_selectivity_info[batch][column][2] = False
                    self.error_flag = True

                elif (test_threshold - Decimal(5)) <= abs(overall_clc_cv) <= test_threshold:
                    message = stats_flag.get_message("cv", **threshold_yellow)
                    result.append(yellow(message, None, self.table_title, self.table_type))

        if not self.error_flag:
            result.append(green('All values within acceptable range', None, self.table_title, self.table_type))

        if "sm" in self.analysis_type and self.valid_required_sample:
            for batch in unique_batches:
                for col, pass_dict in sm_selectivity_info[batch].items():
                    individual_pass, overall_pass = pass_dict.get(1, True), pass_dict.get(2, True)
                    selectivity_type = self.get_selectivity_type(col)
                    if selectivity_type == "LLOQ QC":
                        if individual_pass and overall_pass:
                            self.all_selectivity_info[batch].append(True)
                            message = data_flag.get_message("sm_selectivity_pass", batch=batch,
                                                            selectivity_type=self.get_selectivity_type(col))
                            result.append(green(message, None, self.table_title, self.table_type))
                        elif (individual_pass and not overall_pass) or (not individual_pass and overall_pass):
                            self.all_selectivity_info[batch].append(False)
                            message = data_flag.get_message("sm_selectivity_partial_failed", batch=previous_batch)
                            result.append(yellow(message, None, self.table_title, self.table_type))
                        elif not individual_pass and not overall_pass:
                            self.all_selectivity_info[batch].append(False)
                            message = data_flag.get_message("sm_selectivity_failed", batch=batch)
                            result.append(red(message, None, self.table_title, self.table_type))
                            result.append(red("Less then 6 samples met acceptance criteria", None, self.table_title,
                                              self.table_type))
        if static_df.empty:
            result.append(yellow(f"Could not find overall statistics", None, self.table_title, self.table_type))
        self.result += result

    @staticmethod
    def get_error_dict():
        data_flag = FlagProperties("selectivity_data")
        stats_flag = FlagProperties("selectivity_stats")
        return data_flag, stats_flag

    def check_final_selectivity(self, all_selectivity_info: list, table_nums, table_index):
        result = []
        if (table_nums - 1) == table_index:
            selectivity_info_dict = defaultdict(list)
            for selectivity_info in all_selectivity_info:
                for batch, has_pass in selectivity_info.items():
                    selectivity_info_dict[batch].extend(has_pass)
            for batch, has_pass in selectivity_info_dict.items():
                if all(has_pass) and len(has_pass) == 3:
                    message = self.data_flag.get_message("all_selectivity_pass", batch=batch)
                    result.append(green(message, None, self.table_title, self.table_type))
        return result


class Selectivity(Table, Flags):
    def __init__(self, parsed_table, LLOQ, analytes, multiple_analyte, template_type):
        Table.__init__(self, parsed_table, template_type)
        self.LLOQ = LLOQ
        self.multiple_analyte = multiple_analyte
        self.analyte = get_analyte(analytes, self.tb_title)
        self.blank_selectivity = False
        self.is_threshold = Decimal(self.threshold_values.get_message("is_threshold"))
        self.analyte_threshold = Decimal(self.threshold_values.get_message("analyte_threshold"))
        self.analyte_is_selectivity = False
        self.analyte_selectivity = False
        self.is_selectivity = False
        self.required_sample = 6
        if self.multiple_analyte:
            self.LLOQ = LLOQ.get(self.analyte)
        if "lm" in self.analysis_type:
            self.required_sample = 10
        self.valid_required_sample = True
        self.found_duplicate_run = False
        self.valid_required_sample = True
        self.all_selectivity_info = defaultdict(list)

    def validate_table_format(self):
        if self.template == "ppd":
            self.validate_ppd_format()
        else:
            self.validate_ariadne_format()

    def validate_ariadne_format(self):
        columns = self.table_df.columns
        if "peak" in str(columns).lower() and "area" in str(columns).lower() and "analyte" in str(columns).lower() and (
                "is " in str(columns).lower() or "is_" in str(columns).lower()):
            self.analyte_is_selectivity = True
            required_col = {"Run ID", "Sample", "Analyte Peak Area", "LLOQ Analyte Peak Area", "IS Peak Area",
                            "STDs and QCs Mean IS Peak Area", "%Analyte Response", "%IS Response "}

        elif "peak" in str(columns).lower() and "area" in str(columns).lower() and "analyte" in str(columns).lower():
            self.analyte_selectivity = True
            required_col = {"Run ID", "Sample", "Analyte Peak Area", "LLOQ Analyte Peak Area", "%Analyte Response"}

        elif "peak" in str(columns).lower() and "area" in str(columns).lower() and (
                "is " in str(columns).lower() or "is_" in str(columns).lower()):
            self.is_selectivity = True
            required_col = {"Run ID", "Sample", "IS Peak Area", "STDs and QCs Mean IS Peak Area", "%IS Response "}

        else:
            required_col = {"Run ID", "Sample"}
            self.check_unite()

        self.check_required_column(required_col)

    def validate_ppd_format(self):
        required_col = {"Run ID"}
        self.check_unite()
        self.check_required_column(required_col)

    def process_peak_area_table(self):
        table = utils.fill_val(self.table_df)
        try:
            if self.analyte_selectivity or self.analyte_is_selectivity:
                table["analyte_peak_area"] = table["analyte_peak_area"].apply(utils.parse_decimal)
                table["mean_analyte_area"] = table["mean_analyte_area"].apply(utils.parse_decimal)
                table["per_analyte_response"] = table["per_analyte_response"].apply(utils.parse_decimal)
                table["calc_per_analyte_response"] = (table["analyte_peak_area"] / table["mean_analyte_area"]) * 100

            if self.is_selectivity or self.analyte_is_selectivity:
                table["is_peak_area"] = table["is_peak_area"].apply(utils.parse_decimal)
                table["mean_is_area"] = table["mean_is_area"].apply(utils.parse_decimal)
                table["per_is_response"] = table["per_is_response"].apply(utils.parse_decimal)
                table["calc_per_is_response"] = (table["is_peak_area"] / table["mean_is_area"]) * 100
        except:
            raise Exception
        return table

    def process_table(self):
        table = self.table_df
        static_df = pd.DataFrame()
        calc_static_df = pd.DataFrame()

        try:
            if self.analyte_is_selectivity or self.analyte_selectivity or self.is_selectivity:
                table = self.process_peak_area_table()
                self.data_df = table
            else:
                nominal_column = []
                table, static_df = utils.split_static_table(table)
                all_run_tables, duplicate_id = utils.split_run_df(table)
                final_df = pd.DataFrame()
                self.found_re = utils.find_re_column(table)
                column = list(table.columns)
                column = [re.sub(r"\t\n", "", str(x)) for x in column]
                table_column = list(filter(lambda x: self.conc_unit in str(x), column))
                table.columns = column
                nominal_conc = utils.find_nominal_conc(column)
                if nominal_conc is None:
                    nominal_column = list(filter(lambda x: "nominal" in str(x).lower(), table_column))
                    try:
                        table_column.remove(nominal_column[0])
                    except (ValueError, IndexError):
                        pass

                for run_table in all_run_tables:
                    run_table = run_table["table"]
                    run_table = utils.fill_val(run_table)
                    samples = run_table["sample"].unique()
                    count = 0
                    for sample in samples:
                        if nominal_conc is not None:
                            for key in nominal_conc:
                                sample_df = pd.DataFrame(
                                    columns=["assay_date", "run_id", "sample", "column", "nominal", "conc", "re"])
                                try:
                                    sample_df["assay_date"] = run_table[run_table["sample"] == sample]["assay_date"]
                                except KeyError:
                                    pass
                                sample_df["run_id"] = run_table[run_table["sample"] == sample]["run_id"]
                                sample_df["sample"] = sample
                                sample_df["column"] = key
                                sample_df["nominal"] = nominal_conc[key]
                                values = run_table[run_table["sample"] == sample][key]
                                sample_df["conc"] = values
                                if self.found_re:
                                    run_col = list(run_table.columns)
                                    re_cols = list(filter(lambda x: "re" in str(x), run_col))
                                    if len(re_cols) > 1:
                                        column_index = run_table.columns.get_loc(key)
                                        run_table.columns.values[column_index + 1] = f"re_{count}"

                                        sample_df["re"] = run_table[run_table["sample"] == sample][f"re_{count}"].apply(
                                            utils.parse_decimal)

                                        count += 1
                                    else:
                                        sample_df["re"] = run_table[run_table["sample"] == sample][f"re"].apply(
                                            utils.parse_decimal)
                                final_df = pd.concat([final_df, sample_df]).reset_index(drop=True)
                        else:
                            sample_df = pd.DataFrame(
                                columns=["assay_date", "run_id", "sample", "column", "nominal", "conc", "re"])

                            if len(nominal_column) > 0:
                                sample_df["nominal"] = run_table[run_table["sample"] == sample][
                                    nominal_column[0]].apply(utils.parse_decimal)

                            else:
                                sample_df["nominal"] = None

                            sample_df["assay_date"] = run_table[run_table["sample"] == sample]["assay_date"]
                            sample_df["run_id"] = run_table[run_table["sample"] == sample]["run_id"]
                            sample_df["sample"] = sample
                            sample_df["column"] = table_column[count]
                            sample_df["conc"] = run_table[run_table["sample"] == sample][table_column[count]]
                            if self.found_re:
                                run_col = list(run_table.columns)
                                re_cols = list(filter(lambda x: "re" in str(x), run_col))
                                if len(re_cols) > 1:
                                    column_index = run_table.columns.get_loc(table_column[count])
                                    run_table.columns.values[column_index + 1] = f"re_{count}"

                                    sample_df["re"] = run_table[run_table["sample"] == sample][f"re_{count}"].apply(
                                        utils.parse_decimal)

                                    count += 1
                                else:
                                    sample_df["re"] = run_table[run_table["sample"] == sample][f"re"].apply(
                                        utils.parse_decimal)
                            final_df = pd.concat([final_df, sample_df]).reset_index(drop=True)
                dummy_df = final_df.copy()
                dummy_df["conc"] = dummy_df["conc"].apply(utils.parse_decimal)
                final_df["calc_re"] = utils.calculate_re(dummy_df["conc"], dummy_df["nominal"])

                if not static_df.empty:
                    column = dummy_df["column"].unique()
                    for col in column:
                        values = dummy_df[dummy_df["column"] == col]["conc"].to_list()
                        nominal = dummy_df[dummy_df["column"] == col]["nominal"].to_list()[0]
                        calc_static_df = pd.concat([calc_static_df, utils.build_static_df(values, nominal)])
                    static_df = utils.process_static_df(static_df, self.table_type)
                    static_df.insert(loc=0, column="run_id", value=final_df["run_id"][0])

                    # Concat reported and calculated static dataframes
                    static_df, self.found_sd = utils.concat_static_df(static_df, calc_static_df)

                if self.found_re:
                    final_df["per_diff_re"] = utils.calculate_per_diff(final_df["re"], final_df["calc_re"])
                self.data_df = final_df.fillna("").sort_values(["run_id", "column", "sample"]).reset_index(drop=True)
                self.static_df = static_df
        except Exception as e:
            message = self.error_messages.get_message("data_error")
            self.result.append(red(message, None, self.table_title, self.table_type))
            self.valid_data = False

    def validate_peak_area_table(self, data_flag):
        result = []
        table = self.data_df
        table_type = self.table_type
        table_title = self.table_title
        is_threshold = self.is_threshold
        analyte_threshold = self.analyte_threshold
        run_ids = table["run_id"]
        all_samples = table["sample"]
        per_is_responses = []
        per_analyte_responses = []
        self.found_duplicate_run = False
        sm_selectivity_info = defaultdict(dict)
        if self.analyte_selectivity or self.analyte_is_selectivity:
            per_analyte_responses = table["calc_per_analyte_response"]
        total_samples = len(all_samples)
        if self.is_selectivity or self.analyte_is_selectivity:
            per_is_responses = table["calc_per_is_response"]
        previous_batch = run_ids[0]
        analyte_failed = 0
        is_failed = 0
        for index, run in enumerate(run_ids):
            sample = all_samples[index]

            if "sm" in self.analysis_type and self.valid_required_sample:
                if previous_batch != run:
                    if self.analyte_selectivity or self.analyte_is_selectivity:
                        if (total_samples - analyte_failed) < self.required_sample:
                            sm_selectivity_info[previous_batch]["analyte"] = False
                        else:
                            sm_selectivity_info[previous_batch]["analyte"] = True
                    if self.is_selectivity or self.analyte_is_selectivity:
                        if (total_samples - is_failed) < self.required_sample:
                            sm_selectivity_info[previous_batch]["is"] = False
                        else:
                            sm_selectivity_info[previous_batch]["is"] = True
                    previous_batch = run
                    analyte_failed = 0
                    is_failed = 0
            if self.is_selectivity or self.analyte_is_selectivity:
                is_response = per_is_responses[index]
                if is_response > is_threshold:
                    result.append(
                        red(f"[Batch {run}] Sample: {sample} %IS Response value above {str(is_threshold)}%", None,
                            table_title, table_type))
                    self.error_flag = True
                    is_failed += 1
            if self.analyte_selectivity or self.analyte_is_selectivity:
                analyte_response = per_analyte_responses[index]
                if analyte_response > analyte_threshold:
                    result.append(
                        red(f"[Batch {run}] Sample: {sample} %Analyte Response value above {str(analyte_threshold)}%",
                            None, table_title, table_type))
                    self.error_flag = True
                    analyte_failed += 1
                elif (analyte_threshold - Decimal(5)) <= analyte_response <= analyte_threshold:
                    result.append(yellow(
                        f"[Batch {run}] Sample: {sample} %Analyte Response value above {str(analyte_threshold - Decimal(5))}%",
                        None, table_title, table_type))

        if "sm" in self.analysis_type and self.valid_required_sample:
            if self.analyte_selectivity or self.analyte_is_selectivity:
                if (total_samples - analyte_failed) < self.required_sample:
                    sm_selectivity_info[previous_batch]["analyte"] = False
                else:
                    sm_selectivity_info[previous_batch]["analyte"] = True
            if self.is_selectivity or self.analyte_is_selectivity:
                if (total_samples - is_failed) < self.required_sample:
                    sm_selectivity_info[previous_batch]["is"] = False
                else:
                    sm_selectivity_info[previous_batch]["is"] = True

        if not self.error_flag:
            result.append(green(f"All values within acceptable range", None, table_title, table_type))

        for batch, selectivity_info in sm_selectivity_info.items():
            for selectivity, has_pass in selectivity_info.items():
                if selectivity == "analyte":
                    if has_pass:
                        self.all_selectivity_info[batch].append(True)
                        message = data_flag.get_message("analyte_selectivity_pass", batch=previous_batch)
                        result.append(green(message, None, self.table_title, self.table_type))
                    else:
                        self.all_selectivity_info[batch].append(False)
                        message = data_flag.get_message("sm_selectivity_partial_failed", batch=previous_batch)
                        result.append(yellow(message, None, self.table_title, self.table_type))
                elif selectivity == "is":
                    if has_pass:
                        self.all_selectivity_info[batch].append(True)
                        message = data_flag.get_message("is_selectivity_pass", batch=previous_batch)
                        result.append(green(message, None, self.table_title, self.table_type))
                    else:
                        self.all_selectivity_info[batch].append(False)
                        message = data_flag.get_message("sm_selectivity_partial_failed", batch=previous_batch)
                        result.append(yellow(message, None, self.table_title, self.table_type))
        return result

    def get_selectivity_type(self, column=None):
        table_title = str(self.table_title).lower()
        column = str(column).lower()
        selectivity = ""
        if ("lloq" in table_title and "selectivity" in table_title) or "lloq" in column:
            selectivity = "LLOQ QC"
        elif ("hqc" in table_title and "selectivity" in table_title) or "hqc" in column:
            selectivity = "High QC"
        else:
            selectivity = None
        return selectivity

    def validate_table_data(self):
        valid_count_diff = Decimal(self.valid_difference_values.get_message("difference_n"))
        valid_mean_diff = Decimal(self.valid_difference_values.get_message("difference_mean"))
        valid_re_cv_diff = Decimal(self.valid_difference_values.get_message("difference_re_cv"))
        sm_selectivity_info = defaultdict(lambda: defaultdict(dict))
        failed_sample = 0
        total_bql = 0
        not_bql = 0
        data_flag, stats_flag = self.get_error_dict()
        self.data_flag = data_flag
        result = []
        table = self.data_df
        # check sample repetition
        batch_tables, duplicate_id = utils.split_run_df(self.table_df)
        for batch_table in batch_tables:
            b_table = batch_table["table"]
            duplicate_sample = b_table[b_table["sample"].duplicated()]["sample"]
            for sample in duplicate_sample:
                message = data_flag.get_message("sample_count", sample=sample, required_sample=self.required_sample)
                result.append(red(message, None, self.table_title, self.table_type))

        if self.analyte_is_selectivity or self.is_selectivity or self.analyte_selectivity:
            batch = self.table_df["run_id"].unique()[0]
            if len(self.table_df["sample"].to_list()) < self.required_sample:
                message = data_flag.get_message("required_sample", required_sample=self.required_sample, batch=batch)
                result.append(red(message, None, self.table_title, self.table_type))
                self.valid_required_sample = False
            result += self.validate_peak_area_table(data_flag)
        else:
            unique_batches = table["run_id"].unique()
            for u_batch in unique_batches:
                unique_column = table["column"].unique()

                if 3 <= len(table[(table["column"] == unique_column[0]) & (table["run_id"] == u_batch)][
                                "sample"].unique()) < self.required_sample:
                    message = data_flag.get_message("required_sample", required_sample=self.required_sample,
                                                    batch=u_batch)
                    result.append(yellow(message, None, self.table_title, self.table_type))
                    self.valid_required_sample = False
            threshold = Decimal(self.threshold_values.get_message("re_cv_threshold"))
            if "lm" in self.analysis_type:
                threshold += Decimal(5)

            batches = table["run_id"]
            all_samples = table["sample"].to_list()
            all_columns = table["column"]
            calc_res = table["calc_re"]
            concentrations = table["conc"]
            nominal_concentrations = table["nominal"]
            if self.found_re:
                per_diff_res = table["per_diff_re"]
                reported_res = table["re"]
            previous_col = all_columns[0]
            previous_batch = batches[0]
            previous_blank_batch = batches[0]
            total_samples = 0
            pass_batches = dict.fromkeys(set(batches), 0)
            for index, sample in enumerate(all_samples):
                column = all_columns[index]
                calc_re = calc_res[index]
                batch = batches[index]
                conc = concentrations[index]
                nominal_conc = nominal_concentrations[index]

                if "blank" in str(column).lower() or "unspiked" in str(column).lower() or (
                        "blank" in self.table_title.lower() and "selectivity" in self.table_title.lower()):
                    if previous_blank_batch != batch:
                        if total_bql != 0 and (total_bql - not_bql) / total_bql < 0.8:
                            message = data_flag.get_message("blank_selectivity_failed", batch=previous_blank_batch)
                            result.append(red(message, None, self.table_title, self.table_type))
                        elif self.valid_required_sample:
                            if "lm" in self.analysis_type:
                                message = data_flag.get_message("blank_selectivity_pass", batch=previous_blank_batch)
                                result.append(green(message, None, self.table_title, self.table_type))
                                pass_batches[previous_blank_batch] += 1
                        previous_blank_batch = batch
                        total_bql = 0
                        not_bql = 0
                    total_bql += 1
                    if "bql" not in str(conc).lower():
                        decimal_conc = parse_decimal(conc)
                        if decimal_conc is not None and decimal_conc > self.LLOQ:
                            not_bql += 1
                            message = data_flag.get_message("not_bql", batch=batch, sample=sample, value=conc,
                                                            unit=self.conc_unit, column=column)
                            result.append(red(message, None, self.table_title, self.table_type))
                            self.error_flag = True
                else:
                    total_samples += 1
                    if ("lloq" in str(column).lower() or nominal_conc == self.LLOQ) or "uloq" in str(column).lower():
                        test_threshold = threshold + Decimal(5)
                    else:
                        test_threshold = threshold
                    if previous_col != column:
                        if "lm" in self.analysis_type and self.valid_required_sample:
                            if (total_samples - failed_sample) / total_samples < 0.8:
                                message = data_flag.get_message("selectivity_failed", batch=previous_batch,
                                                                selectivity_type=self.get_selectivity_type(
                                                                    previous_col))
                                result.append(red(message, None, self.table_title, self.table_type))
                            else:
                                message = data_flag.get_message("selectivity_pass", batch=previous_batch,
                                                                selectivity_type=self.get_selectivity_type(
                                                                    previous_col))
                                result.append(green(message, None, self.table_title, self.table_type))
                                pass_batches[previous_batch] += 1
                        elif "sm" in self.analysis_type and self.valid_required_sample:
                            if (total_samples - failed_sample) < self.required_sample:
                                sm_selectivity_info[previous_batch][previous_col][1] = False
                            else:
                                sm_selectivity_info[previous_batch][previous_col][1] = True
                        total_samples = 0
                        failed_sample = 0
                        previous_batch = batch
                        previous_col = column
                    if self.found_re:
                        per_diff_re = per_diff_res[index]
                        reported_re = reported_res[index]
                        re_para = {"batch": batch, "column": column, "sample": sample, "reported_value": reported_re,
                                   "calc_value": utils.format_value(calc_re)}

                        if per_diff_re > valid_re_cv_diff:
                            if abs(calc_re) <= Decimal(1) and abs(reported_re) <= Decimal(1):
                                message = data_flag.get_message("re_rounding", **re_para)
                                result.append(yellow(message, None, self.table_title, self.table_type))

                            else:
                                message = data_flag.get_message("re_error", **re_para)
                                result.append(red(message, None, self.table_title, self.table_type))
                                self.error_flag = True

                    if str(calc_re).strip() != "":
                        if abs(calc_re) > test_threshold:
                            failed_sample += 1
                            message = data_flag.get_message("re", batch=batch, sample=sample, threshold=test_threshold,
                                                            column=column, conc=conc)
                            result.append(red(message, None, self.table_title, self.table_type))
                            self.error_flag = True
                        elif (test_threshold - Decimal(5)) <= abs(calc_re) <= test_threshold:
                            message = data_flag.get_message("re", batch=batch, sample=sample,
                                                            threshold=test_threshold - Decimal(5), column=column,
                                                            conc=conc)
                            result.append(yellow(message, None, self.table_title, self.table_type))

            if "lm" in self.analysis_type and self.valid_required_sample:
                if self.get_selectivity_type(previous_col) is not None:
                    if (total_samples - failed_sample) / total_samples < 0.8:
                        message = data_flag.get_message("selectivity_failed", batch=previous_batch,
                                                        selectivity_type=self.get_selectivity_type(previous_col))
                        result.append(red(message, None, self.table_title, self.table_type))
                    elif self.valid_required_sample:
                        message = data_flag.get_message("selectivity_pass", batch=previous_batch,
                                                        selectivity_type=self.get_selectivity_type(previous_col))
                        result.append(green(message, None, self.table_title, self.table_type))
                        pass_batches[previous_batch] += 1

                if total_bql != 0 and (total_bql - not_bql) / total_bql < 0.8:
                    message = data_flag.get_message("blank_selectivity_failed", batch=previous_blank_batch)
                    result.append(red(message, None, self.table_title, self.table_type))
                elif total_bql != 0 and self.valid_required_sample:
                    message = data_flag.get_message("blank_selectivity_pass", batch=previous_blank_batch)
                    result.append(green(message, None, self.table_title, self.table_type))
                    pass_batches[previous_blank_batch] += 1

                num_col = len(set(all_columns))
                for batch, val in pass_batches.items():
                    if val == num_col:
                        # message = data_flag.get_message("all_selectivity_pass", batch=batch)
                        # result.append(green(message, None, self.table_title, self.table_type))
                        self.all_selectivity_info[batch].append(True)
                    else:
                        self.all_selectivity_info[batch].append(False)

            elif "sm" in self.analysis_type:
                if (total_samples - failed_sample) < 6 or not self.valid_required_sample:
                    sm_selectivity_info[previous_batch][previous_col][1] = False
                else:
                    sm_selectivity_info[previous_batch][previous_col][1] = True

            static_df = self.static_df
            if not static_df.empty:
                all_column = static_df["column"]
                batches = static_df["run_id"]
                reported_mean = static_df["mean"]
                reported_cv = static_df["cv"]
                reported_re = static_df["re"]
                reported_n = static_df["n"]

                calculated_mean = static_df["calc_mean"]
                calculated_cv = static_df["calc_cv"]
                calculated_re = static_df["calc_re"]
                calculated_n = static_df["calc_n"]

                per_diff_mean = static_df["per_diff_mean"]
                per_diff_cv = static_df["per_diff_cv"]
                per_diff_re = static_df["per_diff_re"]
                per_diff_n = static_df["per_diff_n"]

                if self.found_sd:
                    reported_sd = static_df["sd"]
                    calculated_sd = static_df["calc_sd"]
                    per_diff_sd = static_df["per_diff_sd"]
                sm_selectivity_info[batches[0]][all_column[0]][2] = True
                for index, column in enumerate(all_column):
                    batch = batches[index]
                    overall_reported_mean = reported_mean[index]
                    overall_reported_cv = reported_cv[index]
                    overall_reported_re = reported_re[index]
                    overall_reported_n = reported_n[index]

                    overall_clc_mean = calculated_mean[index]
                    overall_clc_cv = calculated_cv[index]
                    overall_clc_re = calculated_re[index]
                    overall_clc_n = calculated_n[index]

                    overall_per_diff_mean = per_diff_mean[index]
                    overall_per_diff_cv = per_diff_cv[index]
                    overall_per_diff_re = per_diff_re[index]
                    overall_per_diff_n = per_diff_n[index]

                    mean_para = {"column": column, "reported_value": overall_reported_mean,
                                 "calc_value": utils.format_value(overall_clc_mean)}
                    cv_para = {"column": column, "reported_value": overall_reported_cv,
                               "calc_value": utils.format_value(overall_clc_cv)}
                    re_para = {"column": column, "reported_value": overall_reported_re,
                               "calc_value": utils.format_value(overall_clc_re)}

                    if self.found_sd:
                        overall_reported_sd = reported_sd[index]
                        overall_clc_sd = calculated_sd[index]
                        overall_per_diff_sd = per_diff_sd[index]

                        sd_para = {"column": column, "reported_value": overall_reported_sd,
                                   "calc_value": utils.format_value(overall_clc_sd)}

                        if overall_per_diff_sd > valid_re_cv_diff:
                            if abs(overall_clc_sd) <= Decimal(1) and abs(overall_reported_sd) <= Decimal(1):
                                message = stats_flag.get_message("sd_rounding", **sd_para)
                                result.append(yellow(message, None, self.table_title, self.table_type))

                            else:
                                message = stats_flag.get_message("sd_error", **sd_para)
                                result.append(red(message, None, self.table_title, self.table_type))
                                self.error_flag = True
                    if "lloq" in str(column).lower() or "uloq" in str(column).lower():
                        test_threshold = threshold + Decimal(5)
                    else:
                        test_threshold = threshold

                    threshold_red = {"column": column, "threshold": test_threshold}
                    threshold_yellow = {"column": column, "threshold": test_threshold - Decimal(5)}

                    if overall_per_diff_n > valid_count_diff:
                        message = stats_flag.get_message("count_error", column=column,
                                                         reported_value=overall_reported_n, calc_value=overall_clc_n)
                        result.append(red(message, None, self.table_title, self.table_type))
                        self.error_flag = True
                    if overall_per_diff_mean > valid_mean_diff:
                        if abs(overall_clc_mean) <= Decimal(1) and abs(overall_reported_mean) <= Decimal(1):
                            message = stats_flag.get_message("mean_rounding", **mean_para)
                            result.append(yellow(message, None, self.table_title, self.table_type))

                        else:
                            message = stats_flag.get_message("mean_error", **mean_para)
                            result.append(red(message, None, self.table_title, self.table_type))
                            self.error_flag = True
                    if overall_per_diff_cv > valid_re_cv_diff:
                        if abs(overall_reported_cv) <= Decimal(1) and abs(overall_clc_cv) <= Decimal(1):
                            message = stats_flag.get_message("cv_rounding", **cv_para)
                            result.append(yellow(message, None, self.table_title, self.table_type))

                        else:
                            message = stats_flag.get_message("cv_error", **cv_para)
                            result.append(red(message, None, self.table_title, self.table_type))
                            self.error_flag = True

                    if overall_per_diff_re > valid_re_cv_diff:
                        if abs(overall_reported_re) <= Decimal(1) and abs(overall_clc_re) <= Decimal(1):
                            message = stats_flag.get_message("re_rounding", **re_para)
                            result.append(yellow(message, None, self.table_title, self.table_type))

                        else:
                            message = stats_flag.get_message("re_error", **re_para)
                            result.append(red(message, None, self.table_title, self.table_type))
                            self.error_flag = True

                    if abs(overall_clc_re) > test_threshold:
                        message = stats_flag.get_message("re", **threshold_red)
                        result.append(red(message, None, self.table_title, self.table_type))
                        sm_selectivity_info[batch][column][2] = False
                        self.error_flag = True

                    elif (test_threshold - Decimal(5)) <= abs(overall_clc_re) <= test_threshold:
                        message = stats_flag.get_message("re", **threshold_yellow)
                        result.append(yellow(message, None, self.table_title, self.table_type))

                    if abs(overall_clc_cv) > test_threshold:
                        message = stats_flag.get_message("cv", **threshold_red)
                        result.append(red(message, None, self.table_title, self.table_type))
                        sm_selectivity_info[batch][column][2] = False
                        self.error_flag = True

                    elif (test_threshold - Decimal(5)) <= abs(overall_clc_cv) <= test_threshold:
                        message = stats_flag.get_message("cv", **threshold_yellow)
                        result.append(yellow(message, None, self.table_title, self.table_type))

            if not self.error_flag:
                result.append(green('All values within acceptable range', None, self.table_title, self.table_type))

            if "sm" in self.analysis_type and self.valid_required_sample:
                for batch in unique_batches:
                    for col, pass_dict in sm_selectivity_info[batch].items():
                        individual_pass, overall_pass = pass_dict.get(1, True), pass_dict.get(2, True)
                        selectivity_type = self.get_selectivity_type(col)
                        if selectivity_type == "LLOQ QC":
                            if individual_pass and overall_pass:
                                self.all_selectivity_info[batch].append(True)
                                message = data_flag.get_message("sm_selectivity_pass", batch=batch,
                                                                selectivity_type=self.get_selectivity_type(col))
                                result.append(green(message, None, self.table_title, self.table_type))
                            elif (individual_pass and not overall_pass) or (not individual_pass and overall_pass):
                                self.all_selectivity_info[batch].append(False)
                                message = data_flag.get_message("sm_selectivity_partial_failed", batch=previous_batch)
                                result.append(yellow(message, None, self.table_title, self.table_type))
                            elif not individual_pass and not overall_pass:
                                self.all_selectivity_info[batch].append(False)
                                message = data_flag.get_message("sm_selectivity_failed", batch=batch)
                                result.append(red(message, None, self.table_title, self.table_type))
                                result.append(red("Less then 6 samples met acceptance criteria", None, self.table_title,
                                                  self.table_type))
            if static_df.empty:
                result.append(yellow(f"Could not find overall statistics", None, self.table_title, self.table_type))
        self.result += result

    @staticmethod
    def get_error_dict():
        data_flag = FlagProperties("selectivity_data")
        stats_flag = FlagProperties("selectivity_stats")
        return data_flag, stats_flag

    def check_final_selectivity(self, all_selectivity_info: list, table_nums, table_index):
        result = []
        if (table_nums - 1) == table_index:
            selectivity_info_dict = defaultdict(list)
            for selectivity_info in all_selectivity_info:
                for batch, has_pass in selectivity_info.items():
                    selectivity_info_dict[batch].extend(has_pass)
            for batch, has_pass in selectivity_info_dict.items():
                if all(has_pass) and len(has_pass) == 3:
                    message = self.data_flag.get_message("all_selectivity_pass", batch=batch)
                    result.append(green(message, None, self.table_title, self.table_type))
        return result


class RegressionModel(Table, Flags):
    def __init__(self, parsed_table, analytes, multiple_analyte, template_type):
        Table.__init__(self, parsed_table, template_type)
        self.failed_batches = []
        self.total_batches = []
        self.pass_batches = []
        self.multiple_analyte = multiple_analyte
        self.r_2_threshold = Decimal(self.threshold_values.get_message("r_2_threshold"))
        self.r_threshold = Decimal(self.threshold_values.get_message("r_2_threshold"))
        if self.analyte is None:
            self.analyte = get_analyte(analytes, self.tb_title)
        self.r_column = False

    def validate_table_format(self):
        error_messages = self.error_messages
        missing_col = set()
        required_col = {"Run ID", "R-Square"}
        missing_r2_col = required_col.intersection(self.missing_col)
        if 'R-Square' in missing_r2_col:
            required_col = {"R"}
            missing_r_col = required_col.intersection(self.missing_col)
            if "R" not in missing_r_col:
                self.r_column = True
            else:
                missing_col = missing_r2_col.union(missing_r_col)
        if missing_col:
            message = error_messages.get_message("missing_col", col_names=", ".join(missing_col))
            self.result.append(red(message, None, self.table_title, self.table_type))
            self.valid_format = False

    def process_table(self):
        try:
            static_df = pd.DataFrame()
            table = self.table_df
            table, static_df = utils.split_static_table(table)
            if self.r_column:
                table["r"] = table["r"].apply(utils.parse_decimal)
            else:
                table["r_2"] = table["r_2"].apply(utils.parse_decimal)
            self.data_df = fill_val(table)
        except:
            message = self.error_messages.get_message("data_error")
            self.result.append(red(message, None, self.table_title, self.table_type))
            self.valid_data = False

    def validate_table_data(self, options):
        result = []
        if self.multiple_analyte:
            analyte_str = f'for analyte {self.analyte}'
        else:
            analyte_str = ''

        data_flag = self.get_error_dict()
        table = self.data_df

        if self.r_column:
            if table['r'].empty or table['r'].dropna().empty:
                result.append(red('No R Values Found', None, self.table_title, self.table_type))
            threshold = self.r_threshold
            r_column = "r"
            failed_msg = "r_failed"
            pass_msg = "r_pass"
        else:
            if table['r_2'].empty or table['r_2'].dropna().empty:
                result.append(red('No R Squared Values Found', None, self.table_title, self.table_type))
            threshold = self.r_2_threshold
            r_column = 'r_2'
            failed_msg = "r2_failed"
            pass_msg = "r2_pass"

        self.total_batches = list(table['run_id'])
        for index, row in table.iterrows():
            r_val = row[r_column]
            batch = row["run_id"]
            try:
                if r_val < threshold:
                    message = data_flag.get_message(failed_msg, batch=batch, analyte_str=analyte_str)
                    result.append(red(message, None, self.table_title, self.table_type))
                    self.failed_batches.append(batch)
            except TypeError:
                continue

        if len(result) == 0:
            message = data_flag.get_message(pass_msg, analyte_str=analyte_str)
            result.append(green(message, None, self.table_title, self.table_type))

        # Test regression model
        headers = list(table.columns)
        regression_model = options.get('regression_model', None)
        analyte_col = list(filter(lambda x: "analyte" in str(x).lower(), headers))

        if 'assay_date' in headers:
            headers.remove('assay_date')
        if 'run_id' in headers:
            headers.remove('run_id')
        if len(analyte_col) > 0:
            headers.remove(analyte_col[0])

        index_header = headers.index(r_column)
        header_len = len(headers[:index_header])

        if regression_model is not None and str(regression_model).strip() != "":
            if header_len == 2:
                if "linear" in regression_model.lower():
                    result.append(
                        green(f'Linear model confirmed {analyte_str}', None, self.table_title, self.table_type))
                else:
                    message = data_flag.get_message("linear_model_mismatch", regression_model=regression_model,
                                                    analyte_str=analyte_str)
                    result.append(red(message, None, self.table_title, self.table_type))
            elif header_len == 3:
                if "quadratic" in regression_model.lower():
                    result.append(
                        green(f'Quadratic model confirmed {analyte_str}', None, self.table_title, self.table_type))
                else:
                    message = data_flag.get_message("quadratic_model_mismatch", regression_model=regression_model,
                                                    analyte_str=analyte_str)
                    result.append(red(message, None, self.table_title, self.table_type))
            elif header_len == 4:
                if "4pl" in regression_model.lower():
                    result.append(green(f'4PL model confirmed {analyte_str}', None, self.table_title, self.table_type))
                else:
                    message = data_flag.get_message("4pl_model_mismatch", regression_model=regression_model,
                                                    analyte_str=analyte_str)
                    result.append(red(message, None, self.table_title, self.table_type))
            elif header_len == 5:
                if "5pl" in regression_model.lower():
                    result.append(green(f'5PL model confirmed {analyte_str}', None, self.table_title, self.table_type))
                else:
                    message = data_flag.get_message("5pl_model_mismatch", regression_model=regression_model,
                                                    analyte_str=analyte_str)
                    result.append(red(message, None, self.table_title, self.table_type))
        else:
            result.append(
                yellow(f'Regression model not provided {analyte_str}', None, self.table_title, self.table_type))

        if len(result) == 0:
            result.append(
                green(f'All values within nominal range {analyte_str}', None, self.table_title, self.table_type))

        self.pass_batches = [x for x in self.total_batches if x not in self.failed_batches]

        self.result += result

    @staticmethod
    def get_error_dict():
        data_flag = FlagProperties("regression_data")
        return data_flag


class Specificity(Table, Flags):
    def __init__(self, parsed_table, multiple_analyte, analytes, template_type):
        Table.__init__(self, parsed_table, template_type)
        self.orientation = ""
        self.multiple_analyte = multiple_analyte
        self.analyte = get_analyte(analytes, self.tb_title)

    def validate_table_format(self):
        error_messages = self.error_messages
        self.orientation = utils.check_orientation(self.table_df)
        if self.orientation == "HZ":
            required_col = {"Run ID"}
            missing_col = required_col.intersection(self.missing_col)
            if missing_col:
                message = error_messages.get_message("missing_col", col_names=", ".join(missing_col))
                self.result.append(red(message, None, self.table_title, self.table_type))
                self.valid_format = False

        elif self.orientation == "VR":
            required_col = {"Sample"}
            missing_col = required_col.intersection(self.missing_col)
            if missing_col:
                message = error_messages.get_message("missing_col", col_names=", ".join(missing_col))
                self.result.append(red(message, None, self.table_title, self.table_type))
                self.valid_format = False

        conc_units = utils.find_units(self.table_df.columns)
        if conc_units == "":
            message = error_messages.get_message("missing_conc_unit")
            self.result.append(red(message, None, self.table_title, self.table_type))
            self.valid_format = False
        else:
            self.conc_unit = conc_units

    def get_error_dict(self):
        if self.orientation == "VR":
            data_flag = FlagProperties("vr_specificity_data")
        elif self.orientation == "HZ":
            data_flag = FlagProperties("hz_specificity_data")

        stats_flag = FlagProperties("hz_specificity_stats")
        return data_flag, stats_flag

    def process_table(self):
        table = self.table_df
        static_df = pd.DataFrame()
        final_df = pd.DataFrame()
        calc_static_df = pd.DataFrame()
        self.found_re = utils.find_re_column(table)
        count = 0
        try:
            if self.orientation == "VR":
                self.found_cv = utils.find_cv_column(table)
                conc_column = list(filter(lambda x: self.conc_unit in str(x), table.columns))
                try:
                    table["nominal"] = table[conc_column[0]].apply(utils.parse_decimal)
                    table["column"] = table[conc_column[0]]
                    table["conc"] = table[conc_column[1]].apply(utils.parse_decimal)
                    if self.found_cv:
                        table["cv"] = table["cv"].apply(utils.parse_decimal)

                    table["calc_re"] = ((table["conc"] - table["nominal"]) / table["nominal"]) * 100

                    if self.found_re:
                        table["re"] = table["re"].apply(utils.parse_decimal)
                        table["per_diff_re"] = utils.calculate_per_diff(table["re"], table["calc_re"])

                    final_df = pd.concat([final_df, table]).reset_index(drop=True)
                except Exception as e:
                    pass
            elif self.orientation == "HZ":
                nominal_conc = utils.find_nominal_conc(table.columns)
                table_df, static_df = utils.split_static_table(table)
                table_df = utils.fill_val(table_df)

                for i, key in enumerate(nominal_conc):
                    run_df = pd.DataFrame(columns=["run_id", "column", "nominal", "conc", "re"])
                    values = table_df[key].apply(utils.parse_decimal)
                    run_df["conc"] = table_df[key]
                    run_df["run_id"] = table_df["run_id"]
                    run_df["column"] = key
                    nominal = nominal_conc[key]
                    run_df["nominal"] = nominal
                    if self.found_re:
                        try:
                            column_index = table_df.columns.get_loc(key)
                            table_df.columns.values[column_index + 1] = f"re_{count}"
                            run_df["re"] = table_df[f"re_{count}"].apply(utils.parse_decimal)
                            count += 1
                        except:
                            run_df["re"] = table_df[f"re"].apply(utils.parse_decimal)

                    final_df = pd.concat([final_df, run_df]).reset_index(drop=True)
                    calc_static_df = pd.concat([calc_static_df, utils.build_static_df(values, nominal)]).reset_index(
                        drop=True)

                dummy_df = final_df.copy()
                dummy_df["conc"] = dummy_df["conc"].apply(utils.parse_decimal)
                final_df["calc_re"] = ((dummy_df["conc"] - dummy_df["nominal"]) / dummy_df["nominal"]) * 100
                if self.found_re:
                    final_df["per_diff_re"] = utils.calculate_per_diff(final_df["re"], final_df["calc_re"])

                static_df = utils.process_static_df(static_df, self.table_type)
                static_df, self.found_sd = utils.concat_static_df(static_df, calc_static_df)

            self.data_df = final_df.fillna("")
            self.static_df = static_df
        except:
            message = self.error_messages.get_message("data_error")
            self.result.append(red(message, None, self.table_title, self.table_type))
            self.valid_data = False

    def validate_table_data(self, LLOQ, ULOQ):
        valid_count_diff = Decimal(self.valid_difference_values.get_message("difference_n"))
        valid_mean_diff = Decimal(self.valid_difference_values.get_message("difference_mean"))
        valid_re_cv_diff = Decimal(self.valid_difference_values.get_message("difference_re_cv"))
        validation_failed_flag_count = 0

        if self.multiple_analyte:
            LLOQ = LLOQ.get(str(self.analyte), 0)
            ULOQ = ULOQ.get(str(self.analyte), 0)

        data_flag, stats_flag = self.get_error_dict()
        result = []
        table = self.data_df
        threshold = Decimal(self.threshold_values.get_message("re_cv_threshold"))
        if "lm" in self.analysis_type:
            threshold = threshold + Decimal(5)
        nominal_values = list(table["nominal"].unique())

        if "lm" in self.analysis_type:
            if len(nominal_values) < 2:
                message = data_flag.get_message("ul_requirement")
                result.append(red(message, None, self.table_title, self.table_type))
                validation_failed_flag_count += 1

            elif LLOQ != 0 and ULOQ != 0:
                if LLOQ not in nominal_values or ULOQ not in nominal_values:
                    message = data_flag.get_message("ul_requirement")
                    result.append(red(message, None, self.table_title, self.table_type))
                    validation_failed_flag_count += 1

        if self.orientation == "VR":
            # LLOQ = min(nominal_values)
            # ULOQ = max(nominal_values)
            samples = table["sample"]
            all_column = table["column"]
            calc_res = table["calc_re"]
            concentrations = table["conc"]
            nominal_values = table["nominal"]
            if self.found_cv:
                cvs = table["cv"]
            if self.found_re:
                reported_res = table["re"]
                per_diff_res = table["per_diff_re"]

            for index, column in enumerate(all_column):
                sample = samples[index]
                calc_re = calc_res[index]
                nominal = nominal_values[index]
                conc = concentrations[index]
                if str(LLOQ) == str(nominal) or str(ULOQ) == str(nominal):
                    test_threshold = threshold + Decimal(5)
                else:
                    test_threshold = threshold
                threshold_red = {"sample": sample, "column": column, "threshold": test_threshold, "conc": conc}
                threshold_yellow = {"sample": sample, "column": column, "threshold": test_threshold - Decimal(5),
                                    "conc": conc}
                if self.found_re:
                    reported_re = reported_res[index]
                    per_diff_re = per_diff_res[index]

                    re_para = {"sample": sample, "column": column, "reported_value": reported_re,
                               "calc_value": utils.format_value(calc_re)}
                    if str(per_diff_re) != "":
                        if per_diff_re > valid_re_cv_diff:
                            if abs(calc_re) <= Decimal(1) and abs(reported_re) <= Decimal(1):
                                message = data_flag.get_message("re_rounding", **re_para)
                                result.append(yellow(message, None, self.table_title, self.table_type))
                            else:
                                message = data_flag.get_message("re_error", **re_para)
                                result.append(red(message, None, self.table_title, self.table_type))
                                validation_failed_flag_count += 1

                if str(calc_re).strip() != "":
                    if abs(calc_re) > test_threshold:
                        message = data_flag.get_message("re", **threshold_red)
                        result.append(red(message, None, self.table_title, self.table_type))
                        validation_failed_flag_count += 1

                    elif (test_threshold - Decimal(5)) <= abs(calc_re) <= test_threshold:
                        message = data_flag.get_message("re", **threshold_yellow)
                        result.append(yellow(message, None, self.table_title, self.table_type))

                if self.found_cv:
                    cv = cvs[index]
                    if str(cv).strip() != "":
                        if abs(cv) > test_threshold:
                            message = data_flag.get_message("cv", **threshold_red)
                            result.append(red(message, None, self.table_title, self.table_type))
                            validation_failed_flag_count += 1

                        elif (test_threshold - Decimal(5)) <= abs(cv) < test_threshold:
                            message = data_flag.get_message("cv", **threshold_yellow)
                            result.append(yellow(message, None, self.table_title, self.table_type))

        elif self.orientation == "HZ":

            all_column = table["column"]
            calc_res = table["calc_re"]
            run_ids = table["run_id"]
            concentrations = table["conc"]
            if self.found_re:
                reported_res = table["re"]
                per_diff_res = table["per_diff_re"]

            for index in range(len(all_column)):
                batch = run_ids[index]
                column = all_column[index]
                calc_re = calc_res[index]
                conc = concentrations[index]

                if "lloq" in str(column).lower() or "uloq" in str(column).lower():
                    test_threshold = threshold + Decimal(5)
                else:
                    test_threshold = threshold

                threshold_red = {"batch": batch, "column": column, "threshold": test_threshold, "conc": conc}
                threshold_yellow = {"batch": batch, "column": column, "threshold": test_threshold - Decimal(5),
                                    "conc": conc}

                if self.found_re:
                    reported_re = reported_res[index]
                    per_diff_re = per_diff_res[index]

                    re_para = {"column": column, "batch": batch, "reported_value": reported_re,
                               "calc_value": utils.format_value(calc_re)}
                    if str(per_diff_re) != "":
                        if per_diff_re > valid_re_cv_diff:
                            if abs(calc_re) <= Decimal(1) and abs(reported_re) <= Decimal(1):
                                message = data_flag.get_message("re_rounding", **re_para)
                                result.append(yellow(message, None, self.table_title, self.table_type))
                            else:
                                message = data_flag.get_message("re_error", **re_para)
                                result.append(red(message, None, self.table_title, self.table_type))
                                validation_failed_flag_count += 1

                if str(calc_re).strip() != "":
                    if abs(calc_re) > test_threshold:
                        message = data_flag.get_message("re", **threshold_red)
                        result.append(red(message, None, self.table_title, self.table_type))
                        validation_failed_flag_count += 1

                    elif (test_threshold - Decimal(5)) <= abs(calc_re) <= test_threshold:
                        message = data_flag.get_message("re", **threshold_yellow)
                        result.append(yellow(message, None, self.table_title, self.table_type))

            static_df = self.static_df
            if not static_df.empty:
                all_column = static_df["column"]
                reported_mean = static_df["mean"]
                reported_cv = static_df["cv"]
                reported_re = static_df["re"]
                reported_n = static_df["n"]

                calculated_mean = static_df["calc_mean"]
                calculated_cv = static_df["calc_cv"]
                calculated_re = static_df["calc_re"]
                calculated_n = static_df["calc_n"]

                per_diff_mean = static_df["per_diff_mean"]
                per_diff_cv = static_df["per_diff_cv"]
                per_diff_re = static_df["per_diff_re"]
                per_diff_n = static_df["per_diff_n"]

                if self.found_sd:
                    reported_sd = static_df["sd"]
                    calculated_sd = static_df["calc_sd"]
                    per_diff_sd = static_df["per_diff_sd"]

                for index, column in enumerate(all_column):

                    overall_reported_mean = reported_mean[index]
                    overall_reported_cv = reported_cv[index]
                    overall_reported_re = reported_re[index]
                    overall_reported_n = reported_n[index]

                    overall_clc_mean = calculated_mean[index]
                    overall_clc_cv = calculated_cv[index]
                    overall_clc_re = calculated_re[index]
                    overall_clc_n = calculated_n[index]

                    overall_per_diff_mean = per_diff_mean[index]
                    overall_per_diff_cv = per_diff_cv[index]
                    overall_per_diff_re = per_diff_re[index]
                    overall_per_diff_n = per_diff_n[index]

                    mean_para = {"column": column, "reported_value": overall_reported_mean,
                                 "calc_value": utils.format_value(overall_clc_mean)}
                    cv_para = {"column": column, "reported_value": overall_reported_cv,
                               "calc_value": utils.format_value(overall_clc_cv)}
                    re_para = {"column": column, "reported_value": overall_reported_re,
                               "calc_value": utils.format_value(overall_clc_re)}

                    if self.found_sd:
                        overall_reported_sd = reported_sd[index]
                        overall_clc_sd = calculated_sd[index]
                        overall_per_diff_sd = per_diff_sd[index]

                        sd_para = {"column": column, "reported_value": overall_reported_sd,
                                   "calc_value": utils.format_value(overall_clc_sd)}

                        if overall_per_diff_sd > Decimal(0.15):
                            if abs(overall_clc_sd) <= Decimal(1) and abs(overall_reported_sd) <= Decimal(1):
                                message = stats_flag.get_message("sd_rounding", **sd_para)
                                result.append(yellow(message, None, self.table_title, self.table_type))

                            else:
                                message = stats_flag.get_message("sd_error", **sd_para)
                                result.append(red(message, None, self.table_title, self.table_type))
                                validation_failed_flag_count += 1

                    if "lloq" in str(column).lower() or "uloq" in str(column).lower():
                        test_threshold = threshold + Decimal(5)
                    else:
                        test_threshold = threshold

                    threshold_red = {"column": column, "threshold": test_threshold}
                    threshold_yellow = {"column": column, "threshold": test_threshold - Decimal(5)}

                    if overall_per_diff_n > valid_count_diff:
                        message = stats_flag.get_message("count_error", column=column,
                                                         reported_value=overall_reported_n, calc_value=overall_clc_n)
                        result.append(red(message, None, self.table_title, self.table_type))
                        validation_failed_flag_count += 1

                    if overall_per_diff_mean > valid_mean_diff:
                        if abs(overall_clc_mean) <= Decimal(1) and abs(overall_reported_mean) <= Decimal(1):
                            message = stats_flag.get_message("mean_rounding", **mean_para)
                            result.append(yellow(message, None, self.table_title, self.table_type))

                        else:
                            message = stats_flag.get_message("mean_error", **mean_para)
                            result.append(red(message, None, self.table_title, self.table_type))
                            validation_failed_flag_count += 1

                    if overall_per_diff_cv > valid_re_cv_diff:
                        if abs(overall_reported_cv) <= Decimal(1) and abs(overall_clc_cv) <= Decimal(1):
                            message = stats_flag.get_message("cv_rounding", **cv_para)
                            result.append(yellow(message, None, self.table_title, self.table_type))

                        else:
                            message = stats_flag.get_message("cv_error", **cv_para)
                            result.append(red(message, None, self.table_title, self.table_type))
                            validation_failed_flag_count += 1

                    if overall_per_diff_re > valid_re_cv_diff:
                        if abs(overall_reported_re) <= Decimal(1) and abs(overall_clc_re) <= Decimal(1):
                            message = stats_flag.get_message("re_rounding", **re_para)
                            result.append(yellow(message, None, self.table_title, self.table_type))

                        else:
                            message = stats_flag.get_message("re_error", **re_para)
                            result.append(red(message, None, self.table_title, self.table_type))
                            validation_failed_flag_count += 1

                    if abs(overall_clc_re) > test_threshold:
                        message = stats_flag.get_message("re", **threshold_red)
                        result.append(red(message, None, self.table_title, self.table_type))
                        validation_failed_flag_count += 1

                    elif (test_threshold - Decimal(5)) <= abs(overall_clc_re) <= test_threshold:
                        message = stats_flag.get_message("re", **threshold_yellow)
                        result.append(yellow(message, None, self.table_title, self.table_type))

                    if abs(overall_clc_cv) > test_threshold:
                        message = stats_flag.get_message("cv", **threshold_red)
                        result.append(red(message, None, self.table_title, self.table_type))
                        validation_failed_flag_count += 1

                    elif (test_threshold - Decimal(5)) <= abs(overall_clc_cv) <= test_threshold:
                        message = stats_flag.get_message("cv", **threshold_yellow)
                        result.append(yellow(message, None, self.table_title, self.table_type))

        if validation_failed_flag_count == 0:
            result.append(green('All values within acceptable range', None, self.table_title, self.table_type))

        self.result += result


class BatchPerformance(Table, Flags):
    def __init__(self, parsed_table, options, expiration_date, analytes, multiple_analyte, template_type):
        Table.__init__(self, parsed_table, template_type)
        self.analysts = list()
        self.verified_batch_table = pd.DataFrame()
        self.options = options
        self.expiration_date = expiration_date
        self.pass_failed_batches = dict()
        self.instrument = True
        self.HPLC_SYSTEM = False
        self.verify_ref_date = True
        self.analytes = analytes
        self.num_analyte = len(self.analytes)
        self.multiple_analyte = multiple_analyte
        self.consecutive_failed = False
        self.missing_analyst = False

    def validate_table_format(self):
        error_messages = self.error_messages
        required_col = {"Run ID", "Batch Acceptance"}
        missing_col = required_col.intersection(self.missing_col)
        if missing_col:
            message = error_messages.get_message("missing_col", col_names=", ".join(missing_col))
            self.result.append(red(message, None, self.table_title, self.table_type))
            self.valid_format = False

        required_col = {"Analyst"}
        missing_col = required_col.intersection(self.missing_col)
        if missing_col:
            # message = error_messages.get_message("missing_col", col_names=", ".join(missing_col))
            message = "Identification of the analyst responsible for a given batch run was not found in the table"
            if "lm" in self.analysis_type and "mv" in self.analysis_type:
                self.result.append(red(message, None, self.table_title, self.table_type))
                self.missing_analyst = True
            else:
                self.result.append(yellow(message, None, self.table_title, self.table_type))
                self.missing_analyst = True
        if self.num_analyte > 1:
            required_col = {"Analyte"}
            missing_col = required_col.intersection(self.missing_col)
            if missing_col:
                message = error_messages.get_message("missing_col", col_names=", ".join(missing_col))
                self.result.append(red(message, None, self.table_title, self.table_type))
                self.valid_format = False

        elif "analyte" in self.table_df.columns and self.num_analyte == 1:
            self.result.append(yellow(
                f"Single analyte name entered by the user; uploaded file found to contain results for multiple analytes"
                f" and may result in an incomplete evaluation by Red Thread™", None, self.table_title, self.table_type))
            self.num_analyte = 2

        if "Instrument ID" in self.missing_col:
            self.instrument = False
            self.HPLC_SYSTEM = True

        if not self.instrument:
            if len({"MS System #", "HPLC System #"}.intersection(set(self.missing_col))) > 0:
                self.instrument = False
            else:
                self.instrument = True

        if not self.instrument:
            self.result.append(
                yellow("Missing column for Instrument ID in Batch Run Summary Table", None, self.table_title,
                       self.table_type))

    @staticmethod
    def get_error_dict():
        data_flag = FlagProperties("BatchPerformance")
        return data_flag

    @staticmethod
    def is_batch_passed(run_status):
        if ("pass" in run_status or "accept" in run_status or "yes" in run_status or (
                "met" in run_status and "crit" in run_status)) and (
                "unaccept" not in run_status and "not" not in run_status):
            return True
        else:
            return False

    def process_table(self):
        table = self.table_df
        final_df = pd.DataFrame()
        self.data_flag = self.get_error_dict()

        if len(self.analytes) > 1 or self.multiple_analyte:
            table["split"] = table["analyte"].apply(check_split_batch_analyte)
            dummy_table = table[table["split"]]
            analytes = dummy_table["analyte"]
            regression_status = dummy_table["regression_status"]
            dummy_table = dummy_table.drop(["analyte", "regression_status"], axis=1)
            split_analytes = []
            split_status = []
            for idx, analyte in enumerate(analytes):
                df = dummy_table[idx: idx + 1]
                analyte = split_batch_analyte(analyte)
                status = split_batch_analyte(regression_status[idx])
                split_analytes += analyte
                analyte_num = len(analyte)
                status_num = len(status)
                if analyte_num == status_num:
                    split_status += status
                elif status_num == 1:
                    split_status += status * analyte_num

                final_df = pd.concat([final_df, *[df] * analyte_num]).reset_index(drop=True)
            final_df["analyte"] = split_analytes
            final_df["regression_status"] = split_status
            final_df = pd.concat([final_df, table[~table["split"]]]).reset_index(drop=True)
        else:
            final_df = table

        if not self.missing_analyst:
            try:
                analyst_names = list(final_df["analyst"].unique())
                for analyst_name in analyst_names:
                    status = final_df[final_df["analyst"] == analyst_name]["regression_status"].to_list()
                    total = len(status)
                    pass_batch = list(filter(self.is_batch_passed, status))
                    failed_batch = total - len(pass_batch)
                    self.analysts.append(
                        {'analyst': analyst_name, 'failed': failed_batch, 'total': total, "pass": pass_batch})

            except KeyError:
                message = self.error_messages.get_message("data_error")
                self.result.append(red(message, None, self.table_title, self.table_type))
                self.valid_data = False

        self.data_df = final_df

    def get_analyst_batch(self):
        analysts = list()
        table = self.verified_batch_table

        if not self.missing_analyst:
            try:
                analyst_names = list(table["analyst"].unique())
                for analyst_name in analyst_names:
                    failed_batch = 0
                    total = 0
                    pass_batch = 0

                    analyst_table = table[table["analyst"] == analyst_name]
                    if self.multiple_analyte:
                        total = analyst_table.shape[0]
                        analytes = analyst_table["analyte"].unique()
                        analytes_batch_str = ""
                        num_analytes = len(analytes)
                        for idx, analyte in enumerate(analytes):
                            status = analyst_table[analyst_table["analyte"] == analyte]["regression_status"].to_list()
                            analyte_batch = len(status)
                            analyte_pass_batch = list(filter(self.is_batch_passed, status))
                            analyte_failed_batch = analyte_batch - len(analyte_pass_batch)
                            pass_batch += len(analyte_pass_batch)
                            failed_batch += analyte_failed_batch

                            analytes_batch_str += f'{analyte_failed_batch} for {analyte}, ' if num_analytes - 1 != idx else f'and {analyte_failed_batch} for {analyte}'

                        analysts.append(
                            {'analyst': analyst_name, 'failed': failed_batch, 'total': total, "pass": pass_batch,
                             "analytes_batch_str": analytes_batch_str})

                    else:
                        status = table[table["analyst"] == analyst_name]["regression_status"].to_list()
                        total = len(status)
                        analyst_pass_batch = list(filter(self.is_batch_passed, status))
                        failed_batch = total - len(analyst_pass_batch)
                        pass_batch = len(analyst_pass_batch)

                        analysts.append(
                            {'analyst': analyst_name, 'failed': failed_batch, 'total': total, "pass": pass_batch})

            except KeyError:
                message = self.error_messages.get_message("data_error")
                self.result.append(red(message, None, self.table_title, self.table_type))
                self.valid_data = False
        return analysts

    def check_ref_date_and_instrument(self):
        result = []
        data_flag = self.data_flag
        table = self.data_df
        table_title = self.table_title
        table_type = self.table_type
        ref_date = utils.format_date(self.options.get("ref_date", None))

        def get_recent_date(assay_date):
            try:
                hours_diff = (datetime.date.today() - utils.parse_date(assay_date)).days * 24
                return hours_diff
            except:
                return False

        try:
            table["assay_date_diff"] = np.vectorize(get_recent_date)(table["assay_date"])
            assay_date = table[table["assay_date_diff"] == table["assay_date_diff"].min()]["assay_date"].to_list()[0]
            if ref_date is not None:
                if (utils.parse_date(ref_date) - utils.parse_date(assay_date)).days > 0:
                    self.verify_ref_date = False
                    message = data_flag.get_message("ref_date_entry_error")
                    result.append(yellow(message, None, table_title, table_type))
        except:
            pass
        if self.instrument:
            batches = table["run_id"].to_list()
            if not self.HPLC_SYSTEM:
                instrument_ids = table["instrument_id"].to_list()
                for index, idx in enumerate(instrument_ids):
                    if str(idx).strip() == "":
                        message = data_flag.get_message("missing_instrument", batch=batches[index])
                        result.append(yellow(message, None, table_title, table_type))
            else:
                hplc_ids = table["hplc_id"].to_list()
                ms_ids = table["ms_id"].to_list()
                for index, _id in enumerate(hplc_ids):
                    if str(_id).strip() == "" or str(ms_ids[index]).strip() == "":
                        message = data_flag.get_message("missing_instrument", batch=batches[index])
                        result.append(yellow(message, None, table_title, table_type))

        return result

    def batch_summary_check_batch_numbers(self):
        result = []
        table_title = self.table_title
        table_type = self.table_type
        table_rows = self.data_df

        ok_result = green('Batch Number Check Ok, No Missing batch found', None, table_title, table_type)
        batch_numbers = list(table_rows['run_id'])

        try:
            num_analyte = len(list(set(list(table_rows['analyte']))))
        except KeyError:
            num_analyte = 1

        if not isinstance(batch_numbers[0], int):
            batch_numbers = utils.get_unique_ids(batch_numbers)

        if len(batch_numbers) == 0:
            return [red('No batch numbers found in table', None, table_title, table_type)]

        batch_numbers = [x for x in batch_numbers if str(x) != ""]

        sorted_batch_numbers = sorted(batch_numbers)

        if (list(sorted_batch_numbers) != list(batch_numbers)) and num_analyte == 1:
            result.append(yellow('Batch numbers out of order', None, table_title, table_type))

        first_batch_number = sorted_batch_numbers[0]
        num_batches = len(sorted_batch_numbers)

        expected_batch_numbers = range(1, first_batch_number + num_batches)

        if list(expected_batch_numbers) != list(sorted_batch_numbers):
            dupes = [item for item, count in collections.Counter(sorted_batch_numbers).items() if count > 1]
            if len(dupes) > 0 and num_analyte == 1:
                result.append(
                    yellow('Duplicate batch numbers: %s' % ', '.join(map(str, dupes)), None, table_title, table_type))

            possible_batches = set(range(1, sorted_batch_numbers[-1] + 1))
            missing = possible_batches.difference(set(sorted_batch_numbers))
            if len(missing) > 0:
                result.append(
                    yellow('Missing batch numbers: %s' % ', '.join(map(str, missing)), None, table_title, table_type))

        if len(result) == 0:
            return [ok_result]
        else:
            return result

    def batch_summary_check_extraction_stability(self):
        extraction_stability = self.options.get('extraction_stability', None)
        data_flag = self.data_flag
        table_title = self.table_title
        table_type = self.table_type
        table_rows = self.data_df

        result = []

        if extraction_stability is None:
            result.append(yellow('No Extraction Stability Provided', None, table_title, table_type))

        try:
            if table_rows['assay_date'].empty or table_rows['assay_date'].dropna().empty:
                result.append(yellow('No Assay Date Found', None, table_title, table_type))
        except:
            result.append(yellow('No Assay Date Found', None, table_title, table_type))

        try:
            if table_rows['extraction_date'].empty or table_rows['extraction_date'].dropna().empty:
                result.append(yellow('No Extraction Date Found', None, table_title, table_type))
        except:
            result.append(yellow('No Extraction Date Found', None, table_title, table_type))

        if len(result) > 0:
            return result

        def test_stability(extraction_date, assay_date, stability):
            try:
                hours_diff = (utils.parse_date(assay_date) - utils.parse_date(extraction_date)).days * 24
                if hours_diff > stability:
                    return True
                else:
                    return False
            except:
                return False

        def test_extraction(extraction_date, assay_date):
            try:
                hours_diff = (utils.parse_date(extraction_date) - utils.parse_date(assay_date)).days * 24
                if hours_diff > 0:
                    return True
                else:
                    return False
            except:
                return False

        missing_batches = table_rows[table_rows['extraction_date'].replace('', np.nan).isnull()]["run_id"].tolist()

        expired_batches = table_rows[
            table_rows.apply(lambda x: test_stability(x.extraction_date, x.assay_date, extraction_stability), axis=1)][
            "run_id"].tolist()

        extraction_batches = \
            table_rows[table_rows.apply(lambda x: test_extraction(x.extraction_date, x.assay_date), axis=1)][
                "run_id"].tolist()

        for batch in missing_batches:
            result.append(yellow(f'[Batch {batch}] missing extraction date', None, table_title, table_type))

        for batch in expired_batches:
            message = data_flag.get_message("stability_period_exceeded", batch=batch,
                                            extraction_stability=extraction_stability)
            result.append(red(message, None, table_title, table_type))

        for batch in extraction_batches:
            message = data_flag.get_message("extraction_after_analysis", batch=batch)
            result.append(red(message, None, table_title, table_type))
        if len(result) == 0:
            result.append(
                green('All batches assayed within stability period of ' + str(extraction_stability) + " hours", None,
                      table_title, table_type))

        return result

    def stability_date_check(self):
        ref_date = self.options.get('ref_date', None) or None
        lts_days = self.options.get('lts80', None) or None
        data_flag = self.data_flag
        table_title = self.table_title
        table_type = self.table_type
        table_rows = self.data_df

        result = []

        if ref_date is None or str(ref_date).strip() == "":
            result.append(yellow('No Reference Date Provided', None, table_title, table_type))

        if lts_days is None or str(lts_days).strip() == "":
            result.append(yellow('No Long-Term Stability Days Provided', None, table_title, table_type))

        try:
            if table_rows['assay_date'].empty or table_rows['assay_date'].dropna().empty:
                result.append(yellow('No Assay Date Found', None, table_title, table_type))
        except:
            pass

        if len(result) > 0:
            return result

        def test_stability(r_date, assay_date, lts_day):
            r_date = utils.parse_date(r_date)
            assay_date = utils.parse_date(assay_date)
            lts_day = int(lts_day)
            try:
                diff = (assay_date - r_date).days
                if diff > lts_day or diff < 0:
                    return True
                else:
                    return False
            except:
                return False

        if self.verify_ref_date:
            expired_batches = \
                table_rows[table_rows.apply(lambda x: test_stability(ref_date, x.assay_date, lts_days), axis=1)][
                    "run_id"].tolist()

            if len(expired_batches) > 0:
                message = data_flag.get_message("expired_batch")
                result.append(red(message, None, table_title, table_type))

            if len(result) == 0:
                result.append(
                    green(f'All batches assayed within established long-term stability', None, table_title, table_type))
        else:
            message = data_flag.get_message("unable_verify_batch_in_tp")
            result.append(red(message, None, table_title, table_type))
        return result

    def batch_summary_check_expiry(self):
        table_title = self.table_title
        table_type = self.table_type
        table_rows = self.data_df
        data_flag = self.data_flag
        result = []

        if table_rows.empty:
            return []

        try:
            if table_rows['assay_date'].empty or table_rows['assay_date'].dropna().empty:
                result.append(yellow('No Assay Date Found', None, table_title, table_type))
        except KeyError:
            result.append(yellow('No Assay Date Found', None, table_title, table_type))

        missing_batches = table_rows[table_rows['assay_date'].replace("", np.nan).isnull()]["run_id"].tolist()

        for batch in missing_batches:
            message = data_flag.get_message("missing_assay_date", batch=batch)
            result.append(yellow(message, None, table_title, table_type))

        if self.expiration_date is None:
            result.append(yellow('Could not find expiry dates', None, 'Reference Standards', 'Reference Standards'))
            return result

        def test_expiration(exp_date, run_date):
            try:
                run_date = utils.parse_date(run_date)
                diff = (exp_date - run_date).days
                if diff < 0:
                    return True
                else:
                    return False
            except:
                return False

        expired_batches = \
            table_rows[table_rows['assay_date'].apply(lambda x: test_expiration(self.expiration_date, x))][
                "run_id"].tolist()

        for batch in expired_batches:
            message = data_flag.get_message("after_expiry", batch=batch,
                                            expiration_date=utils.format_date(self.expiration_date))
            result.append(red(message, None, table_title, table_type))

        if len(result) == 0:
            message = data_flag.get_message("before_expiry", expiration_date=utils.format_date(self.expiration_date))
            result.append(green(message, None, table_title, table_type))

        return result

    def batch_summary_check_accept_reject_ratio(self, total_batches, pass_failed_batches):
        results = []
        extra_batch = []
        total = []
        verified_batch_table = pd.DataFrame()
        verified_batch_status = []
        data_flag = self.data_flag
        table_title = self.table_title
        table_type = self.table_type
        self.pass_failed_batches = pass_failed_batches
        num_analyte = self.num_analyte
        for analyte in self.analytes:
            if self.multiple_analyte:
                bc_regression_batches = list(set(total_batches["bc_batches"][analyte]).intersection(
                    set(total_batches["regression_batches"][analyte])))
                analyte_str = f"for analyte {analyte}"
                pass_batches = pass_failed_batches["bc_pass_batches"][analyte] + \
                               pass_failed_batches["regression_pass_batches"][analyte]
                fail_batches = pass_failed_batches["bc_failed_batches"][analyte] + \
                               pass_failed_batches["regression_failed_batches"][analyte]

                qc_batches = total_batches["qc_batches"][analyte]
                qc_pass_batches = pass_failed_batches["qc_pass_batches"][analyte]
                qc_failed_batches = pass_failed_batches["qc_failed_batches"][analyte]

                table_rows = self.data_df[self.data_df["analyte"].str.lower() == analyte.lower()]
            else:
                bc_regression_batches = list(
                    set(total_batches["bc_batches"]).intersection(set(total_batches["regression_batches"])))
                analyte_str = ""
                pass_batches = pass_failed_batches["bc_pass_batches"] + pass_failed_batches["regression_pass_batches"]
                fail_batches = pass_failed_batches["bc_failed_batches"] + pass_failed_batches[
                    "regression_failed_batches"]
                table_rows = self.data_df
                qc_batches = total_batches["qc_batches"]
                qc_pass_batches = pass_failed_batches["qc_pass_batches"]
                qc_failed_batches = pass_failed_batches["qc_failed_batches"]
            if num_analyte > 1 and table_rows.empty:
                results.append(
                    red(f"Could not find analyte {analyte} in the Batch Performance table, unable to verify batch accepted/rejected status",
                        None, table_title, table_type))
                continue

            table_rows["run_id"] = table_rows["run_id"].apply(utils.parse_int)
            table_rows = table_rows.sort_values("run_id").reset_index(drop=True)
            table_rows["run_id"] = table_rows["run_id"].astype(str)
            reported_batches = list(table_rows["run_id"].unique())
            if "mv" in self.analysis_type:
                if self.multiple_analyte:
                    ap_batches = total_batches["ap_batches"][analyte]
                    ap_pass_batches = pass_failed_batches["ap_pass_batches"][analyte]
                    ap_failed_batches = pass_failed_batches["ap_failed_batches"][analyte]
                else:
                    ap_batches = total_batches["ap_batches"]
                    ap_pass_batches = pass_failed_batches["ap_pass_batches"]
                    ap_failed_batches = pass_failed_batches["ap_failed_batches"]

                ap_batches = sorted([utils.parse_int(x) for x in ap_batches])
                ap_pass_batches = sorted([utils.parse_int(x) for x in ap_pass_batches])
                ap_failed_batches = sorted([utils.parse_int(x) for x in ap_failed_batches])
                add_str = "was omitted from the calibration standards table or regression model/curve parameter table"
            else:
                add_str = "was omitted from the calibration standards table or regression model/curve parameter table or QC samples table"

            qc_batches = sorted([utils.parse_int(x) for x in qc_batches])
            qc_pass_batches = sorted([utils.parse_int(x) for x in qc_pass_batches])
            qc_failed_batches = sorted([utils.parse_int(x) for x in qc_failed_batches])

            pass_batches.extend(qc_pass_batches)
            fail_batches.extend(qc_failed_batches)

            bc_regression_batches = sorted([utils.parse_int(x) for x in bc_regression_batches])

            fail_batches = sorted(set([utils.parse_int(x) for x in fail_batches]))
            pass_batches = sorted(set([utils.parse_int(x) for x in pass_batches]))
            pass_batches = [x for x in pass_batches if x not in fail_batches]

            reported_batches = sorted([x if str(x) == "NA" else utils.parse_int(x) for x in reported_batches])

            if len(bc_regression_batches) > 0:
                extra_batch = sorted([x for x in reported_batches if x not in bc_regression_batches])
                batch_iter = iter(extra_batch)
                next_batch = next(batch_iter, 0)
                next_batch = next(batch_iter, 0)
                if len(extra_batch) > 0:
                    for batch in extra_batch:
                        batch_status = table_rows[table_rows["run_id"] == str(batch)]["regression_status"].to_list()
                        for run_status in batch_status:
                            status = self.check_status(run_status)
                            if status.strip() != "":
                                message = data_flag.get_message("unable_to_verify_batch", batch=batch, status=status,
                                                                add_str=add_str, analyte_str=analyte_str)
                                results.append(red(message, None, table_title, table_type))
                                if next_batch == batch + 1 and status == 'rejected':
                                    self.consecutive_failed = True
                                next_batch = next(batch_iter, 0)
                        table_rows = table_rows[table_rows["run_id"] != str(batch)].reset_index(drop=True)
            passed = 0
            failed = 0

            verified_batch_table = pd.concat([verified_batch_table, table_rows]).reset_index(drop=True)

            if not table_rows.empty:
                regression_statuses = table_rows['regression_status'].apply(lambda x: str(x).lower()).tolist()

                batches = table_rows["run_id"].to_list()
                batches = [parse_int(x) for x in batches]

                for i, run_status in enumerate(regression_statuses):
                    status = self.check_status(run_status)
                    added_str = ""
                    batch = batches[i]
                    try:
                        if "sa" in self.analysis_type:
                            if batch in qc_batches:

                                if status == "accepted":
                                    if batch in fail_batches:
                                        failed += 1
                                        message = data_flag.get_message("with_qc_rejected", batch=batch,
                                                                        analyte_str=analyte_str)
                                        results.append(red(message, None, table_title, table_type))
                                        verified_batch_status.append("failed")
                                    else:
                                        passed = passed + 1
                                        verified_batch_status.append("passed")

                                elif status == "rejected":
                                    if batch in pass_batches:
                                        passed += 1
                                        message = data_flag.get_message("with_qc_accepted", batch=batch,
                                                                        analyte_str=analyte_str)
                                        results.append(red(message, None, table_title, table_type))
                                        verified_batch_status.append("passed")
                                    else:
                                        failed = failed + 1
                                        verified_batch_status.append("failed")
                            else:
                                message = data_flag.get_message("omitted_batch", batch=batch, status=status,
                                                                analyte_str=analyte_str)
                                verified_batch_status.append(np.nan)
                                results.append(red(message, None, table_title, table_type))
                        elif "mv" in self.analysis_type:
                            if batch in ap_batches:
                                if status == "accepted":
                                    if batch in fail_batches:
                                        failed += 1
                                        message = data_flag.get_message("without_qc_rejected", batch=batch,
                                                                        analyte_str=analyte_str)
                                        results.append(red(message, None, table_title, table_type))
                                        verified_batch_status.append("failed")
                                    else:
                                        passed = passed + 1
                                        verified_batch_status.append("passes")

                                elif status == "rejected":
                                    if batch in pass_batches:
                                        passed += 1
                                        message = data_flag.get_message("without_qc_accepted", batch=batch,
                                                                        analyte_str=analyte_str)
                                        results.append(red(message, None, table_title, table_type))
                                        verified_batch_status.append("passed")
                                    else:
                                        failed = failed + 1
                                        verified_batch_status.append("failed")

                            elif batch in qc_batches:
                                if status == "accepted":
                                    if batch in fail_batches:
                                        failed += 1
                                        message = data_flag.get_message("with_qc_rejected", batch=batch,
                                                                        analyte_str=analyte_str)
                                        results.append(red(message, None, table_title, table_type))
                                        verified_batch_status.append("failed")
                                    else:
                                        passed = passed + 1
                                        verified_batch_status.append("passed")

                                elif status == "rejected":
                                    if batch in pass_batches:
                                        passed += 1
                                        message = data_flag.get_message("with_qc_accepted", batch=batch,
                                                                        analyte_str=analyte_str)
                                        results.append(red(message, None, table_title, table_type))
                                        verified_batch_status.append("passed")
                                    else:
                                        failed = failed + 1
                                        verified_batch_status.append("failed")
                            else:
                                message = data_flag.get_message("omitted_batch", batch=batch, status=status,
                                                                analyte_str=analyte_str)
                                results.append(red(message, None, table_title, table_type))
                                verified_batch_status.append(np.nan)
                    except:
                        continue

                total = passed + failed

                if total > 0:
                    accept_ratio = passed / total

                    if accept_ratio > 0.90:
                        results.append(green(
                            f'{round(accept_ratio * 100, 1)}% batch acceptance status confirmed for verifiable batches {analyte_str}',
                            None, table_title, table_type))
                    elif accept_ratio > 0.50:
                        results.append(yellow(
                            f'{round(accept_ratio * 100, 1)}% batch acceptance status confirmed for verifiable batches {analyte_str}',
                            None, table_title, table_type))
                    else:
                        results.append(
                            red(f'{round(accept_ratio * 100, 1)}% batch acceptance status confirmed for verifiable batches {analyte_str}',
                                None, table_title, table_type))
                else:
                    results.append(
                        red(f'Could not detect regression statuses for batches {analyte_str}', None, table_title,
                            table_type))

        try:
            verified_batch_table["regression_status"] = verified_batch_status
            verified_batch_table = verified_batch_table.dropna()
        except Exception as e:
            print(e)
            pass

        self.verified_batch_table = verified_batch_table
        return results

    @staticmethod
    def check_status(run_status):
        run_status = str(run_status).lower()
        status = ""
        if ("pass" in run_status or "accept" in run_status or "yes" in run_status or (
                "met" in run_status and "crit" in run_status)) and (
                "unaccept" not in run_status and "not" not in run_status):
            status = "accepted"
        elif "fail" in run_status or "reject" in run_status or "unaccept" in run_status or (
                "no" in run_status and "na" not in run_status and "n/a" not in run_status) or "not" in run_status:
            status = "rejected"
        return status

    def check_in_failed_batches(self, batch):
        batch = str(batch)
        added_str = ""
        failed_batch = 0
        pass_failed_batch = self.pass_failed_batches
        if batch in pass_failed_batch["bc_failed_batch"]:
            added_str += "calibration curve"
        if batch in pass_failed_batch["regression_failed_batch"]:
            if len(added_str) > 1:
                added_str += " and regression model/curve parameter"

    def batch_summary_check_accept_reject_sequential(self, pass_failed_batches):
        overall_result = []
        table_title = self.table_title
        table_type = self.table_type
        num_analyte = self.num_analyte
        for analyte in self.analytes:
            result = []
            if self.multiple_analyte:
                analyte_str = f"for analyte {analyte}"
                table_rows = self.data_df[self.data_df["analyte"].str.lower() == analyte.lower()]
                fail_batches = pass_failed_batches["bc_failed_batches"][analyte] + \
                               pass_failed_batches["regression_failed_batches"][analyte] + \
                               pass_failed_batches["qc_failed_batches"][analyte]
            else:
                analyte_str = ""
                table_rows = self.data_df
                fail_batches = pass_failed_batches["bc_failed_batches"] + pass_failed_batches[
                    "regression_failed_batches"] + pass_failed_batches["qc_failed_batches"]

            ok_result = green(f'Rejection Batch Check Completed, No consecutive batches failed {analyte_str}', None,
                              table_title, table_type)

            def has_failed(x):
                x = str(x).lower()
                if ("f" in x or "rej" in x or "unaccepted" in x) or (
                        "n" in x and "na" not in x and "n/a" not in x) or "not" in x:
                    return True
                else:
                    return False

            batch_numbers = table_rows['run_id'].tolist()
            regression_statuses = table_rows['regression_status'].apply(lambda x: has_failed(x)).tolist()

            if len(fail_batches) > 0:
                fail_batches = sorted(set([utils.parse_int(x) for x in fail_batches]))
                for idx, batch in enumerate(funcy.butlast(fail_batches)):
                    next_batch = fail_batches[idx + 1]
                    if next_batch == (batch + 1):
                        result.append(
                            red(f'Batch {str(batch)} and Batch {str(next_batch)} failed consecutively {analyte_str}',
                                None, table_title, table_type))
            else:

                for idx, batch_status in enumerate(funcy.butlast(regression_statuses)):
                    my_status = batch_status
                    next_status = regression_statuses[idx + 1]
                    if my_status and next_status:
                        result.append(
                            red(f'Batch {str(batch_numbers[idx])} and Batch {str(batch_numbers[idx + 1])} failed consecutively {analyte_str}',
                                None, table_title, table_type))

            if len(result) == 0 and not self.consecutive_failed:
                result.append(ok_result)
            overall_result += result
        return overall_result

    def check_validated_run_length(self, run_length):
        table_title = self.table_title
        table_type = self.table_type
        table_rows = self.data_df
        results = []

        if table_rows.empty:
            return results

        def has_samples(raw_string):
            try:
                check_string = str(raw_string).lower().strip()
                string_array = check_string.split(" ")

                samples = 0
                samples_index = -1
                for index, item in enumerate(string_array):
                    if "sample" in item or "inject" in item:
                        samples_index = index - 1
                        if samples == 0 and samples_index != -1:
                            samples = utils.parse_signed_int(string_array[samples_index])
                            if samples is None:
                                samples = 0
                if samples == 0:
                    return False
                else:
                    return True
            except:
                return False

        def find_num_samples(raw_string):
            try:
                check_string = str(raw_string).lower().strip()
                string_array = check_string.split(" ")

                samples = 0
                samples_index = -1
                for index, item in enumerate(string_array):
                    if "sample" in item or "inject" in item:
                        samples_index = index - 1
                        if samples == 0 and samples_index != -1:
                            samples = utils.parse_signed_int(string_array[samples_index])
                            if samples is None:
                                samples = 0
                            if samples < 1:
                                samples = 0
                if samples == 0:
                    return False
                else:
                    return samples
            except:
                return False

        try:
            table_rows['has_samples'] = table_rows['description'].apply(lambda x: has_samples(x))

            has_samples_df = table_rows[table_rows['has_samples']]
            if has_samples_df.empty:
                return results
            has_samples_df['found_samples'] = has_samples_df['description'].apply(lambda x: find_num_samples(x))
            has_samples_df['failed'] = pd.DataFrame([abs(has_samples_df['found_samples']) > run_length]).transpose()

            failed = has_samples_df[has_samples_df['failed']]
            if failed.empty:
                return results

            failed_list = failed["run_id"].tolist()

            for batch_number in failed_list:
                message = str(
                    "Batch " + str(batch_number).strip() + " has more samples than the validated run length of " + str(
                        run_length).strip() + " samples").strip()
                results.append(red(message, None, table_title, table_type))
        except:
            pass

        return results


class MatrixEffect(Table, Flags):
    def __init__(self, parsed_table, analytes, template_type):
        Table.__init__(self, parsed_table, template_type)
        self.required_sample = Decimal(6)
        self.analyte = get_analyte(analytes, self.tb_title)

    def validate_table_format(self):
        error_messages = self.error_messages
        required_col = {"QC Level", "Run ID", "Matrix Sample Peak Area Ratio", "Solution Sample Mean Peak Area Ratio",
                        "Matrix Factor Individual", "Matrix Factor Mean", "Sample", "Matrix Factor %CV"}

        missing_col = required_col.intersection(self.missing_col)
        if missing_col:
            message = error_messages.get_message("missing_col", col_names=", ".join(missing_col))
            self.result.append(red(message, None, self.table_title, self.table_type))
            self.valid_format = False

        elif ["Matrix Factor SD"] in self.missing_col:
            self.result.append(
                yellow(f"Could not identify in Table: {' ,'.join(self.missing_col)}", None, self.table_title,
                       self.table_type))

    def process_table(self):
        calc_static_df = pd.DataFrame()
        static_df = pd.DataFrame()
        try:
            table = utils.fill_val(self.table_df)
            table["matrix_peak_area"] = table["matrix_peak_area"].apply(utils.parse_decimal)
            table["mean_solution_peak_area"] = table["mean_solution_peak_area"].apply(utils.parse_decimal)
            table["matrix_factor_individual"] = table["matrix_factor_individual"].apply(utils.parse_decimal)
            table["calc_matrix_factor"] = table["matrix_peak_area"] / table["mean_solution_peak_area"]

            table["per_diff_factor"] = abs(
                (table["matrix_factor_individual"] - table["calc_matrix_factor"]) / table["matrix_factor_individual"])

            qc_level = table["qc_level"].unique()
            for level in qc_level:
                values = table[table["qc_level"] == level]["calc_matrix_factor"]
                calc_static_df = pd.concat([calc_static_df, utils.build_static_df(values, 0)]).reset_index(drop=True)
            calc_static_df["qc_level"] = qc_level

            run_id = list(table["run_id"].unique())
            run_id = run_id * len(qc_level)
            try:
                mean = table["matrix_factor_mean"]
                mean = [x for x in mean if str(x).strip() != ""]
            except:
                mean = [Decimal(0)] * len(qc_level)
            try:
                cv = table["matrix_factor_cv"]
                cv = [x for x in cv if str(x).strip() != ""]
            except:
                cv = [Decimal(0)] * len(qc_level)
            try:
                sd = table["matrix_factor_sd"]
                sd = [x for x in sd if str(x).strip() != ""]
            except:
                sd = [Decimal(0)] * len(qc_level)

            res = [Decimal(0)] * len(qc_level)
            n = [Decimal(0)] * len(qc_level)
            static_df["mean"] = mean
            static_df["sd"] = sd
            static_df["cv"] = cv
            static_df["re"] = res
            static_df["n"] = n
            static_df["run_id"] = run_id

            static_df, self.found_sd = utils.concat_static_df(static_df, calc_static_df)

            self.data_df = table
            self.static_df = static_df
        except:
            message = self.error_messages.get_message("data_error")
            self.result.append(red(message, None, self.table_title, self.table_type))
            self.valid_data = False

    def validate_table_data(self):
        validation_failed_flag_count = 0
        valid_count_diff = Decimal(self.valid_difference_values.get_message("difference_n"))
        valid_mean_diff = Decimal(self.valid_difference_values.get_message("difference_mean"))
        valid_re_cv_diff = Decimal(self.valid_difference_values.get_message("difference_re_cv"))

        data_flag, stats_flag = self.get_error_dict()
        result = []
        table = self.data_df
        threshold = Decimal(self.threshold_values.get_message("re_cv_threshold"))
        if "lm" in self.analysis_type:
            threshold += Decimal(5)

        batches = table["run_id"]
        reported_factors = table["matrix_factor_individual"]
        calc_factors = table["calc_matrix_factor"]
        per_diff_factors = table["per_diff_factor"]
        levels = table["qc_level"]
        samples = table["sample"]

        for index, batch in enumerate(batches):
            reported_factor = reported_factors[index]
            calc_factor = calc_factors[index]
            per_diff_factor = per_diff_factors[index]
            level = levels[index]
            sample = samples[index]

            if per_diff_factor > Decimal(0.15):
                message = data_flag.get_message("factor_error", sample=sample, batch=batch, level=level,
                                                reported_value=reported_factor,
                                                calc_value=utils.format_value(calc_factor, 3))
                result.append(red(message, None, self.table_title, self.table_type))

        static_df = self.static_df
        if static_df is not None:
            mean = static_df["mean"].to_list()
            calc_cvs = static_df["calc_cv"].to_list()
            levels = static_df["qc_level"]
            batches = static_df["run_id"]

            reported_means = static_df["mean"]
            reported_cvs = static_df["cv"]

            calculated_means = static_df["calc_mean"]
            calculated_cvs = static_df["calc_cv"]
            calculated_ns = static_df["calc_n"]

            per_diff_means = static_df["per_diff_mean"]
            per_diff_cvs = static_df["per_diff_cv"]

            if self.found_sd:
                reported_sds = static_df["sd"]
                calculated_sds = static_df["calc_sd"]
                per_diff_sds = static_df["per_diff_sd"]

            if len(mean) == 2:
                if ((mean[0] - mean[1]) / mean[0]) > valid_mean_diff:
                    message = stats_flag.get_message("mean_diff", level1=levels[0], level2=levels[1])
                    result.append(yellow(message, None, self.table_title, self.table_type))

            for index, level in enumerate(levels):
                cv = calc_cvs[index]
                batch = batches[index]
                per_diff_cv = per_diff_cvs[index]
                per_diff_mean = per_diff_means[index]
                reported_cv = reported_cvs[index]
                calculated_cv = calculated_cvs[index]
                reported_mean = reported_means[index]
                calculated_mean = calculated_means[index]
                calculated_n = calculated_ns[index]

                if self.found_sd:
                    per_diff_sd = per_diff_sds[index]
                    reported_sd = reported_sds[index]
                    calculated_sd = calculated_sds[index]
                    if per_diff_sd > Decimal(0.15):
                        message = stats_flag.get_message("sd_error", batch=batch, level=level,
                                                         reported_value=reported_sd,
                                                         calc_value=str(round(float(calculated_sd), 3)))
                        result.append(red(message, None, self.table_title, self.table_type))

                if per_diff_cv > valid_re_cv_diff:
                    message = stats_flag.get_message("cv_error", batch=batch, level=level, reported_value=reported_cv,
                                                     calc_value=str(round(float(calculated_cv), 3)))
                    result.append(red(message, None, self.table_title, self.table_type))

                if per_diff_mean > valid_mean_diff:
                    message = stats_flag.get_message("mean_error", batch=batch, level=level,
                                                     reported_value=reported_mean,
                                                     calc_value=str(round(float(calculated_mean), 3)))

                    result.append(red(message, None, self.table_title, self.table_type))

                if calculated_n < self.required_sample:
                    message = stats_flag.get_message("required_sample", batch=batch, level=level,
                                                     samples=self.required_sample)
                    result.append(red(message, None, self.table_title, self.table_type))
                    validation_failed_flag_count += 1

                if cv > threshold:
                    message = stats_flag.get_message("cv", batch=batch, level=level, threshold=threshold)
                    result.append(red(message, None, self.table_title, self.table_type))
                    validation_failed_flag_count += 1

                elif (threshold - Decimal(5)) <= cv <= threshold:
                    message = stats_flag.get_message("cv", batch=batch, level=level, threshold=threshold - Decimal(5))
                    result.append(yellow(message, None, self.table_title, self.table_type))

        if validation_failed_flag_count == 0:
            result.append(green(f"All values within acceptable range", None, self.table_title, self.table_type))

        self.result += result

    @staticmethod
    def get_error_dict():
        data_flag = FlagProperties("matrix_data")
        stats_flag = FlagProperties("matrix_stats")
        return data_flag, stats_flag


class CarryOver(Table, Flags):
    def __init__(self, parsed_table, template_type):
        Table.__init__(self, parsed_table, template_type)
        self.is_threshold = Decimal(self.threshold_values.get_message("is_threshold"))
        self.analyte_threshold = Decimal(self.threshold_values.get_message("analyte_threshold"))

    def validate_table_format(self):
        error_messages = self.error_messages
        required_col = {"Run ID", "Analyte %Response", "IS %Response", "Blank Analyte Peak Area", "Blank IS Peak Area",
                        "STD Analyte Peak Area", "Accepted Batch STDs and QCs Mean IS Peak Area"}
        missing_col = required_col.intersection(self.missing_col)
        if missing_col:
            message = error_messages.get_message("missing_col", col_names=", ".join(missing_col))
            self.result.append(red(message, None, self.table_title, self.table_type))
            self.valid_format = False

    def process_table(self):
        table = self.table_df
        try:
            table = utils.fill_val(table)
            table["blank_analyte_area"] = table["blank_analyte_area"].apply(utils.parse_decimal)
            table["blank_is_area"] = table["blank_is_area"].apply(utils.parse_decimal)
            table["std_analyte_area"] = table["std_analyte_area"].apply(utils.parse_decimal)
            table["batch_qc_is_mean_area"] = table["batch_qc_is_mean_area"].apply(utils.parse_decimal)
            table["analyte_response"] = table["analyte_response"].apply(utils.parse_decimal)
            table["is_response"] = table["is_response"].apply(utils.parse_decimal)
            table["calc_analyte_response"] = (table["blank_analyte_area"] / table["std_analyte_area"]) * 100
            table["calc_is_response"] = (table["blank_is_area"] / table["batch_qc_is_mean_area"]) * 100
            table["analyte_response_per_diff"] = utils.calculate_per_diff(table["analyte_response"],
                                                                          table["calc_analyte_response"])
            table["is_response_per_diff"] = utils.calculate_per_diff(table["is_response"], table["calc_is_response"])

            self.data_df = table
        except:
            message = self.error_messages.get_message("data_error")
            self.result.append(red(message, None, self.table_title, self.table_type))
            self.valid_data = False

    def validate_table_data(self):
        valid_is_diff = Decimal(self.valid_difference_values.get_message("difference_is"))
        valid_analyte_diff = Decimal(self.valid_difference_values.get_message("difference_analyte"))

        result = []
        data_flag = self.get_error_dict()
        table = self.data_df
        is_threshold = self.is_threshold
        analyte_threshold = self.analyte_threshold

        batches = table["run_id"]
        calculated_analyte_responses = table["calc_analyte_response"]
        calculated_is_responses = table["calc_is_response"]
        reported_is_responses = table["is_response"]
        reported_analyte_responses = table["analyte_response"]
        analyte_peak_areas = table["blank_analyte_area"]
        is_peak_areas = table["blank_is_area"]
        is_response_per_diffs = table["is_response_per_diff"]
        analyte_response_per_diffs = table["analyte_response_per_diff"]
        red_flags = 0
        try:
            for index, batch in enumerate(batches):
                calc_is_response = calculated_is_responses[index]
                calc_analyte_response = calculated_analyte_responses[index]
                is_peak_area = is_peak_areas[index]
                analyte_peak_area = analyte_peak_areas[index]
                is_diff = is_response_per_diffs[index]
                analyte_diff = analyte_response_per_diffs[index]
                reported_is_response = reported_is_responses[index]
                reported_analyte_response = reported_analyte_responses[index]

                if analyte_diff > valid_analyte_diff:
                    message = data_flag.get_message("analyte_error", batch=batch,
                                                    reported_value=reported_analyte_response,
                                                    calc_value=utils.format_value(calc_analyte_response),
                                                    peak_area=analyte_peak_area)
                    result.append(red(message, None, self.table_title, self.table_type))

                if is_diff > valid_is_diff:
                    message = data_flag.get_message("is_error", batch=batch, reported_value=reported_is_response,
                                                    calc_value=utils.format_value(calc_is_response),
                                                    peak_area=analyte_peak_area)
                    result.append(red(message, None, self.table_title, self.table_type))

                if calc_analyte_response > analyte_threshold:
                    message = data_flag.get_message("analyte_response", batch=batch, threshold=analyte_threshold,
                                                    peak_area=analyte_peak_area)
                    result.append(red(message, None, self.table_title, self.table_type))
                    red_flags += 1
                elif (analyte_threshold - Decimal(5)) <= calc_analyte_response <= analyte_threshold:
                    message = data_flag.get_message("analyte_response", batch=batch,
                                                    threshold=analyte_threshold - Decimal(5),
                                                    peak_area=analyte_peak_area)
                    result.append(yellow(message, None, self.table_title, self.table_type))

                if calc_is_response > is_threshold:
                    message = data_flag.get_message("is_response", batch=batch, threshold=is_threshold,
                                                    peak_area=is_peak_area)
                    result.append(red(message, None, self.table_title, self.table_type))
                    red_flags += 1
                elif (is_threshold - Decimal(2)) < calc_is_response < is_threshold:
                    message = data_flag.get_message("is_response", batch=batch, threshold=is_threshold - Decimal(2),
                                                    peak_area=is_peak_area)
                    result.append(yellow(message, None, self.table_title, self.table_type))
        except:
            pass
        if len(result) == 0:
            result.append(green(f"All values within acceptable range", None, self.table_title, self.table_type))
        if red_flags == 0:
            message = data_flag.get_message("no_carry_over")
            result.append(green(message, None, self.table_title, self.table_type))
        self.result += result

    @staticmethod
    def get_error_dict():
        data_flag = FlagProperties("co_data")
        return data_flag


class DilutionIntegrity(Table, Flags):
    def __init__(self, parsed_table, analytes, template_type):
        Table.__init__(self, parsed_table, template_type)
        self.replicate = Decimal(5)
        self.analyte = get_analyte(analytes, self.tb_title)
        self.nominal_conc = []
        self.column_info = {}
        self.column_failed = {}

    def validate_table_format(self):
        error_messages = self.error_messages
        required_col = {"Run ID"}
        missing_col = required_col.intersection(self.missing_col)
        if missing_col:
            message = error_messages.get_message("missing_col", col_names=", ".join(missing_col))
            self.result.append(red(message, None, self.table_title, self.table_type))
            self.valid_format = False

        conc_units = utils.find_units(self.table_df.columns)
        if conc_units == "":
            message = error_messages.get_message("missing_conc_unit")
            self.result.append(red(message, None, self.table_title, self.table_type))
            self.valid_format = False
        else:
            self.conc_unit = conc_units

    def process_table(self):
        table = self.table_df
        final_df = pd.DataFrame()
        calc_static_df = pd.DataFrame()
        static_df = pd.DataFrame()
        try:
            result = utils.validate_unit(table.columns, self.table_title, self.table_type)
            if len(result) > 0:
                self.result += result
                self.valid_data = False
            else:
                self.found_re = utils.find_re_column(table)
                column = list(table.columns)
                nominal_conc = utils.find_nominal_conc(column)
                dilution_factor = utils.extract_dilution_factor(column)
                table_df, static_df = utils.split_static_table(table)
                final_df = pd.DataFrame()
                calc_static_df = pd.DataFrame()

                while (len(nominal_conc) - 1) == len(dilution_factor):
                    dilution_factor.append(Decimal(1))
                count = 0
                for i, key in enumerate(nominal_conc):
                    run_df = pd.DataFrame(columns=["run_id", "column", "nominal", "dilution", "conc", "re"])
                    values = table_df[key].apply(utils.parse_decimal)
                    run_df["conc"] = table_df[key]
                    run_df["run_id"] = table_df["run_id"]
                    run_df["column"] = key
                    nominal = nominal_conc[key]
                    run_df["nominal"] = nominal
                    run_df["dilution"] = dilution_factor[i]
                    if self.found_re:
                        try:
                            column_index = table_df.columns.get_loc(key)
                            table_df.columns.values[column_index + 1] = f"re_{count}"
                            run_df["re"] = table_df[f"re_{count}"].apply(utils.parse_decimal)
                            count += 1
                        except:
                            run_df["re"] = table_df["re"].apply(utils.parse_decimal)
                    run_df = run_df.fillna("")
                    final_df = pd.concat([final_df, run_df]).reset_index(drop=True)
                    calc_static_df = pd.concat([calc_static_df, utils.build_static_df(values, nominal)]).reset_index(
                        drop=True)

                final_df = utils.fill_val(final_df)
                dummy_df = final_df.copy()
                dummy_df["conc"] = dummy_df["conc"].apply(utils.parse_decimal)
                final_df["calc_re"] = ((dummy_df["conc"] - dummy_df["nominal"]) / dummy_df["nominal"]) * 100
                if self.found_re:
                    final_df["per_diff_re"] = utils.calculate_per_diff(final_df["re"], final_df["calc_re"])

                # Static Calculation
                if not static_df.empty:
                    static_df = utils.process_static_df(static_df, self.table_type)
                    static_df.insert(loc=0, column="run_id", value=final_df["run_id"][0])

                    # Concat reported and calculated static dataframes
                    static_df, self.found_sd = utils.concat_static_df(static_df, calc_static_df)

            self.data_df = final_df.fillna("")
            self.static_df = static_df
            self.nominal_conc = list(final_df["nominal"].unique())
        except Exception as e:
            message = self.error_messages.get_message("data_error")
            self.result.append(red(message, None, self.table_title, self.table_type))
            self.valid_data = False

    def validate_table_data(self):
        valid_count_diff = Decimal(self.valid_difference_values.get_message("difference_n"))
        valid_re_cv_diff = Decimal(self.valid_difference_values.get_message("difference_re_cv"))
        valid_mean_diff = Decimal(self.valid_difference_values.get_message("difference_mean"))

        result = []
        data_flag, stats_flag = self.get_error_dict()
        table = self.data_df
        threshold = Decimal(self.threshold_values.get_message("re_cv_threshold"))
        if "lm" in self.analysis_type:
            test_threshold = threshold + Decimal(5)
        else:
            test_threshold = threshold

        batches = table["run_id"]
        all_columns = table["column"]
        all_conc = table["conc"]
        all_calc_res = table["calc_re"]
        if self.found_re:
            reported_res = table["re"]
            per_diff_res = table["per_diff_re"]
        previous_column = all_columns[0]
        red_count = 0
        failed_sample = 0
        total_sample = 0
        is_failed = False
        for index, column in enumerate(all_columns):
            batch = batches[index]
            conc = all_conc[index]
            column = all_columns[index]
            calc_re = all_calc_res[index]

            if previous_column != column:
                if Decimal(total_sample - failed_sample) < self.replicate:
                    red_count += 1
                    is_failed = True
                self.column_info[previous_column] = red_count
                self.column_failed[previous_column] = is_failed
                previous_column = column
                total_sample = 0
                failed_sample = 0
                red_count = 0
                is_failed = False

            total_sample += 1
            if self.found_re:
                reported_re = reported_res[index]
                per_diff_re = per_diff_res[index]
                re_para = {"batch": batch, "column": column, "reported_value": reported_re,
                           "calc_value": utils.format_value(calc_re), "conc": conc}
                if str(per_diff_re) != "":
                    if per_diff_re > valid_re_cv_diff:
                        if abs(calc_re) <= Decimal(1) and abs(reported_re) <= Decimal(1):
                            message = data_flag.get_message("re_rounding", **re_para)
                            result.append(yellow(message, None, self.table_title, self.table_type))
                        else:
                            message = data_flag.get_message("re_error", **re_para)
                            result.append(red(message, None, self.table_title, self.table_type))

            if str(calc_re).strip() != "":
                if abs(calc_re) > test_threshold:
                    message = data_flag.get_message("re", batch=batch, column=column, threshold=test_threshold,
                                                    conc=conc)
                    result.append(red(message, None, self.table_title, self.table_type))
                    failed_sample += 1

                elif (test_threshold - Decimal(5)) <= abs(calc_re) <= test_threshold:
                    message = data_flag.get_message("re", batch=batch, column=column,
                                                    threshold=test_threshold - Decimal(5), conc=conc)
                    result.append(yellow(message, None, self.table_title, self.table_type))

        if Decimal(total_sample - failed_sample) < self.replicate:
            red_count += 1
            is_failed = True
        self.column_info[previous_column] = red_count
        self.column_failed[previous_column] = is_failed

        static_df = self.static_df
        if not static_df.empty:
            all_column = static_df["column"]
            reported_mean = static_df["mean"]

            reported_cv = static_df["cv"]
            reported_re = static_df["re"]
            reported_n = static_df["n"]

            calculated_mean = static_df["calc_mean"]
            calculated_sd = static_df["calc_sd"]
            calculated_cv = static_df["calc_cv"]
            calculated_re = static_df["calc_re"]
            calculated_n = static_df["calc_n"]

            per_diff_mean = static_df["per_diff_mean"]

            per_diff_cv = static_df["per_diff_cv"]
            per_diff_re = static_df["per_diff_re"]
            per_diff_n = static_df["per_diff_n"]

            batches = static_df["run_id"].to_list()

            if self.found_sd:
                reported_sd = static_df["sd"]
                per_diff_sd = static_df["per_diff_sd"]
            previous_column = all_columns[0]
            red_count = 0
            for index, column in enumerate(all_column):
                batch = batches[index]
                overall_reported_mean = reported_mean[index]
                overall_reported_cv = reported_cv[index]
                overall_reported_re = reported_re[index]
                overall_reported_n = reported_n[index]

                overall_clc_mean = calculated_mean[index]
                overall_clc_cv = calculated_cv[index]
                overall_clc_re = calculated_re[index]
                overall_clc_n = calculated_n[index]

                overall_per_diff_mean = per_diff_mean[index]
                overall_per_diff_cv = per_diff_cv[index]
                overall_per_diff_re = per_diff_re[index]
                overall_per_diff_n = per_diff_n[index]

                if previous_column != column:
                    self.column_info[previous_column] += red_count
                    if self.column_failed[previous_column] and red_count == 0:
                        message = stats_flag.get_message("further_investigate", batch=batch, replicate=self.replicate,
                                                         column=previous_column)
                        result.append(yellow(message, None, self.table_title, self.table_type))

                    red_count = 0
                    previous_column = column

                replicate_para = {"column": column, "replicate": self.replicate}
                count_error_para = {"column": column, "reported_value": overall_reported_n, "calc_value": overall_clc_n}
                mean_error_para = {"column": column, "reported_value": overall_reported_mean,
                                   "calc_value": utils.format_value(overall_clc_mean)}
                re_error_para = {"column": column, "reported_value": overall_reported_re,
                                 "calc_value": utils.format_value(overall_clc_re)}
                cv_error_para = {"column": column, "reported_value": overall_reported_cv,
                                 "calc_value": utils.format_value(overall_clc_cv)}
                red_threshold_para = {"column": column, 'threshold': test_threshold}
                yellow_threshold_para = {"column": column, 'threshold': test_threshold - Decimal(5)}

                dilution_failed_message = stats_flag.get_message("dilution_fail", column=column, batch=batch)

                if self.found_sd:
                    overall_reported_sd = reported_sd[index]
                    overall_clc_sd = calculated_sd[index]
                    overall_per_diff_sd = per_diff_sd[index]

                    sd_error_para = {"column": column, "reported_value": overall_reported_sd,
                                     "calc_value": utils.format_value(overall_clc_sd)}
                    if overall_per_diff_sd > valid_re_cv_diff:
                        if abs(overall_clc_sd) <= Decimal(1) and abs(overall_reported_sd) <= Decimal(1):
                            message = stats_flag.get_message("sd_rounding", **sd_error_para)
                            result.append(yellow(message, None, self.table_title, self.table_type))

                        else:
                            message = stats_flag.get_message("sd_error", **sd_error_para)
                            result.append(red(message, None, self.table_title, self.table_type))

                if overall_clc_n < self.replicate:
                    message = stats_flag.get_message("replicate", **replicate_para)
                    result.append(red(message, None, self.table_title, self.table_type))
                    red_count += 1

                if overall_per_diff_n > valid_count_diff:
                    message = stats_flag.get_message("count_error", **count_error_para)
                    result.append(red(message, None, self.table_title, self.table_type))

                if str(overall_per_diff_mean).strip() != "":
                    if overall_per_diff_mean > valid_mean_diff:
                        if abs(overall_clc_mean) <= Decimal(1) and abs(overall_reported_mean <= Decimal(1)):
                            message = stats_flag.get_message("mean_rounding", **mean_error_para)
                            result.append(yellow(message, None, self.table_title, self.table_type))

                        else:
                            message = stats_flag.get_message("mean_error", **mean_error_para)
                            result.append(red(message, None, self.table_title, self.table_type))

                if str(overall_per_diff_cv).strip() != "":
                    if overall_per_diff_cv > valid_re_cv_diff:
                        if abs(overall_reported_cv) <= Decimal(1) and abs(overall_clc_cv) <= Decimal(1):
                            message = stats_flag.get_message("cv_rounding", **cv_error_para)
                            result.append(yellow(message, None, self.table_title, self.table_type))

                        else:
                            message = stats_flag.get_message("cv_error", **cv_error_para)
                            result.append(red(message, None, self.table_title, self.table_type))
                            red_count += 1

                if str(overall_per_diff_re).strip() != "":
                    if overall_per_diff_re > valid_re_cv_diff:
                        if abs(overall_reported_re) <= Decimal(1) and abs(overall_clc_re) <= Decimal(1):
                            message = stats_flag.get_message("re_rounding", **re_error_para)
                            result.append(yellow(message, None, self.table_title, self.table_type))

                        else:
                            message = stats_flag.get_message("re_error", **re_error_para)
                            result.append(red(message, None, self.table_title, self.table_type))
                            red_count += 1

                if abs(overall_clc_re) > test_threshold:
                    message = stats_flag.get_message("re", **red_threshold_para)
                    result.append(red(message, None, self.table_title, self.table_type))
                    result.append(red(dilution_failed_message, None, self.table_title, self.table_type))
                    red_count += 1

                elif (test_threshold - Decimal(5)) <= abs(overall_clc_re) <= test_threshold:
                    message = stats_flag.get_message("re", **yellow_threshold_para)
                    result.append(yellow(message, None, self.table_title, self.table_type))

                if abs(overall_clc_cv) > test_threshold:
                    message = stats_flag.get_message("cv", **red_threshold_para)
                    result.append(red(message, None, self.table_title, self.table_type))
                    result.append(red(dilution_failed_message, None, self.table_title, self.table_type))
                    red_count += 1

                elif (test_threshold - Decimal(5)) <= abs(overall_clc_cv) <= test_threshold:
                    message = stats_flag.get_message("cv", **yellow_threshold_para)
                    result.append(yellow(message, None, self.table_title, self.table_type))

            self.column_info[previous_column] += red_count

            if self.column_failed[previous_column] and red_count == 0:
                message = stats_flag.get_message("further_investigate", batch=batches[-1], replicate=self.replicate,
                                                 column=previous_column)
                result.append(yellow(message, None, self.table_title, self.table_type))

        for col, val in self.column_info.items():
            if val == 0:
                message = stats_flag.get_message("dil_demonstrated", batch=batches[0], column=col)
                self.result.append(green(message, None, self.table_title, self.table_type))
        if len(result) == 0:
            result.append(green('All values within nominal range', None, self.table_title, self.table_type))

        self.result += result

    @staticmethod
    def get_error_dict():
        data_flag = FlagProperties("dil_integrity_data")
        stats_flag = FlagProperties("dil_integrity_stats")
        return data_flag, stats_flag


class ISR(Table, Flags):
    def __init__(self, parsed_table, samples_table, total_samples, LLOQ, units, analytes, multiple_analyte,
                 multi_analyte_total_samples, template_type):
        Table.__init__(self, parsed_table, template_type)
        self.isr_threshold = Decimal(20)
        self.total_isr = 0
        self.isr_passed = 0
        self.total_samples = total_samples
        self.samples_table = samples_table
        self.analyte = get_analyte(analytes, self.tb_title)
        self.multiple_analyte = multiple_analyte
        self.multi_analyte_total_samples = multi_analyte_total_samples
        if isinstance(LLOQ, dict):
            LLOQ = {str(key).lower(): val for key, val in LLOQ.items()}
            self.LLOQ = LLOQ.get(self.analyte.lower())
        else:
            self.LLOQ = LLOQ
        if self.multiple_analyte:
            self.units = units.get(self.analyte, "")
        else:
            self.units = units

    @staticmethod
    def get_error_dict():
        treatment_flag = FlagProperties("isr_treatment")
        pct_diff_flag = FlagProperties("isr_percent_test")
        return treatment_flag, pct_diff_flag

    def validate_table_format(self):
        error_messages = self.error_messages
        table = self.table_df
        table_column = table.columns
        time_point = list(filter(lambda x: ("time" in str(x).lower()) & ("point" in str(x).lower()), table_column))
        if len(time_point) > 0:
            time_point_val = table[time_point[0]].to_list()
            day_val = []
            hour_val = []
            for item in time_point_val:
                values = utils.split_time_point(item)
                day_val.append(values[0])
                hour_val.append(values[1])
            table["day"] = day_val
            table["hour"] = hour_val
            self.table_df, self.missing_col = utils.format_table_data_frame(table, self.table_type)
        required_col = {"Subject", "Final Original Concentration", "ISR Concentration", "%Diff", "Day", "Hour"}
        missing_col = required_col.intersection(self.missing_col)
        if missing_col:
            message = error_messages.get_message("missing_col", col_names=", ".join(missing_col))
            self.result.append(red(message, None, self.table_title, self.table_type))
            self.valid_format = False

    def process_table(self):
        table = self.table_df
        try:
            cols = ["subject", "day", "hour"]
            table['unique_id'] = table[cols].apply(lambda row: ' '.join(row.values.astype(str)), axis=1)
            analyte_col = list(filter(lambda x: 'analyte' in str(x), list(table.columns)))
            if len(analyte_col) == 0:
                table["analyte"] = np.nan

            table["original_string"] = table["original_conc"]
            table["isr_string"] = table["isr_conc"]

            table["original_conc"] = table["original_conc"].apply(utils.parse_decimal)
            table["isr_conc"] = table["isr_conc"].apply(utils.parse_decimal)
            table["string_per_diff"] = table["percent_diff"]
            table["percent_diff"] = table["percent_diff"].apply(utils.parse_decimal)
            self.data_df = table
        except:
            message = self.error_messages.get_message("data_error")
            self.result.append(red(message, None, self.table_title, self.table_type))
            self.valid_data = False

    def validate_table_data(self):
        treatment_flag, pct_diff_flag = self.get_error_dict()
        self.treatment_groups(treatment_flag)
        self.isr_tests_percent_diff(pct_diff_flag)

    def treatment_groups(self, data_flag):
        table_rows = self.data_df
        mismatch_samples = list()
        if self.multiple_analyte:
            sample_table = self.samples_table.get(self.analyte, pd.DataFrame())
        else:
            sample_table = self.samples_table
        result = []
        treatment_flag = False

        unique_counts = table_rows['unique_id'].value_counts()
        if unique_counts.max() > 1:
            samples = list()
            for item in unique_counts.iteritems():
                if item[1] > 1:
                    samples.append(item[0])
            sample_str = ', '.join(samples)
            message = data_flag.get_message("duplicate_sample", sample=sample_str)
            result.append(yellow(message, None, self.table_title, self.table_type))

        # Check for Treatment
        table_rows = table_rows.dropna(axis=1, how='all')
        cols = list(table_rows.columns)
        if 'treatment' in cols:
            treatment_flag = True
        else:
            message = data_flag.get_message("no_treatment")
            result.append(yellow(message, None, self.table_title, self.table_type))

        if not sample_table.empty:
            for index, row in table_rows.iterrows():
                isr_id = row['unique_id']
                u_id = str(isr_id).split(" ")
                isr_conc = row["original_string"]

                found = False
                try:
                    sample_conc = sample_table[sample_table["unique_id"] == isr_id]["value"].to_list()

                    if len(sample_conc) > 0:
                        sample_conc = sample_conc[0]
                        mismatch_conc_message = data_flag.get_message("conc_mismatch", subject=u_id[0], day=u_id[1],
                                                                      hour=" ".join(u_id[2:]), isr_conc=isr_conc,
                                                                      sample_conc=sample_conc)
                        if str(sample_conc).isalpha() and str(isr_conc).isalpha():
                            if str(sample_conc).strip() != str(isr_conc).strip():
                                result.append(red(mismatch_conc_message, None, self.table_title, self.table_type))
                        else:
                            if Decimal(sample_conc) != Decimal(isr_conc):
                                result.append(red(mismatch_conc_message, None, self.table_title, self.table_type))
                    else:
                        mismatch_samples.append(isr_id)
                        message = data_flag.get_message("no_id", subject=u_id[0], day=u_id[1], hour=" ".join(u_id[2:]))
                        result.append(red(message, None, self.table_title, self.table_type))
                except Exception as e:
                    print("ISR Error", e)
                    pass
                try:
                    if treatment_flag:
                        isr_treatment = row["treatment"]
                        sample_treatment = sample_table[sample_table["unique_id"] == isr_id]["treatment"].to_list()
                        if len(sample_treatment) > 0:
                            mismatch_treatment_message = data_flag.get_message("treatment_mismatch", subject=u_id[0],
                                                                               day=u_id[1], hour=" ".join(u_id[2:]),
                                                                               isr_treatment=isr_treatment,
                                                                               sample_treatment=sample_treatment)

                            sample_treatment = sample_treatment[0]
                            if str(sample_treatment).strip() != str(isr_treatment).strip():
                                result.append(red(mismatch_treatment_message, None, self.table_title, self.table_type))
                        else:
                            mismatch_samples.append(isr_id)
                            message = data_flag.get_message("no_id", subject=u_id[0], day=u_id[1],
                                                            hour=" ".join(u_id[2:]))
                            result.append(red(message, None, self.table_title, self.table_type))
                except Exception as e:
                    pass

        if len(mismatch_samples) == 0:
            result.append(green('All ISR samples found in sample table', None, self.table_title, self.table_type))
        self.result += result

    def isr_tests_percent_diff(self, data_flag):
        if self.multiple_analyte:
            total_samples = self.multi_analyte_total_samples.get(self.analyte, 0)
        else:
            total_samples = self.total_samples

        LLOQ = self.LLOQ
        units = self.units
        table_rows = self.data_df
        result = []

        test_threshold = self.isr_threshold
        if "lm" in self.analysis_type:
            test_threshold += Decimal(10)

        try:
            limit = 3 * LLOQ
        except TypeError:
            limit = 0

        def check_bql(item):
            item = str(item).lower()
            if "bql" in item or "blq" in item or "<" in item:
                return True
            else:
                return False

        def check_aql(item):
            item = str(item).lower()
            if "aql" in item or "alq" in item or ">" in item:
                return True
            else:
                return False

        def check_na(item):
            item = str(item).lower()
            if "na" in item or "n/a" in item or "nc" in item:
                return True
            else:
                return False

        table_rows['orig_below'] = pd.DataFrame([table_rows['original_conc'] < limit]).transpose()
        table_rows['reassay_below'] = pd.DataFrame([table_rows['isr_conc'] < limit]).transpose()

        table_rows['orig_blq'] = table_rows['original_string'].apply(lambda x: check_bql(x))
        table_rows['reassay_blq'] = table_rows['isr_string'].apply(lambda x: check_bql(x))

        table_rows['orig_aql'] = table_rows['original_string'].apply(lambda x: check_aql(x))
        table_rows['reassay_aql'] = table_rows['isr_string'].apply(lambda x: check_aql(x))

        table_rows["per_diff_na"] = table_rows["string_per_diff"].apply(lambda x: check_na(x))

        table_rows['calc_percent'] = 200 * (table_rows['isr_conc'] - table_rows['original_conc']) / (
                table_rows['isr_conc'] + table_rows['original_conc'])
        table_rows['calc_percent'] = table_rows['calc_percent'].fillna("")
        blank_pct = table_rows["calc_percent"].to_list()
        blank_pct = [x for x in blank_pct if str(x).strip() == ""]
        table_rows['calc_percent'] = table_rows['calc_percent'].replace("", Decimal(0))
        table_rows['failed'] = pd.DataFrame([abs(table_rows['calc_percent']) > test_threshold]).transpose()

        try:
            table_rows['mismatch'] = (table_rows['percent_diff'].fillna(Decimal(0)) - table_rows['calc_percent'].fillna(
                Decimal(0))) / (table_rows['percent_diff'] + Decimal(0.00000001))
        except Exception as e:
            print("error", e)
            pass
        table_rows['mismatch'] = table_rows['mismatch'].apply(lambda x: abs(round(x, 2))).fillna(Decimal(0))
        table_rows['mismatch'] = pd.DataFrame([table_rows['mismatch'] > Decimal(.1)]).transpose()

        orig_below = table_rows[table_rows['orig_below']]['unique_id'].tolist()
        for u_id in orig_below:
            u_id = u_id.split(" ")
            message = data_flag.get_message("original_below_lloq", subject=u_id[0], day=u_id[1],
                                            hour=" ".join(u_id[2:]), LLOQ=LLOQ, conc_unit=units)
            result.append(red(message, None, self.table_title, self.table_type))

        reassay_below = table_rows[table_rows['reassay_below']]['unique_id'].tolist()
        for u_id in reassay_below:
            u_id = u_id.split(" ")
            message = data_flag.get_message("reassay_below_lloq", subject=u_id[0], day=u_id[1], hour=" ".join(u_id[2:]),
                                            LLOQ=LLOQ, conc_unit=units)
            result.append(red(message, None, self.table_title, self.table_type))

        reassay_blq = table_rows[table_rows['reassay_blq']]['unique_id'].tolist()
        for u_id in reassay_blq:
            u_id = u_id.split(" ")
            message = data_flag.get_message("reassay_bql", subject=u_id[0], day=u_id[1], hour=" ".join(u_id[2:]))
            result.append(red(message, None, self.table_title, self.table_type))

        orig_blq = table_rows[table_rows['orig_aql']]['unique_id'].tolist()
        for u_id in orig_blq:
            u_id = u_id.split(" ")
            message = data_flag.get_message("original_aql", subject=u_id[0], day=u_id[1], hour=" ".join(u_id[2:]))
            result.append(red(message, None, self.table_title, self.table_type))

        reassay_blq = table_rows[table_rows['reassay_aql']]['unique_id'].tolist()
        for u_id in reassay_blq:
            u_id = u_id.split(" ")
            message = data_flag.get_message("reassay_aql", subject=u_id[0], day=u_id[1], hour=" ".join(u_id[2:]))
            result.append(red(message, None, self.table_title, self.table_type))

        orig_blq = table_rows[table_rows['orig_blq']]['unique_id'].tolist()
        for u_id in orig_blq:
            u_id = u_id.split(" ")
            message = data_flag.get_message("original_bql", subject=u_id[0], day=u_id[1], hour=" ".join(u_id[2:]))
            result.append(red(message, None, self.table_title, self.table_type))

        per_diff_na = table_rows[table_rows["per_diff_na"]]
        unique_ids = per_diff_na["unique_id"].to_list()
        for u_id in unique_ids:
            try:
                per_diff = per_diff_na[per_diff_na["unique_id"] == u_id]["string_per_diff"].to_list()[0]
                u_id = u_id.split(" ")
                message = data_flag.get_message("na_pct_diff", subject=u_id[0], day=u_id[1], hour=" ".join(u_id[2:]),
                                                per_diff=per_diff)
                result.append(yellow(message, None, self.table_title, self.table_type))
            except:
                pass

        mismatch = table_rows[table_rows['mismatch']]
        mismatch = mismatch.applymap(str)
        unique_ids = mismatch["unique_id"]
        for u_id in unique_ids:
            try:
                cal_pcr = mismatch[mismatch["unique_id"] == u_id]["calc_percent"].to_list()[0]
                per_diff = mismatch[mismatch["unique_id"] == u_id]["percent_diff"].to_list()[0]
                u_id = u_id.split(" ")
                message = data_flag.get_message("mismatch_pct_diff", subject=u_id[0], day=u_id[1],
                                                hour=" ".join(u_id[2:]), per_diff=per_diff, reported_value=per_diff,
                                                calc_value=utils.format_value(cal_pcr))
                result.append(red(message, None, self.table_title, self.table_type))
            except:
                pass

        failed = table_rows[table_rows['failed']]['unique_id'].tolist()
        for u_id in failed:
            u_id = str(u_id).split(" ")
            message = data_flag.get_message("pct_diff", subject=u_id[0], day=u_id[1], hour=" ".join(u_id[2:]),
                                            threshold=test_threshold)
            result.append(red(message, None, self.table_title, self.table_type))

        failed_samples = len(failed) + len(blank_pct)
        samples = len(table_rows['unique_id'].tolist())

        if samples > 0:
            pass_pct = 100 * (samples - failed_samples) / samples
        else:
            pass_pct = 0

        if samples > 0 and total_samples > 1:
            pct_isr = 100 * Decimal(samples) / Decimal(total_samples)
            pct_isr = round(pct_isr, 1)
        else:
            pct_isr = 0

        if (pct_isr < 10) and pct_isr > 0:
            message = data_flag.get_message("isr_chosen1", pct=pct_isr)
            result.append(red(message, None, self.table_title, self.table_type))
        elif pct_isr == 0:
            pass
        else:
            message = data_flag.get_message("isr_chosen2", pct=pct_isr)
            result.append(green(message, None, self.table_title, self.table_type))

        sample_pass_message = data_flag.get_message("sample_pass", pct=round(pass_pct))
        if pass_pct < 66.67:
            result.append(red(sample_pass_message, None, self.table_title, self.table_type))
        else:
            result.append(green(sample_pass_message, None, self.table_title, self.table_type))

        self.result += result
        self.total_isr = samples
        self.isr_passed = samples - failed_samples


class Sample(Table, Flags):
    def __init__(self, parsed_table, reported_sample, analytes, multiple_analyte, LLOQ, units, template_type):
        Table.__init__(self, parsed_table, template_type)
        self.multiple_analyte = multiple_analyte
        if self.analyte is None:
            self.analyte = get_analyte(analytes, self.tb_title)
        self.samples = 0
        self.reported_sample = reported_sample
        if self.multiple_analyte and self.analyte is not None:
            self.unit = units.get(self.analyte, "")
        else:
            self.unit = units
        if isinstance(LLOQ, dict):
            LLOQ = {str(key).lower(): val for key, val in LLOQ.items()}
            self.LLOQ = LLOQ.get(self.analyte.lower())
        else:
            self.LLOQ = LLOQ

    def validate_table_format(self):
        error_messages = self.error_messages
        table = self.table_df
        table_column = table.columns
        time_point = list(filter(lambda x: ("time" in str(x).lower()) & ("point" in str(x).lower()), table_column))
        if len(time_point) > 0:
            time_point_val = table[time_point[0]].to_list()
            day_val = []
            hour_val = []
            for item in time_point_val:
                values = utils.split_time_point(item)
                day_val.append(values[0])
                hour_val.append(values[1])
            table["day"] = day_val
            table["hour"] = hour_val
            self.table_df, self.missing_col = utils.format_table_data_frame(table, self.table_type)
        required_col = {"Reported Concentration", "Subject ", "Day", "Hour"}
        missing_col = required_col.intersection(self.missing_col)
        if missing_col:
            message = error_messages.get_message("missing_col", col_names=", ".join(missing_col))
            self.result.append(red(message, None, self.table_title, self.table_type))
            self.valid_format = False

    @staticmethod
    def get_error_dict():
        data_flag = FlagProperties("sample_data")
        return data_flag

    def process_table(self):
        table = self.table_df
        try:
            cols = ["subject", "day", "hour"]
            table['unique_id'] = table[cols].apply(lambda row: ' '.join(row.values.astype(str)), axis=1)
            analyte_col = list(filter(lambda x: 'analyte' in str(x), list(table.columns)))
            if len(analyte_col) == 0:
                table["analyte"] = np.nan

            self.samples = len(table['unique_id'].to_list())
            self.data_df = table
        except:
            message = self.error_messages.get_message("data_error")
            self.result.append(red(message, None, self.table_title, self.table_type))
            self.valid_data = False

    def check_samples(self):
        data_flag = self.get_error_dict()
        if self.multiple_analyte:
            analyte_str = f"for analyte {self.analyte}"
        else:
            analyte_str = ""
        try:
            total_found = self.samples
            if self.reported_sample == total_found:
                message = data_flag.get_message("sample_match", analyte_str=analyte_str)
                self.result += [green(message, None, "Samples Table", "Samples Table")]
            elif self.reported_sample > 1:
                if (abs(self.reported_sample - total_found) / self.reported_sample) < .03:
                    message = data_flag.get_message("sample_match", analyte_str=analyte_str)
                    self.result += [green(message, None, "Samples Table", "Samples Table")]
                else:
                    if self.reported_sample > 1:
                        message = data_flag.get_message("sample_mismatch", reported_sample=self.reported_sample,
                                                        found_sample=total_found, analyte_str=analyte_str)
                        self.result += [red(message, None, "Samples Table", "Samples Table")]
        except:
            pass

    def get_concentration_units(self):
        columns = self.data_df.columns
        conc_column = [col for col in columns if 'unit' in str(col).lower()]
        if conc_column:
            conc_units = list(self.data_df[conc_column[0]].unique())
            conc_units.append(self.conc_unit)
            return conc_units
        return [self.conc_unit]

    def check_bloq(self):
        result = []
        LLOQ = self.LLOQ
        units = self.unit
        data_flag = self.get_error_dict()
        table_title = self.table_title
        table_type = self.table_type
        table_rows = self.data_df

        if self.multiple_analyte:
            analyte_str = f"for analyte {self.analyte}"
        else:
            analyte_str = ""

        if len(units) > 0:
            units = " " + str(units)

        if table_rows.empty:
            return [red('Data not found in table', None, table_title, table_type)]

        def compare_to_lloq(value, LLOQ):
            try:
                string_value = str(value).lower()
                if "llq" in string_value or "blq" in string_value or "lloq" in string_value or "lq" in string_value or "bql" in string_value or "bloq" in string_value:
                    return False

                if Decimal(value) < LLOQ:
                    return True
                else:
                    return False
            except:
                return False

        try:
            table_rows['bloq'] = table_rows['value'].fillna('BLQ').apply(lambda x: compare_to_lloq(x, LLOQ))
            bloq_rows = table_rows[table_rows['bloq']]

            bloq_rows = bloq_rows.applymap(str)
            unique_id = bloq_rows["unique_id"].to_list()
            values = bloq_rows["value"].to_list()
            for i, u_id in enumerate(unique_id):
                u_id = str(u_id).split(" ")
                message = data_flag.get_message("below_lloq_conc", sample=u_id[0], day=u_id[1], hour=u_id[2],
                                                conc=values[i], units=units, LLOQ=LLOQ, analyte_str=analyte_str)
                result.append(red(message, None, table_title, table_type))

            if len(result) == 0:
                message = data_flag.get_message("above_lloq_conc", analyte_str=analyte_str, units=units, LLOQ=LLOQ)
                result.append(green(message, None, table_title, table_type))
        except:
            pass

        self.result += result


class ReanalysisAndReassay(Table, Flags):
    def __init__(self, parsed_table, total_samples, sample_table, options, analytes, multiple_analyte,
                 multi_analyte_total_samples, template_type):
        Table.__init__(self, parsed_table, template_type)
        self.multiple_analyte = multiple_analyte
        self.total_samples = total_samples
        self.options = options
        self.sample_table = sample_table
        self.analyte = get_analyte(analytes, self.tb_title)
        self.multi_analyte_total_samples = multi_analyte_total_samples

    def validate_table_format(self):
        error_messages = self.error_messages
        table = self.table_df
        table_column = table.columns
        time_point = list(filter(lambda x: ("time" in str(x).lower()) & ("point" in str(x).lower()), table_column))
        if len(time_point) > 0:
            time_point_val = table[time_point[0]].to_list()
            day_val = []
            hour_val = []
            for item in time_point_val:
                values = utils.split_time_point(item)
                day_val.append(values[0])
                hour_val.append(values[1])
            table["day"] = day_val
            table["hour"] = hour_val
            self.table_df, self.missing_col = utils.format_table_data_frame(table, self.table_type)
        required_col = {"Reported conc", "Subject ", "Day", "Hour", "Re-assay Concentration", "Reassay Reason"}
        missing_col = required_col.intersection(self.missing_col)
        if missing_col:
            message = error_messages.get_message("missing_col", col_names=", ".join(missing_col))
            self.result.append(red(message, None, self.table_title, self.table_type))
            self.valid_format = False

    def process_table(self):
        try:
            table = self.table_df
            cols = ["subject", "day", "hour"]
            table['unique_id'] = table[cols].apply(lambda row: ' '.join(row.values.astype(str)), axis=1)
            analyte_col = list(filter(lambda x: 'analyte' in str(x), list(table.columns)))
            if len(analyte_col) == 0:
                table["analyte"] = np.nan

            self.data_df = table
        except:
            message = self.error_messages.get_message("data_error")
            self.result.append(red(message, None, self.table_title, self.table_type))
            self.valid_data = False

    def validate_table_data(self):
        if self.multiple_analyte:
            total_samples = self.multi_analyte_total_samples.get(self.analyte, 0)
        else:
            total_samples = self.total_samples

        table = self.data_df
        data_flag = self.get_error_dict()
        if self.multiple_analyte:
            sample_table = self.sample_table.get(self.analyte, pd.DataFrame())
        else:
            sample_table = self.sample_table
        result = []
        # TODO Why we have more than one sample table
        if not sample_table.empty:
            unique_id = list(sample_table['unique_id'])
            if len(unique_id) != total_samples:
                message = data_flag.get_message("no_sample_id", id=unique_id[0])
                result.append(yellow(message, None, self.table_title, self.table_type))

            # Check -  All reassay Id are there in Sample Table
            mismatch_samples = list()
            for index, row in table.iterrows():
                reassay_id = row['unique_id']
                reassay_conc = row["reported_conc"]
                if isinstance(reassay_conc, pd.Series):
                    reassay_conc = reassay_conc.tolist()[-1]

                u_id = str(reassay_id).split(" ")
                try:
                    sample_conc = sample_table[sample_table["unique_id"] == reassay_id]["value"].to_list()
                    if len(sample_conc) > 0:
                        sample_conc = sample_conc[0]
                        mismatch_conc_message = data_flag.get_message("conc_mismatch", subject=u_id[0], day=u_id[1],
                                                                      hour=" ".join(u_id[2:]),
                                                                      reassay_conc=reassay_conc,
                                                                      sample_conc=sample_conc)
                        if str(sample_conc).isalpha() and str(reassay_conc).isalpha():
                            if str(sample_conc).strip() != str(reassay_conc).strip():
                                result.append(red(mismatch_conc_message, None, self.table_title, self.table_type))
                        else:
                            if Decimal(sample_conc) != Decimal(reassay_conc):
                                result.append(red(mismatch_conc_message, None, self.table_title, self.table_type))

                    else:
                        mismatch_samples.append(reassay_id)
                        message = data_flag.get_message("no_id", subject=u_id[0], day=u_id[1], hour=" ".join(u_id[2:]))
                        result.append(red(message, None, self.table_title, self.table_type))
                except Exception as e:
                    pass

            if len(mismatch_samples) == 0:
                result.append(
                    green('All reassayed samples found in sample table', None, self.table_title, self.table_type))

        no_reason = 0
        more_appearances = 0
        freeze_thaws = self.options.get('freeze_thaws', None) or None
        no_freeze_thaw_message = data_flag.get_message("no_freeze_thaw")
        try:
            if freeze_thaws is None or freeze_thaws < 2:
                freeze_thaws = 2

                result.append(yellow(no_freeze_thaw_message, None, self.table_title, self.table_type))
        except:
            freeze_thaws = 2
            result.append(yellow(no_freeze_thaw_message, None, self.table_title, self.table_type))
        no_reason_message = data_flag.get_message("no_reason")
        try:
            if table['reason'].empty or table['reason'].dropna().empty:
                return [red(no_reason_message, None, self.table_title, self.table_type)]
        except:
            return [red(no_reason_message, None, self.table_title, self.table_type)]

        ids = list(table['unique_id'])
        if len(ids) == 0:
            message = data_flag.get_message("no_identifiers")
            return [red(message, None, self.table_title, self.table_type)]

        dupes = [item for item, count in collections.Counter(ids).items() if count > (freeze_thaws - 1)]

        if len(dupes) == 0:
            message = data_flag.get_message("less_reassay_analyzed", freeze_thaws=freeze_thaws)
            result.append(green(message, None, self.table_title, self.table_type))
        else:
            for u_id in dupes:
                message = data_flag.get_message("more_reassay_analyzed", freeze_thaws=freeze_thaws, id=u_id)
                result.append(red(message, None, self.table_title, self.table_type))

        def no_reason(reason):
            try:
                reason_str = str(reason).strip()
                if len(reason_str) == 0 or reason is None:
                    return True
                else:
                    return False
            except:
                return True

        missing_reason = table[table['reason'].apply(lambda x: no_reason(x))]['unique_id'].tolist()

        if len(missing_reason) == 0:
            message = data_flag.get_message("missing_reason")
            result.append(green(message, None, self.table_title, self.table_type))
        else:
            for u_id in missing_reason:
                message = data_flag.get_message("missing_reason2", id=u_id)
                result.append(red(message, None, self.table_title, self.table_type))

        try:
            num_samples = len(list(set(ids)))

            if total_samples > 0:
                percent_reassayed = round(100 * num_samples / total_samples, 2)
                message = data_flag.get_message("total_reassay", pct=percent_reassayed)
                if percent_reassayed < 20:
                    result.append(green(message, None, self.table_title, self.table_type))
                elif percent_reassayed < 50:
                    result.append(yellow(message, None, self.table_title, self.table_type))
                elif percent_reassayed <= 100:
                    result.append(red(message, None, self.table_title, self.table_type))
                else:
                    pass
        except:
            pass
        self.result += result

    @staticmethod
    def get_error_dict():
        data_flag = FlagProperties("reassay_data")
        return data_flag


class HookEffect(Table, Flags):
    def __init__(self, parsed_table, multiple_analyte, analytes, template_type):
        Table.__init__(self, parsed_table, template_type)
        self.replicate = Decimal(5)
        self.found_stats_re = False
        self.found_mean_response = False
        self.multiple_analyte = multiple_analyte
        self.analyte = get_analyte(analytes, self.tb_title)

    def validate_table_format(self):
        error_messages = self.error_messages
        required_col = {"Run ID"}
        missing_col = required_col.intersection(self.missing_col)
        if missing_col:
            message = error_messages.get_message("missing_col", col_names=", ".join(missing_col))
            self.result.append(red(message, None, self.table_title, self.table_type))
            self.valid_format = False

        conc_units = utils.find_units(self.table_df.columns)
        if conc_units == "":
            message = error_messages.get_message("missing_conc_unit")
            self.result.append(red(message, None, self.table_title, self.table_type))
            self.valid_format = False
            self.conc_unit = conc_units

    @staticmethod
    def get_error_dict():
        data_flag = FlagProperties("hook_data")
        stats_flag = FlagProperties("hook_stats")
        return data_flag, stats_flag

    def process_table(self):
        table = self.table_df
        self.found_re = utils.find_re_column(table)
        column = list(table.columns)
        nominal_conc = utils.find_nominal_conc(column)
        dilution_factor = utils.extract_dilution_factor(column)
        try:
            table_df, static_df = utils.split_static_table(table)
            final_df = pd.DataFrame()
            calc_static_df = pd.DataFrame()
            final_static_df = pd.DataFrame()

            while (len(nominal_conc) - 1) == len(dilution_factor):
                dilution_factor.append(Decimal(1))
            count = 0
            for i, key in enumerate(nominal_conc):
                run_df = pd.DataFrame(columns=["run_id", "column", "nominal", "dilution", "conc", "re"])
                values = table_df[key].apply(utils.parse_decimal_2)
                dilution = dilution_factor[i]
                run_df["conc"] = values
                run_df["run_id"] = table_df["run_id"]
                run_df["column"] = key
                nominal = nominal_conc[key]
                run_df["nominal"] = nominal
                run_df["dilution"] = dilution
                if self.found_re:
                    try:
                        column_index = table_df.columns.get_loc(key)
                        table_df.columns.values[column_index + 1] = f"re_{count}"
                        run_df["re"] = table_df[f"re_{count}"].apply(utils.parse_decimal)
                        count += 1
                    except:
                        run_df["re"] = table_df[f"re"].apply(utils.parse_decimal)
                run_df = run_df.fillna("")
                final_df = pd.concat([final_df, run_df]).reset_index(drop=True)
                c_static_df = utils.build_static_df(values, nominal)
                c_static_df["dilution"] = dilution
                calc_static_df = pd.concat([calc_static_df, c_static_df]).reset_index(drop=True)

            final_df = utils.fill_val(final_df)
            dummy_df = final_df.copy()
            dummy_df["conc"] = dummy_df["conc"].apply(utils.parse_decimal_2_1)
            final_df["calc_re"] = ((dummy_df["conc"] - dummy_df["nominal"]) / dummy_df["nominal"]) * 100
            if self.found_re:
                final_df["per_diff_re"] = utils.calculate_per_diff(final_df["re"], final_df["calc_re"])

            unique_nominal = list(final_df["nominal"].unique())
            unique_column = list(final_df["column"].unique())
            if len(unique_nominal) != len(unique_column):
                final_df["nominal"] = final_df["nominal"] / final_df["dilution"]
                self.found_mean_response = True
            column = self.table_df.columns
            raw_mean_col = list(filter(lambda x: "mean" in str(x).lower(), column))
            if len(raw_mean_col) > 0:
                self.found_mean_response = True
            if not static_df.empty:
                # Static Calculation
                static_df = utils.drop_blank_col(static_df)
                static_df = static_df.T
                index_column = list(static_df.index)
                index_column = list(filter(lambda x: self.conc_unit in str(x), index_column))
                static_df = utils.remove_header(static_df)
                static_df.insert(loc=0, column="run_id", value=final_df["run_id"][0])
                static_df.insert(loc=1, column="column", value=index_column)
                static_df = static_df.reset_index(drop=True)
                static_df, missing_cols = utils.format_table_data_frame(static_df, self.table_type)
                self.found_stats_re = utils.find_re_column(static_df)

                # Concat reported and calculated static dataframes
                final_static_df, self.found_sd = utils.concat_static_df(static_df, calc_static_df)
                if self.found_mean_response:
                    final_static_df["nominal"] = final_static_df["nominal"] / final_static_df["dilution"]
                final_static_df = final_static_df.sort_values("dilution").reset_index(drop=True)

            final_df = final_df.sort_values("dilution").reset_index(drop=True)
            self.data_df = final_df.fillna("")
            self.static_df = final_static_df.fillna("")
            self.found_stats_re = utils.find_stats_re_column(final_static_df)
        except:
            message = self.error_messages.get_message("data_error")
            self.result.append(red(message, None, self.table_title, self.table_type))
            self.valid_data = False

    def validate_table_data(self, LLOQ, ULOQ):
        valid_count_diff = Decimal(self.valid_difference_values.get_message("difference_n"))
        valid_re_cv_diff = Decimal(self.valid_difference_values.get_message("difference_re_cv"))
        valid_mean_diff = Decimal(self.valid_difference_values.get_message("difference_mean"))

        if self.multiple_analyte:
            ULOQ = ULOQ.get(str(self.analyte), 0)
            LLOQ = LLOQ.get(str(self.analyte), 0)

        data_flag, stats_flag = self.get_error_dict()
        result = []
        table = self.data_df
        threshold = Decimal(self.threshold_values.get_message("re_cv_threshold"))
        if "lm" in self.analysis_type:
            test_threshold = threshold + Decimal(5)
        else:
            test_threshold = threshold
        dilution_factors = list(table["dilution"].unique())
        previous_mean = []
        found_aql = False
        hook_effect = False
        AQL_PATTERN = re.compile(r"aql|alq|>uloq")
        BQL_PATTERN = re.compile(r"bql|blq|>lloq")

        hook_message = stats_flag.get_message("hook")
        for index, dilution in enumerate(dilution_factors):
            conc_value = table[table["dilution"] == dilution]["conc"].to_list()
            nominal_value = table[table["dilution"] == dilution]["nominal"].to_list()[0]
            if index == 0:
                if (AQL_PATTERN.search(str(conc_value).lower())) or nominal_value > ULOQ != 0:
                    found_aql = True
                    try:
                        previous_mean.append(statistics.mean(conc_value))
                    except:
                        previous_mean.append(ULOQ)
                else:
                    previous_mean.append(statistics.mean(conc_value))
            elif index != 0:
                if AQL_PATTERN.search(str(conc_value).lower()) and not found_aql:
                    result.append(red(hook_message, None, self.table_title, self.table_type))
                    hook_effect = True
                    break
                elif AQL_PATTERN.search(str(conc_value).lower()):
                    previous_mean.append(ULOQ)
                    continue
                elif BQL_PATTERN.search(str(conc_value).lower()):
                    calc_mean = LLOQ
                else:
                    conc_value = [utils.parse_decimal_2_1(x) for x in conc_value if
                                  utils.parse_decimal_2_1(x) is not None]
                    calc_mean = statistics.mean(conc_value)

                hook_mean = [x for x in previous_mean if x < calc_mean]
                if len(hook_mean) > 0:
                    result.append(red(hook_message, None, self.table_title, self.table_type))
                    hook_effect = True
                    break
                else:
                    previous_mean.append(calc_mean)

        table["below_lloq"] = table["nominal"] < LLOQ
        table["above_uloq"] = table["nominal"] > ULOQ

        below_lloq_df = table[table["below_lloq"]]
        above_uloq_df = table[table["above_uloq"]]

        mid_df = table[(~table["below_lloq"]) & (~table["above_uloq"])]

        if len(mid_df["column"].unique()) < 3:
            if self.found_mean_response:
                message = data_flag.get_message("calibration_range_raw_mean")
                result.append(yellow(message, None, self.table_title, self.table_type))
            else:
                message = data_flag.get_message("calibration_range")
                result.append(yellow(message, None, self.table_title, self.table_type))

        if below_lloq_df.empty:
            if self.found_mean_response:
                message = data_flag.get_message("lloq_raw_mean", LLOQ=LLOQ, unit=self.conc_unit)
                result.append(red(message, None, self.table_title, self.table_type))
            else:
                message = data_flag.get_message("lloq", LLOQ=LLOQ, unit=self.conc_unit)
                result.append(red(message, None, self.table_title, self.table_type))
        if above_uloq_df.empty:
            if self.found_mean_response:
                message = data_flag.get_message("uloq_raw_mean", ULOQ=ULOQ, unit=self.conc_unit)
                result.append(red(message, None, self.table_title, self.table_type))
            else:
                message = data_flag.get_message("uloq", ULOQ=ULOQ, unit=self.conc_unit)
                result.append(red(message, None, self.table_title, self.table_type))

        all_columns = table["column"]

        if len(set(all_columns)) < 5:
            message = data_flag.get_message("dilution_level")
            result.append(red(message, None, self.table_title, self.table_type))

        if not self.found_mean_response:
            batches = table["run_id"]

            all_conc = table["conc"]
            all_calc_res = table["calc_re"]
            if self.found_re:
                reported_res = table["re"]
                per_diff_res = table["per_diff_re"]
            for index in range(len(all_columns)):
                batch = batches[index]
                conc = all_conc[index]
                column = all_columns[index]
                calc_re = all_calc_res[index]
                if self.found_re:
                    reported_re = reported_res[index]
                    per_diff_re = per_diff_res[index]

                    if str(per_diff_re) != "":
                        if per_diff_re > valid_re_cv_diff:
                            if abs(calc_re) <= Decimal(1) and abs(reported_re) <= Decimal(1):
                                message = data_flag.get_message("re_rounding", batch=batch, reported_value=reported_re,
                                                                calc_value=utils.format_value(calc_re), column=column,
                                                                conc=conc)
                                result.append(yellow(message, None, self.table_title, self.table_type))
                            else:
                                message = data_flag.get_message("re_error", batch=batch, reported_value=reported_re,
                                                                calc_value=utils.format_value(calc_re), column=column,
                                                                conc=conc)
                                result.append(red(message, None, self.table_title, self.table_type))

                if str(calc_re).strip() != "":
                    if abs(calc_re) > test_threshold:
                        message = data_flag.get_message("re", column=column, threshold=test_threshold, batch=batch,
                                                        conc=conc)
                        result.append(red(message, None, self.table_title, self.table_type))

                    elif (test_threshold - Decimal(5)) <= abs(calc_re) <= test_threshold:
                        message = data_flag.get_message("re", column=column, threshold=test_threshold - Decimal(5),
                                                        batch=batch, conc=conc)
                        result.append(yellow(message, None, self.table_title, self.table_type))

        static_df = self.static_df
        if not static_df.empty:
            all_column = static_df["column"]
            reported_n = static_df["n"]

            reported_mean = static_df["mean"]
            calculated_mean = static_df["calc_mean"]
            per_diff_mean = static_df["per_diff_mean"]

            reported_cv = static_df["cv"]
            calculated_cv = static_df["calc_cv"]
            per_diff_cv = static_df["per_diff_cv"]

            if self.found_stats_re:
                reported_re = static_df["re"]
                per_diff_re = static_df["per_diff_re"]
            calculated_re = static_df["calc_re"]

            calculated_n = static_df["calc_n"]
            per_diff_n = static_df["per_diff_n"]

            calculated_sd = static_df["calc_sd"]
            if self.found_sd:
                reported_sd = static_df["sd"]
                per_diff_sd = static_df["per_diff_sd"]

            for index, column in enumerate(all_column):

                overall_reported_mean = reported_mean[index]
                overall_reported_cv = reported_cv[index]
                overall_reported_n = reported_n[index]

                overall_clc_mean = calculated_mean[index]
                overall_clc_cv = calculated_cv[index]
                overall_clc_re = calculated_re[index]
                overall_clc_n = calculated_n[index]

                overall_per_diff_mean = per_diff_mean[index]
                overall_per_diff_cv = per_diff_cv[index]
                overall_per_diff_n = per_diff_n[index]

                mean_error_para = {"column": column, "reported_value": overall_reported_mean,
                                   "calc_value": utils.format_value(overall_clc_mean)}
                cv_error_para = {"column": column, "reported_value": overall_reported_cv,
                                 "calc_value": utils.format_value(overall_clc_cv)}
                if self.found_stats_re:
                    overall_reported_re = reported_re[index]
                    overall_per_diff_re = per_diff_re[index]
                    re_error_para = {"column": column, "reported_value": overall_reported_re,
                                     "calc_value": utils.format_value(overall_clc_re)}
                    if str(overall_per_diff_re).strip() != "":
                        if overall_per_diff_re > valid_re_cv_diff:
                            if abs(overall_reported_re) <= Decimal(1) and abs(overall_clc_re) <= Decimal(1):
                                message = stats_flag.get_message("re_rounding", **re_error_para)
                                result.append(yellow(message, None, self.table_title, self.table_type))

                            else:
                                message = stats_flag.get_message("re_error", **re_error_para)
                                result.append(red(message, None, self.table_title, self.table_type))

                if self.found_sd:
                    overall_reported_sd = reported_sd[index]
                    overall_clc_sd = calculated_sd[index]
                    overall_per_diff_sd = per_diff_sd[index]

                    sd_error_para = {"column": column, "reported_value": overall_reported_sd,
                                     "calc_value": utils.format_value(overall_clc_sd)}
                    if overall_per_diff_sd > valid_re_cv_diff:
                        if abs(overall_clc_sd) <= Decimal(1) and abs(overall_reported_sd) <= Decimal(1):
                            message = stats_flag.get_message("sd_rounding", **sd_error_para)
                            result.append(yellow(message, None, self.table_title, self.table_type))

                        else:
                            message = stats_flag.get_message("sd_error", **sd_error_para)
                            result.append(red(message, None, self.table_title, self.table_type))

                if str(overall_per_diff_mean).strip() != "":
                    if overall_per_diff_mean > valid_mean_diff:
                        if abs(overall_clc_mean) <= Decimal(1) and abs(overall_reported_mean <= Decimal(1)):
                            message = stats_flag.get_message("mean_rounding", **mean_error_para)
                            result.append(yellow(message, None, self.table_title, self.table_type))

                        else:
                            message = stats_flag.get_message("mean_error", **mean_error_para)
                            result.append(red(message, None, self.table_title, self.table_type))

                if str(overall_per_diff_cv).strip() != "":
                    if overall_per_diff_cv > valid_re_cv_diff:
                        if abs(overall_reported_cv) <= Decimal(1) and abs(overall_clc_cv) <= Decimal(1):
                            message = stats_flag.get_message("cv_rounding", **cv_error_para)
                            result.append(yellow(message, None, self.table_title, self.table_type))

                        else:
                            message = stats_flag.get_message("cv_error", **cv_error_para)
                            result.append(red(message, None, self.table_title, self.table_type))

                if overall_clc_n < self.replicate:
                    message = stats_flag.get_message("required_sample", column=column, replicate=self.replicate)
                    result.append(red(message, None, self.table_title, self.table_type))

                if overall_per_diff_n > valid_count_diff:
                    message = stats_flag.get_message("count_error", column=column, reported_value=overall_reported_n,
                                                     calc_value=overall_clc_n)
                    result.append(red(message, None, self.table_title, self.table_type))

                if abs(overall_clc_cv) > test_threshold:
                    message = stats_flag.get_message("cv", column=column, threshold=test_threshold)
                    result.append(red(message, None, self.table_title, self.table_type))

                elif (test_threshold - Decimal(5)) <= abs(overall_clc_cv) <= test_threshold:
                    message = stats_flag.get_message("cv", column=column, threshold=test_threshold - Decimal(5))
                    result.append(yellow(message, None, self.table_title, self.table_type))

        if not hook_effect:
            message = stats_flag.get_message("no_hook")
            result.append(green(message, None, self.table_title, self.table_type))

        if len(result) == 0:
            result.append(green(f"All values within acceptable range", None, self.table_title, self.table_type))

        self.result += result


class BloodStability(Table, Flags):

    def __init__(self, parsed_table, analytes, template_type):
        Table.__init__(self, parsed_table, template_type)
        self.analyte = get_analyte(analytes, self.tb_title)
        if 'sm' in self.analysis_type:
            self.required_sample = Decimal(3)
        else:
            self.required_sample = Decimal(0)

    def validate_table_format(self):
        error_messages = self.error_messages
        required_col = {"Run ID", "%Difference", "Peak Area Ratio", "Mean Peak Area Ratio", "Time",
                        "Concentration Level"}
        missing_col = required_col.intersection(self.missing_col)
        if missing_col:
            message = error_messages.get_message("missing_col", col_names=", ".join(missing_col))
            self.result.append(red(message, None, self.table_title, self.table_type))
            self.valid_format = False

    @staticmethod
    def get_error_dict():
        data_flag = FlagProperties("blood_stability")
        return data_flag

    def process_table(self):
        table = self.table_df
        self.found_cv = utils.find_cv_column(table)
        self.found_sd = utils.find_sd_column(table)
        calc_static_df = pd.DataFrame()
        static_df = pd.DataFrame()
        try:
            table = utils.fill_val(table)
            table["peak_area"] = table["peak_area"].apply(utils.parse_decimal)
            qc_levels = list(table["qc_level"].unique())
            table_times = list(table["time"].unique())
            for level in qc_levels:
                df = table[table["qc_level"] == level]
                run_id = list(df["run_id"].unique())
                table_times = list(df["time"].unique())
                sub_calc_static_df = pd.DataFrame()
                sub_static_df = pd.DataFrame()
                for time in table_times:
                    values = df[df["time"] == time]["peak_area"]
                    sub_calc_static_df = pd.concat([sub_calc_static_df, utils.build_static_df(values, 0)]).reset_index(
                        drop=True)
                sub_calc_static_df["time"] = table_times
                sub_calc_static_df["qc_level"] = level
                sub_calc_static_df["run_id"] = run_id[0]
                calc_static_df = pd.concat([calc_static_df, sub_calc_static_df]).reset_index(drop=True)

                mean = df["mean_peak_area"]
                mean = [x for x in mean if str(x).strip() != ""]
                if self.found_cv:
                    cv = df["cv"]
                    cv = [x for x in cv if str(x).strip() != ""]
                else:
                    cv = [Decimal(0)] * len(table_times)
                if self.found_sd:
                    sd = df["sd"]
                    sd = [x for x in sd if str(x).strip() != ""]
                else:
                    sd = [Decimal(0)] * len(table_times)
                res = [Decimal(0)] * len(table_times)
                n = [Decimal(0)] * len(table_times)
                sub_static_df["mean"] = mean
                sub_static_df["sd"] = sd
                sub_static_df["cv"] = cv
                sub_static_df["re"] = res
                sub_static_df["n"] = n
                dummy_df = sub_static_df.copy()
                dummy_df["mean"] = dummy_df["mean"].apply(utils.parse_decimal)
                sub_static_df["per_diff"] = (abs(dummy_df["mean"] - Decimal(mean[0])) / Decimal(mean[0])) * 100
                static_df = pd.concat([static_df, sub_static_df]).reset_index(drop=True)
            static_df, found_sd = utils.concat_static_df(static_df, calc_static_df)
            self.static_df = static_df
        except:
            message = self.error_messages.get_message("data_error")
            self.result.append(red(message, None, self.table_title, self.table_type))
            self.valid_data = False

    @staticmethod
    def check_required_level(levels):
        levels = str(levels).lower()
        if ("low" in levels and "high" in levels) or ("lqc" in levels and "hqc" in levels):
            return True
        else:
            return False

    def validate_table_data(self):
        valid_count_diff = Decimal(self.valid_difference_values.get_message("difference_n"))
        valid_re_cv_diff = Decimal(self.valid_difference_values.get_message("difference_re_cv"))
        valid_mean_diff = Decimal(self.valid_difference_values.get_message("difference_mean"))

        stats_flag = self.get_error_dict()
        result = []
        static_df = self.static_df
        threshold = Decimal(self.threshold_values.get_message("re_cv_threshold"))
        if "lm" in self.analysis_type:
            threshold += Decimal(5)

        if not static_df.empty:
            calc_cvs = static_df["calc_cv"].to_list()
            levels = static_df["qc_level"]
            batches = static_df["run_id"]
            times = static_df["time"]

            reported_means = static_df["mean"]
            reported_cvs = static_df["cv"]

            calculated_means = static_df["calc_mean"]
            calculated_cvs = static_df["calc_cv"]
            calculated_ns = static_df["calc_n"]

            per_diff_means = static_df["per_diff_mean"]
            per_diff_cvs = static_df["per_diff_cv"]

            per_diffs = static_df["per_diff"]

            # check required level
            required_levels = self.check_required_level(levels)
            if not required_levels:
                message = stats_flag.get_message("required_level", batch=batches[0])
                self.result.append(red(message, None, self.table_title, self.table_type))

            if self.found_sd:
                reported_sds = static_df["sd"]
                calculated_sds = static_df["calc_sd"]
                per_diff_sds = static_df["per_diff_sd"]

            for index, level in enumerate(levels):
                calc_cv = calc_cvs[index]
                batch = batches[index]
                per_diff_cv = per_diff_cvs[index]
                per_diff_mean = per_diff_means[index]
                reported_cv = reported_cvs[index]
                calculated_cv = calculated_cvs[index]
                reported_mean = reported_means[index]
                calculated_mean = calculated_means[index]
                calculated_n = calculated_ns[index]
                per_diff = per_diffs[index]
                time = times[index]

                if self.found_sd:
                    per_diff_sd = per_diff_sds[index]
                    reported_sd = reported_sds[index]
                    calculated_sd = calculated_sds[index]
                    if per_diff_sd > valid_re_cv_diff:
                        message = stats_flag.get_message("sd_error", batch=batch, level=level, time=time,
                                                         reported_value=reported_sd,
                                                         calc_value=utils.format_value(calculated_sd))
                        result.append(red(message, None, self.table_title, self.table_type))

                if per_diff_cv > valid_re_cv_diff:
                    message = stats_flag.get_message("cv_error", batch=batch, level=level, time=time,
                                                     reported_value=reported_cv,
                                                     calc_value=utils.format_value(calculated_cv))
                    result.append(red(message, None, self.table_title, self.table_type))

                if per_diff_mean > valid_mean_diff:
                    message = stats_flag.get_message("mean_error", batch=batch, level=level, time=time,
                                                     reported_value=reported_mean,
                                                     calc_value=utils.format_value(calculated_mean))
                    result.append(red(message, None, self.table_title, self.table_type))

                if calc_cv > threshold:
                    message = stats_flag.get_message("cv", batch=batch, level=level, time=time, threshold=threshold)
                    result.append(red(message, None, self.table_title, self.table_type))
                elif (threshold - Decimal(5)) <= calc_cv <= threshold:
                    message = stats_flag.get_message("cv", batch=batch, level=level, time=time,
                                                     threshold=threshold - Decimal(5))
                    result.append(yellow(message, None, self.table_title, self.table_type))

                if per_diff > threshold:
                    message = stats_flag.get_message("difference", batch=batch, level=level, time=time,
                                                     threshold=threshold)
                    result.append(red(message, None, self.table_title, self.table_type))
                elif (threshold - Decimal(5)) <= per_diff <= threshold:
                    message = stats_flag.get_message("difference", batch=batch, level=level, time=time,
                                                     threshold=threshold - Decimal(5))
                    result.append(yellow(message, None, self.table_title, self.table_type))

                if calculated_n < self.required_sample:
                    message = stats_flag.get_message("required_sample", batch=batch, level=level, time=time,
                                                     required_sample=self.required_sample)
                    result.append(red(message, None, self.table_title, self.table_type))
        if len(result) == 0:
            result.append(green(f"All values within acceptable range", None, self.table_title, self.table_type))

        self.result += result

    def validate_table_title(self):
        table_title = self.table_title
        if "blood" in str(table_title).lower():
            pass
        else:
            self.result.append(
                yellow(f"Table title does not specify [Blood] as the biological matrix", None, table_title,
                       self.table_type))


class StockStability(Table, Flags):

    def __init__(self, parsed_table, analytes, template_type):
        Table.__init__(self, parsed_table, template_type)
        self.found_diff = False
        self.per_diff_threshold = Decimal(10)
        self.sample_count = Decimal(3)
        self.analyte = get_analyte(analytes, self.tb_title)

    def validate_table_format(self):
        error_messages = self.error_messages
        required_col = {"Run ID"}
        missing_col = required_col.intersection(self.missing_col)
        if missing_col:
            message = error_messages.get_message("missing_col", col_names=", ".join(missing_col))
            self.result.append(red(message, None, self.table_title, self.table_type))
            self.valid_format = False

    @staticmethod
    def get_error_dict():
        data_flag = FlagProperties("stock_stability")
        return data_flag

    def process_table(self):
        table = self.table_df
        calc_static_df = pd.DataFrame()
        try:
            table, static_df = utils.split_static_table(table)
            tab_column = list(table.columns)
            column = list(filter(lambda x: "peak" in str(x).lower() or "area" in str(x).lower(), tab_column))
            if not column:
                column = utils.find_nominal_conc(tab_column)
                if column:
                    column = column.keys()
                else:
                    raise Exception

            for col in column:
                values = table[col].apply(utils.parse_decimal).to_list()
                calc_static_df = pd.concat([calc_static_df, utils.build_static_df(values, 0)]).reset_index(drop=True)
            calc_static_df["column"] = column

            static = static_df.T
            static = utils.remove_header(static)
            static.insert(loc=0, column='run_id', value=table['run_id'][0])
            static, missing_cols = utils.format_table_data_frame(static, self.table_type)
            if "%Difference" not in missing_cols:
                self.found_diff = True

            if self.found_diff:
                static["diff"] = static["diff"].apply(utils.parse_decimal)

            static_df, self.found_sd = utils.concat_static_df(static, calc_static_df)
            self.static_df = static_df
            self.found_cv = utils.find_cv_column(static_df)
            self.found_mean = utils.find_mean_column(static_df)
            self.found_re = utils.find_stats_re_column(static_df)

        except:
            message = self.error_messages.get_message("data_error")
            self.result.append(red(message, None, self.table_title, self.table_type))
            self.valid_data = False

    def validate_table_data(self):
        validation_failed_flag_count = 0

        valid_count_diff = Decimal(self.valid_difference_values.get_message("difference_n"))
        valid_re_cv_diff = Decimal(self.valid_difference_values.get_message("difference_re_cv"))
        valid_mean_diff = Decimal(self.valid_difference_values.get_message("difference_mean"))

        stats_flag = self.get_error_dict()
        result = []
        threshold = Decimal(self.threshold_values.get_message("re_cv_threshold"))
        if "lm" in self.analysis_type:
            test_threshold = threshold + Decimal(5)
        else:
            test_threshold = threshold
        static_df = self.static_df

        if not static_df.empty:
            if self.found_diff:
                difference = static_df["diff"].to_list()
                for diff in difference:
                    if str(diff).strip() != "" and diff is not None:
                        if abs(diff) > self.per_diff_threshold:
                            message = stats_flag.get_message("difference", threshold=self.per_diff_threshold)
                            result.append(red(message, None, self.table_title, self.table_type))
                            validation_failed_flag_count += 1
            else:
                message = stats_flag.get_message("no_difference", )
                result.append(red(message, None, self.table_title, self.table_type))

            all_column = static_df["column"]
            reported_mean = static_df["mean"]
            calculated_mean = static_df["calc_mean"]
            per_diff_mean = static_df["per_diff_mean"]
            batches = static_df["run_id"]

            if self.found_cv:
                reported_cv = static_df["cv"]
                per_diff_cv = static_df["per_diff_cv"]
            calculated_cv = static_df["calc_cv"]

            if self.found_re:
                reported_re = static_df["re"]
                per_diff_re = static_df["per_diff_re"]
            calculated_re = static_df["calc_re"]

            if self.found_n:
                reported_n = static_df["n"]
                per_diff_n = static_df["per_diff_n"]
            calculated_n = static_df["calc_n"]

            if self.found_sd:
                reported_sd = static_df["sd"]
                per_diff_sd = static_df["per_diff_sd"]
            calculated_sd = static_df["calc_sd"]

            for index, column in enumerate(all_column):
                batch = batches[index]
                overall_reported_mean = reported_mean[index]
                overall_clc_mean = calculated_mean[index]
                overall_clc_cv = calculated_cv[index]
                overall_clc_re = calculated_re[index]
                overall_clc_n = calculated_n[index]
                overall_clc_sd = calculated_sd[index]
                overall_per_diff_mean = per_diff_mean[index]

                mean_error_para = {"column": column, "reported_value": overall_reported_mean,
                                   "calc_value": utils.format_value(overall_clc_mean)}

                red_threshold = {"column": column, "threshold": test_threshold}
                yellow_threshold = {"column": column, "threshold": test_threshold - Decimal(5)}

                if self.found_cv:
                    overall_reported_cv = reported_cv[index]
                    overall_per_diff_cv = per_diff_cv[index]
                    cv_error_para = {"column": column, "reported_value": overall_reported_cv,
                                     "calc_value": utils.format_value(overall_clc_cv)}

                    if overall_per_diff_cv > valid_re_cv_diff:
                        if abs(overall_reported_cv) <= Decimal(1) and abs(overall_clc_cv) <= Decimal(1):
                            message = stats_flag.get_message("cv_rounding", **cv_error_para)
                            result.append(yellow(message, None, self.table_title, self.table_type))

                        else:
                            message = stats_flag.get_message("cv_error", **cv_error_para)
                            result.append(red(message, None, self.table_title, self.table_type))

                if self.found_re:
                    overall_reported_re = reported_re[index]
                    overall_per_diff_re = per_diff_re[index]
                    re_error_para = {"column": column, "reported_value": overall_reported_re,
                                     "calc_value": utils.format_value(overall_clc_re)}

                    if overall_per_diff_re > valid_re_cv_diff:
                        if abs(overall_reported_re) <= Decimal(1) and abs(overall_clc_re) <= Decimal(1):
                            message = stats_flag.get_message("re_rounding", **re_error_para)
                            result.append(yellow(message, None, self.table_title, self.table_type))

                        else:
                            message = stats_flag.get_message("re_error", **re_error_para)
                            result.append(red(message, None, self.table_title, self.table_type))

                if self.found_n:
                    overall_reported_n = reported_n[index]
                    overall_per_diff_n = per_diff_n[index]
                    count_error_para = {"column": column, "reported_value": overall_reported_n,
                                        "calc_value": overall_clc_n}
                    if overall_per_diff_n > valid_count_diff:
                        message = stats_flag.get_message("count_error", **count_error_para)
                        result.append(red(message, None, self.table_title, self.table_type))

                if self.found_sd:
                    overall_reported_sd = reported_sd[index]
                    overall_per_diff_sd = per_diff_sd[index]
                    sd_error_para = {"column": column, "reported_value": overall_reported_sd,
                                     "calc_value": utils.format_value(overall_clc_sd)}
                    if overall_per_diff_sd > valid_re_cv_diff:
                        if abs(overall_clc_sd) <= Decimal(1) and abs(overall_reported_sd) <= Decimal(1):
                            message = stats_flag.get_message("sd_rounding", **sd_error_para)
                            result.append(yellow(message, None, self.table_title, self.table_type))

                        else:
                            message = stats_flag.get_message("sd_error", **sd_error_para)
                            result.append(red(message, None, self.table_title, self.table_type))

                if overall_per_diff_mean > valid_mean_diff:
                    if abs(overall_clc_mean) <= Decimal(1) and abs(overall_reported_mean) <= Decimal(1):
                        message = stats_flag.get_message("mean_rounding", **mean_error_para)
                        result.append(yellow(message, None, self.table_title, self.table_type))

                    else:
                        message = stats_flag.get_message("mean_error", **mean_error_para)
                        result.append(red(message, None, self.table_title, self.table_type))

                try:
                    if abs(overall_clc_re) > test_threshold:
                        message = stats_flag.get_message("re", **red_threshold)
                        result.append(red(message, None, self.table_title, self.table_type))
                        validation_failed_flag_count += 1

                    elif (test_threshold - Decimal(5)) <= abs(overall_clc_re) <= test_threshold:
                        message = stats_flag.get_message("re", **yellow_threshold)
                        result.append(yellow(message, None, self.table_title, self.table_type))
                except:
                    pass

                if abs(overall_clc_cv) > test_threshold:
                    message = stats_flag.get_message("cv", **red_threshold)
                    result.append(red(message, None, self.table_title, self.table_type))
                    validation_failed_flag_count += 1

                elif (test_threshold - Decimal(5)) <= abs(overall_clc_cv) <= test_threshold:
                    message = stats_flag.get_message("cv", **yellow_threshold)
                    result.append(yellow(message, None, self.table_title, self.table_type))

                if overall_clc_n < self.sample_count:
                    required_count_para = {"batch": batch, "column": column, "sample_count": self.sample_count}
                    message = stats_flag.get_message("required_count", **required_count_para)
                    result.append(red(message, None, self.table_title, self.table_type))
                    validation_failed_flag_count += 1

            if validation_failed_flag_count == 0:
                result.append(green(f"All values within acceptable range", None, self.table_title, self.table_type))
        self.result += result


class Interference(Table, Flags):
    def __init__(self, parsed_table, analytes, template_type):
        Table.__init__(self, parsed_table, template_type)
        self.interference = 0
        self.interference_threshold = Decimal(5)
        self.analyte_interference = False
        self.is_interference = False
        self.found_std_area = False
        self.analyte_threshold = Decimal(self.threshold_values.get_message("analyte_threshold"))
        self.is_threshold = Decimal(self.threshold_values.get_message("is_threshold"))
        self.required_sample = 6
        if "lm" in self.analysis_type:
            self.required_sample = 10
        self.found_duplicate_run = False
        self.analyte = get_analyte(analytes, self.tb_title)

    def check_interference_type(self):
        column = self.table_df.columns
        analyte_column = [x for x in column if
                          ("%" in str(x).lower() and "analyte" in str(x).lower()) or "analyte_response" in str(
                              x).lower()]
        if len(analyte_column) > 0:
            self.analyte_interference = True
        else:
            self.is_interference = True

    def validate_table_format(self):
        error_messages = self.error_messages
        self.check_interference_type()
        if self.analyte_interference:
            required_col = {"Run ID", "Sample", "Analyte Peak Area", "Analyte LLOQ (STD1) Peak Area",
                            "%Analyte Response"}
        else:
            required_col = {"Run ID", "Sample", "Interference Sample IS Peak Area",
                            "Mean Interference Sample IS Peak Area", "Mean STD and QC IS Peak Area", "%Interference"}
        missing_col = required_col.intersection(self.missing_col)
        if missing_col:
            message = error_messages.get_message("missing_col", col_names=", ".join(missing_col))
            self.result.append(red(message, None, self.table_title, self.table_type))
            self.valid_format = False

    @staticmethod
    def get_error_dict():
        data_flag = FlagProperties("interference_test")
        return data_flag

    def process_table(self):
        table = self.table_df
        try:
            if self.analyte_interference:
                table["analyte_area"] = np.vectorize(utils.parse_decimal)(table["analyte_area"])
                table["lloq_analyte_area"] = np.vectorize(utils.parse_decimal)(table["lloq_analyte_area"])
                table["analyte_response"] = np.vectorize(utils.parse_decimal)(table["analyte_response"])
                table["calc_analyte_response"] = (table["analyte_area"] / table["lloq_analyte_area"]) * 100
            else:
                table["is_area"] = np.vectorize(utils.parse_decimal)(table["is_area"])
                table["is_area_mean"] = np.vectorize(utils.parse_decimal)(table["is_area_mean"])
                table["std_qc_is_area_mean"] = np.vectorize(utils.parse_decimal)(table["std_qc_is_area_mean"])
                table["interference"] = np.vectorize(utils.parse_decimal)(table["interference"])
                table["calc_is_mean"] = Decimal(table["is_area"].mean())
                try:
                    table["std_qc_is_area"] = np.vectorize(utils.parse_decimal)(table["std_qc_is_area"])
                    table["calc_std_qc_mean"] = Decimal(table["std_qc_is_area"].mean())
                    table["calc_interference"] = (table["calc_is_mean"] / table["calc_std_qc_mean"]) * 100
                    self.found_std_area = True
                except KeyError:
                    table["calc_interference"] = (table["calc_is_mean"] / table["std_qc_is_area_mean"]) * 100
                    self.result.append(yellow(f"Could not identify column: STD and QC IS Peak Area in the table", None,
                                              self.table_title, self.table_type))
            self.data_df = table.fillna("")
        except:
            message = self.error_messages.get_message("data_error")
            self.result.append(red(message, None, self.table_title, self.table_type))
            self.valid_data = False

    def validate_table_data(self):
        data_flag = self.get_error_dict()
        result = []
        table = self.data_df

        if self.analyte_interference:
            if len(self.table_df["sample"].to_list()) < self.required_sample:
                message = data_flag.get_message("required_sample", required_sample=self.required_sample)
                result.append(red(message, None, self.table_title, self.table_type))

            batch_tables, duplicate_id = utils.split_run_df(utils.fill_val(self.table_df))
            for batch_table in batch_tables:
                b_table = batch_table["table"]
                duplicate_sample = b_table[b_table["sample"].duplicated()]["sample"]
                for sample in duplicate_sample:
                    message = data_flag.get_message("sample_more_than_once", sample=sample)
                    result.append(red(message, None, self.table_title, self.table_type))

            analyte_response = table["calc_analyte_response"]
            samples = table["sample"]
            total_samples = table["sample"].count()
            batches = table["run_id"]
            failed_sample = 0
            for index, response in enumerate(analyte_response):
                sample = samples[index]
                batch = batches[index]
                if str(response) != "":
                    if response > self.analyte_threshold:
                        failed_sample += 1
                        message = data_flag.get_message("analyte_response", batch=batch, sample=sample,
                                                        threshold=self.analyte_threshold)
                        result.append(red(message, None, self.table_title, self.table_type))

        elif self.is_interference:
            batches = table["run_id"]
            calc_is_means = table["calc_is_mean"].unique()
            reported_is_area_means = table["is_area_mean"].unique()
            reported_interferences = table["interference"]
            calc_interferences = table["calc_interference"].unique()
            if self.found_std_area:
                calc_std_qc_means = table["calc_std_qc_mean"].unique()
                reported_std_qc_means = table["std_qc_is_area_mean"]
            for index, calc_is_mean in enumerate(calc_is_means):
                batch = batches[index]
                reported_is_area_mean = reported_is_area_means[index]
                calc_interference = calc_interferences[index]
                reported_interference = reported_interferences[index]

                try:
                    if calc_interference == 0 and reported_interference == 0:
                        pass
                    else:
                        if abs((calc_interference - reported_interference) / calc_interference) > Decimal(0.15):
                            message = data_flag.get_message("interference_error", batch=batch,
                                                            reported_value=reported_interference,
                                                            calc_value=utils.format_value(calc_interference))
                            result.append(red(message, None, self.table_title, self.table_type))
                except ZeroDivisionError:
                    if abs((reported_interference - calc_interference) / reported_interference) > Decimal(0.15):
                        message = data_flag.get_message("interference_error", batch=batch,
                                                        reported_value=reported_interference,
                                                        calc_value=utils.format_value(calc_interference))
                        result.append(red(message, None, self.table_title, self.table_type))

                try:
                    if calc_is_mean == 0 and reported_is_area_mean == 0:
                        pass
                    else:
                        if abs((calc_is_mean - reported_is_area_mean) / calc_is_mean) > Decimal(0.15):
                            message = data_flag.get_message("mean_interference_is_area_error", batch=batch,
                                                            reported_value=reported_is_area_mean,
                                                            calc_value=utils.format_value(calc_is_mean))
                            result.append(red(message, None, self.table_title, self.table_type))

                except ZeroDivisionError:
                    if abs((reported_is_area_mean - calc_is_mean) / reported_is_area_mean) > Decimal(0.15):
                        message = data_flag.get_message("mean_interference_is_area_error", batch=batch,
                                                        reported_value=reported_is_area_mean,
                                                        calc_value=utils.format_value(calc_is_mean))
                        result.append(red(message, None, self.table_title, self.table_type))

                if self.found_std_area:
                    calc_std_qc_mean = calc_std_qc_means[index]
                    reported_std_qc_mean = reported_std_qc_means[index]
                    if abs((calc_std_qc_mean - reported_std_qc_mean) / calc_std_qc_mean) > Decimal(0.15):
                        message = data_flag.get_message("mean_std_qc_area_error", batch=batch,
                                                        reported_value=reported_std_qc_mean,
                                                        calc_value=utils.format_value(calc_std_qc_mean))
                        result.append(red(message, None, self.table_title, self.table_type))

                if str(calc_interference) != "":
                    if calc_interference > self.is_threshold:
                        message = data_flag.get_message("interference_above", batch=batch, threshold=self.is_threshold)
                        result.append(red(message, None, self.table_title, self.table_type))

                    elif (self.is_threshold - Decimal(1)) < calc_interference <= self.is_threshold:
                        message = data_flag.get_message("interference_above", batch=batch,
                                                        threshold=self.is_threshold - Decimal(1))
                        result.append(yellow(message, None, self.table_title, self.table_type))
                    elif calc_interference <= (self.is_threshold - Decimal(1)):
                        message = data_flag.get_message("interference_below", batch=batch,
                                                        threshold=self.is_threshold - Decimal(1))
                        result.append(green(message, None, self.table_title, self.table_type))
        if len(result) == 0:
            result.append(green(f"All values within acceptable range", None, self.table_title, self.table_type))

        self.result += result


class ValidationSummary(Table):
    def __init__(self, parsed_table, possible_anticoagulants, possible_sample_matrices, possible_sample_preps,
                 possible_sample_species, possible_regressions, dil_integrity_dict, ULOQ, LLOQ, qc_header_values,
                 cs_header_values, ap_stats, units, options, multiple_analytes, analytes, template_type):
        Table.__init__(self, parsed_table, template_type)
        self.multiple_analytes = multiple_analytes
        self.analytes = analytes
        self.possible_anticoagulants = possible_anticoagulants
        self.possible_sample_matrices = possible_sample_matrices
        self.possible_sample_species = possible_sample_species
        self.possible_sample_preps = possible_sample_preps
        self.possible_regressions = possible_regressions
        self.dil_integrity_dict = dil_integrity_dict
        self.qc_header_values = qc_header_values
        self.cs_header_values = cs_header_values
        self.options = options
        self.units = units
        self.ap_stats = ap_stats
        self.ULOQ = ULOQ
        self.LLOQ = LLOQ
        self.found_sample_matrix = []
        self.found_species = []
        self.found_anticoagulant = []
        self.found_sample_prep = []
        self.found_regression = ''
        self.found_cs_concentrations = []
        self.found_LLOQ = 0
        self.found_ULOQ = 0
        self.found_qc_concentrations = []
        self.found_freeze_thaws = None
        self.found_short_stability = None
        self.found_frozen_stability = None
        self.found_long_stability = None
        self.found_extract_stability = None
        self.found_run_length = 0
        self.found_intra_precisions = {}
        self.found_intra_accuracies = {}
        self.found_inter_precisions = {}
        self.found_inter_accuracies = {}
        self.found_inter_tes = {}
        self.found_dilution_integrity_conc = None
        self.found_dilution_linearity_conc = None
        self.analyte = get_analyte(analytes, self.tb_title)

    def get_dil_integrity_conc(self, value):
        if self.units:
            value_array = value.split(self.units.lower())
            if len(value_array) >= 2:
                raw_conc_val = value_array[0]
                conc_val = utils.parse_last_decimal(raw_conc_val)
                return conc_val
        return

    def process_table(self):
        table = self.table_df
        for index, row in table.iterrows():
            row = row.to_list()
            search_text = str(row[0]).lower()
            value_text = str(row[-1]).lower()
            value_text2 = str(row[1]).lower()
            value_text2 = str(row[1]).lower().replace("-", " ")
            search_array = search_text.split()
            value_array = value_text.split()
            value_array2 = value_text2.split()
            for item in self.possible_anticoagulants:
                if len(item.split(" ")) > 1:
                    if item.lower() in search_text or item.lower() in value_text or item.lower() in value_text2:
                        self.found_anticoagulant.append(item)
                elif item.lower() == 'heparin' and (
                        'sodium heparin' in self.found_anticoagulant or 'na heparin' in self.found_anticoagulant):
                    pass
                elif item.lower() in search_array or item.lower() in value_array or item.lower() in value_array2:
                    self.found_anticoagulant.append(item)

            for item in self.possible_sample_matrices:
                if len(item.split(" ")) > 1:
                    if item.lower() in search_text or item.lower() in value_text or item.lower() in value_text2:
                        self.found_sample_matrix.append(item)
                elif item.lower() in search_array or item.lower() in value_array or item.lower() in value_array2:
                    self.found_sample_matrix.append(item)

            for item in self.possible_sample_preps:
                if len(item.split(" ")) > 1:
                    if item.lower() in search_text or item.lower() in value_text or item.lower() in value_text2:
                        self.found_sample_prep.append(item)
                elif item.lower() in search_array or item.lower() in value_array or item.lower() in value_array2:
                    self.found_sample_prep.append(item)

            for item in self.possible_sample_species:
                if len(item.split(" ")) > 1:
                    if item.lower() in search_text or item.lower() in value_text or item.lower() in value_text2:
                        self.found_species.append(item)
                elif item.lower() in search_array or item.lower() in value_array or item.lower() in value_array2:
                    self.found_species.append(item)

            if len(self.found_regression) == 0:
                for item in self.possible_regressions:
                    if item.lower() in search_array or item.lower() in value_array or item.lower() in value_array2:
                        self.found_regression = item

            if self.found_freeze_thaws is None:
                if (
                        "freeze" in search_text and "thaw" in search_text) or "f/t" in search_text or "f / t" in search_text or "ft" in search_text:
                    search_list = ['freeze', 'thaw', 'f/t', 'f / t', 'ft', 'freeze/thaw', "cycle"]
                    if self.found_freeze_thaws is None:
                        index_list = [j for j, x in enumerate(value_array) for y in search_list if y in x]
                        if len(index_list) > 0:
                            end_index = max(index_list)
                            for z in value_array[:end_index + 1]:
                                if self.found_freeze_thaws is None:
                                    self.found_freeze_thaws = utils.parse_signed_int(z)
                                    if self.found_freeze_thaws is not None:
                                        if self.found_freeze_thaws < 0:
                                            self.found_freeze_thaws = None
                                else:
                                    break

            if self.found_run_length < 1:
                for i, cell in enumerate(row):
                    cell_string = str(cell).lower().strip()
                    if self.found_run_length < 1:
                        if "run" in cell_string and "length" in cell_string:
                            search_row = row[i:]
                            for search_cell in search_row:
                                search_cell_array = str(search_cell).lower().strip().split()
                                search_list = ['run', 'length']
                                index_list = []
                                if self.found_run_length < 1:
                                    index_list = [j for j, x in enumerate(search_cell_array) for y in search_list if
                                                  y in x]
                                    if len(index_list) > 0:
                                        begin_index = max(index_list)
                                        for z in search_cell_array[begin_index:]:
                                            if self.found_run_length < 1:
                                                self.found_run_length = utils.parse_signed_int(z)
                                                if self.found_run_length is None:
                                                    self.found_run_length = 0
                                                elif self.found_run_length < 1:
                                                    self.found_run_length = 0

            if self.found_dilution_integrity_conc is None:
                if "dilution" in search_text and "integrity" in search_text:
                    self.found_dilution_integrity_conc = self.get_dil_integrity_conc(value_text)
                    for val in value_array:
                        if self.found_dilution_integrity_conc is None:
                            value = utils.parse_decimal(val)
                            if value is not None:
                                self.found_dilution_integrity_conc = value
                        else:
                            break

            if self.found_dilution_linearity_conc is None:
                if "dilution" in search_text and "linearity" in search_text:
                    for val in value_array:
                        if self.found_dilution_linearity_conc is None:
                            value = utils.parse_decimal(val)
                            if value is not None:
                                self.found_dilution_linearity_conc = value
                        else:
                            break

            if self.found_short_stability is None:
                if ("bench" in search_text and "stab" in search_text) or (
                        "short" in search_text and "stab" in search_text) or (
                        "bench" in search_text and "top" in search_text) or (
                        "thaw" in search_text and "matrix" in search_text):
                    if self.found_short_stability is None:
                        hour_index = [j for j, x in enumerate(value_array) if "hour" in x or "hr" in x]
                        if len(hour_index) > 0:
                            for v in value_array[:hour_index[0] + 1]:
                                if self.found_short_stability is None:
                                    value = utils.parse_float(v)
                                    if value is not None and value >= 0:
                                        self.found_short_stability = value
                                else:
                                    break

            if self.found_long_stability is None:
                if ("long" in search_text and "stab" in search_text) or (
                        "froze" in search_text and "stab" in search_text) or (
                        "long" in search_text and "term" in search_text):
                    if self.found_long_stability is None:
                        day_index = [j for j, x in enumerate(value_array) if "day" in x]
                        if len(day_index) > 0:
                            for v in value_array[:day_index[0] + 1]:
                                if self.found_long_stability is None:
                                    value = utils.parse_signed_int(v)
                                    if value is not None and value >= 0:
                                        self.found_long_stability = value
                                else:
                                    break

            if self.found_extract_stability is None:
                if ("extract" in search_text and "stab" in search_text) or (
                        "proces" in search_text and "stab" in search_text):
                    if self.found_extract_stability is None:
                        hour_index = [j for j, x in enumerate(value_array) if "hour" in x or "hr" in x]
                        if len(hour_index) > 0:
                            for v in value_array[:hour_index[0] + 1]:
                                if self.found_extract_stability is None:
                                    value = utils.parse_float(v)
                                    if value is not None and value >= 0:
                                        self.found_extract_stability = value
                                else:
                                    break

            if ("cs" in search_text and "conc" in search_text) or (
                    "st" in search_text and "curve" in search_text) or (
                    "st" in search_text and "conc" in search_text) or (
                    "calibration" in search_text and "range" in search_text) or (
                    "analytical" in search_text and "range" in search_text):
                try:
                    self.found_cs_concentrations += [utils.parse_decimal(s) for s in value_text.split() if
                                                     utils.parse_decimal(s) is not None]
                except:
                    continue
            elif ("q" in search_text and "conc" in search_text) or ("q" in search_text and "levels" in search_text):
                try:
                    self.found_qc_concentrations += [utils.parse_decimal(s) for s in value_text.split() if
                                                     utils.parse_decimal(s) is not None]
                except:
                    continue
            elif ("lloq" in search_text and "%" not in search_text) or (
                    "low" in search_text and "quant" in search_text):
                try:
                    self.found_LLOQ = utils.parse_decimal(value_text)
                except:
                    continue
            elif ("uloq" in search_text and "%" not in search_text) or (
                    "up" in search_text and "quant" in search_text):
                try:
                    self.found_ULOQ = utils.parse_decimal(value_text)
                except:
                    continue
            elif "frozen storage" in search_text:
                try:
                    temp_array = [utils.parse_int(s) for s in value_text.split() if utils.parse_int(s) is not None]
                    self.found_frozen_stability = temp_array[0]
                except:
                    continue
            elif "lloq" in search_text and (
                    "intra" in search_text or "within" in search_text) and "ac" in search_text:
                try:
                    self.found_intra_accuracies["lloq_accuracy"] = [utils.parse_decimal(s) for s in value_array if
                                                                    utils.parse_decimal(s) is not None]
                except Exception as e:
                    continue
            elif "uloq" in search_text and (
                    "intra" in search_text or "within" in search_text) and "ac" in search_text:
                try:
                    self.found_intra_accuracies["uloq_accuracy"] = [utils.parse_decimal(s) for s in
                                                                    value_text.split() if
                                                                    utils.parse_decimal(s) is not None]
                except Exception as e:
                    continue
            elif ("intra" in search_text or "within" in search_text) and "ac" in search_text:
                try:
                    self.found_intra_accuracies["lmh_accuracy"] = [utils.parse_decimal(s) for s in
                                                                   value_text.split() if
                                                                   utils.parse_decimal(s) is not None]
                except Exception as e:
                    continue
            elif "lloq" in search_text and (
                    "inter" in search_text or "overall" in search_text or "total" in search_text) and "ac" in search_text:
                try:
                    self.found_inter_accuracies["lloq_accuracy"] = [utils.parse_decimal(s) for s in
                                                                    value_text.split() if
                                                                    utils.parse_decimal(s) is not None]
                except:
                    continue
            elif "uloq" in search_text and (
                    "inter" in search_text or "overall" in search_text or "total" in search_text) and "ac" in search_text:
                try:
                    self.found_inter_accuracies["uloq_accuracy"] = [utils.parse_decimal(s) for s in
                                                                    value_text.split() if
                                                                    utils.parse_decimal(s) is not None]
                except:
                    continue
            elif (
                    "inter" in search_text or "overall" in search_text or "total" in search_text) and "ac" in search_text:
                try:
                    self.found_inter_accuracies["lmh_accuracy"] = [utils.parse_decimal(s) for s in
                                                                   value_text.split() if
                                                                   utils.parse_decimal(s) is not None]
                except:
                    continue
            elif "lloq" in search_text and (
                    "intra" in search_text or "within" in search_text) and "pre" in search_text:
                try:
                    self.found_intra_precisions["lloq_precision"] = [utils.parse_decimal(s) for s in
                                                                     value_text.split() if
                                                                     utils.parse_decimal(s) is not None]

                except:
                    continue
            elif "uloq" in search_text and (
                    "intra" in search_text or "within" in search_text) and "pre" in search_text:
                try:
                    self.found_intra_precisions["uloq_precision"] = [utils.parse_decimal(s) for s in
                                                                     value_text.split() if
                                                                     utils.parse_decimal(s) is not None]

                except:
                    continue
            elif ("intra" in search_text or "within" in search_text) and "pre" in search_text:
                try:
                    self.found_intra_precisions["lmh_precision"] = [utils.parse_decimal(s) for s in
                                                                    value_text.split() if
                                                                    utils.parse_decimal(s) is not None]

                except:
                    continue
            elif "lloq" in search_text and (
                    "inter" in search_text or "overall" in search_text or "total" in search_text) and "pre" in search_text:
                try:
                    self.found_inter_precisions["lloq_precision"] = [utils.parse_decimal(s) for s in
                                                                     value_text.split() if
                                                                     utils.parse_decimal(s) is not None]
                except:
                    continue
            elif "uloq" in search_text and (
                    "inter" in search_text or "overall" in search_text or "total" in search_text) and "pre" in search_text:
                try:
                    self.found_inter_precisions["uloq_precision"] = [utils.parse_decimal(s) for s in
                                                                     value_text.split() if
                                                                     utils.parse_decimal(s) is not None]
                except:
                    continue
            elif (
                    "inter" in search_text or "overall" in search_text or "total" in search_text) and "pre" in search_text:
                try:
                    self.found_inter_precisions["lmh_precision"] = [utils.parse_decimal(s) for s in
                                                                    value_text.split() if
                                                                    utils.parse_decimal(s) is not None]
                except:
                    continue
            elif "lloq" in search_text and (
                    "inter" in search_text or "overall" in search_text or "total" in search_text) and (
                    "error" in search_text or "%te" in search_text):
                try:
                    self.found_inter_tes["lloq_te"] = [utils.parse_decimal(s) for s in value_text.split() if
                                                       utils.parse_decimal(s) is not None]
                except:
                    continue
            elif "uloq" in search_text and (
                    "inter" in search_text or "overall" in search_text or "total" in search_text) and (
                    "error" in search_text or "%te" in search_text):
                try:
                    self.found_inter_tes["uloq_te"] = [utils.parse_decimal(s) for s in value_text.split() if
                                                       utils.parse_decimal(s) is not None]
                except:
                    continue
            elif ("inter" in search_text or "overall" in search_text or "total" in search_text) and (
                    "error" in search_text or "%te" in search_text) and "range" in search_text:
                try:
                    self.found_inter_tes["lmh_te"] = [utils.parse_decimal(s) for s in value_text.split() if
                                                      utils.parse_decimal(s) is not None]
                except:
                    continue

        if len(self.found_regression) > 0:
            if "quad" in self.found_regression:
                self.found_regression = "quadratic"
            elif "lin" in self.found_regression:
                self.found_regression = "linear"
            elif "pl" in self.found_regression:
                self.found_regression = self.found_regression.replace("pl", "PL")

        if len(self.found_cs_concentrations) > 0:
            self.found_cs_concentrations = list(set(self.found_cs_concentrations))

    @staticmethod
    def get_error_dict():
        data_flag = FlagProperties("validation_summary")
        return data_flag

    def validate_table_data(self):
        data_flag = self.get_error_dict()
        result = []
        if self.multiple_analytes:
            self.ap_stats = self.ap_stats.get(self.analyte, {})
            self.dil_integrity_dict = self.dil_integrity_dict.get(self.analyte, {})
            self.qc_header_values = self.qc_header_values.get(self.analyte, [])
            self.cs_header_values = self.cs_header_values.get(self.analyte, [])
            self.LLOQ = self.LLOQ.get(self.analyte, Decimal(0))
            self.ULOQ = self.ULOQ.get(self.analyte, Decimal(0))
            self.units = self.units.get(self.analyte, "")

        ap_intra_precisions = self.ap_stats.get("intra_precisions", None)
        ap_intra_accuracies = self.ap_stats.get("intra_accuracies", None)

        ap_inter_precisions = self.ap_stats.get("inter_precisions", None)
        ap_inter_accuracies = self.ap_stats.get("inter_accuracies", None)
        ap_inter_tes = self.ap_stats.get("inter_tes", None)

        dil_integrity_conc = self.dil_integrity_dict.get("nominal", None)
        qc_header_values = self.qc_header_values
        cs_header_values = self.cs_header_values

        LLOQ = self.LLOQ
        ULOQ = self.ULOQ

        sample_matrix = self.options.get("sample_matrix", None)
        anticoagulant = self.options.get("anticoagulant", None)
        sample_prep = self.options.get("samplePrep", None)
        species = self.options.get("species", None)
        regression_weighting = self.options.get("regression_model", None)
        extraction_stability = self.options.get("extraction_stability", None)
        short_stability = self.options.get("benchtop_stability", None)
        long_stability = self.options.get("lts80", None)
        freeze_thaws = self.options.get("freeze_thaws", None)

        units = ""
        if len(self.units) > 0:
            units = " " + str(self.units)

        try:
            if dil_integrity_conc is not None and self.found_dilution_integrity_conc != 0 and len(
                    dil_integrity_conc) != 0:
                if dil_integrity_conc[0] != self.found_dilution_integrity_conc:
                    message = data_flag.get_message("dil_integrity_conc",
                                                    reported_value=self.found_dilution_integrity_conc,
                                                    user_value=dil_integrity_conc[0], unit=units)
                    result.append(red(message, None, self.table_title, self.table_type))
        except Exception as e:
            pass

        try:
            if len(sample_matrix) > 0 and sample_matrix is not None:
                if len(self.found_sample_matrix) > 0 and self.found_sample_matrix is not None:
                    for item in set(self.found_sample_matrix):
                        if len(item) > 0 and item is not None:
                            message = data_flag.get_message("sample_matrix", reported_value=item,
                                                            user_value=sample_matrix)
                            if str(item).lower().strip() != str(sample_matrix).lower().strip():
                                result.append(red(message, None, self.table_title, self.table_type))
                else:
                    message = data_flag.get_message("no_sample_matrix")
                    result.append(yellow(message, None, self.table_title, self.table_type))
        except:
            pass

        try:
            if len(anticoagulant) > 0 and anticoagulant is not None:
                if len(self.found_anticoagulant) > 0 and self.found_anticoagulant is not None:
                    for item in self.found_anticoagulant:
                        item = item.replace("-", "")
                        if len(item) > 0 and item is not None:
                            message = data_flag.get_message("anticoagulant", reported_value=item,
                                                            user_value=anticoagulant)
                            if str(item).lower().strip() != str(anticoagulant).lower().strip():
                                result.append(red(message, None, self.table_title, self.table_type))
                else:
                    message = data_flag.get_message("no_anticoagulant")
                    result.append(yellow(message, None, self.table_title, self.table_type))
        except:
            pass

        try:
            if len(sample_prep) > 0 and sample_prep is not None:
                if len(self.found_sample_prep) > 0 and self.found_sample_prep is not None:
                    for item in self.found_sample_prep:
                        if len(item) > 0 and item is not None:
                            message = data_flag.get_message("sample_pre", reported_value=item, user_value=sample_prep)
                            if str(item).lower().strip() != str(sample_prep).lower().strip():
                                result.append(red(message, None, self.table_title, self.table_type))
                else:
                    message = data_flag.get_message("no_sample_pre")
                    result.append(yellow(message, None, self.table_title, self.table_type))
        except:
            pass

        try:
            if len(species) > 0 and species is not None:
                if len(self.found_species) > 0 and self.found_species is not None:
                    for item in set(self.found_species):
                        if len(item) > 0 and item is not None:
                            if str(item).lower().strip() != str(species).lower().strip():
                                message = data_flag.get_message("species", reported_value=item, user_value=species)
                                result.append(red(message, None, self.table_title, self.table_type))
                else:
                    message = data_flag.get_message("no_species")
                    result.append(yellow(message, None, self.table_title, self.table_type))
        except:
            pass

        try:
            if len(regression_weighting) > 0 and regression_weighting is not None:
                if len(self.found_regression) > 0 and self.found_regression is not None:
                    if str(self.found_regression).lower().strip() != str(regression_weighting).lower().strip():
                        message = data_flag.get_message("regression", reported_value=self.found_regression,
                                                        user_value=regression_weighting)
                        result.append(red(message, None, self.table_title, self.table_type))
                else:
                    message = data_flag.get_message("no_regression")
                    result.append(yellow(message, None, self.table_title, self.table_type))
        except:
            pass

        try:
            if long_stability is not None:
                if self.found_long_stability is not None:
                    if Decimal(self.found_long_stability) != Decimal(long_stability):
                        message = data_flag.get_message("long_term_stability", reported_value=self.found_long_stability,
                                                        user_value=long_stability)
                        result.append(red(message, None, self.table_title, self.table_type))
                else:
                    message = data_flag.get_message("no_long_term_stability")
                    result.append(yellow(message, None, self.table_title, self.table_type))
        except:
            pass

        try:
            if extraction_stability is not None:
                if self.found_extract_stability is not None:
                    if Decimal(self.found_extract_stability) != Decimal(extraction_stability):
                        message = data_flag.get_message("extraction_stability",
                                                        reported_value=self.found_extract_stability,
                                                        user_value=extraction_stability)
                        result.append(red(message, None, self.table_title, self.table_type))
                else:
                    message = data_flag.get_message("no_extraction_stability")
                    result.append(yellow(message, None, self.table_title, self.table_type))
        except:
            pass

        try:
            if freeze_thaws is not None:
                if self.found_freeze_thaws is not None:
                    if self.found_freeze_thaws != freeze_thaws:
                        message = data_flag.get_message("freeze_thaw", reported_value=self.found_freeze_thaws,
                                                        user_value=freeze_thaws)
                        result.append(red(message, None, self.table_title, self.table_type))
                else:
                    message = data_flag.get_message("no_freeze_thaw")
                    result.append(yellow(message, None, self.table_title, self.table_type))
        except:
            pass

        try:
            if short_stability is not None:
                if self.found_short_stability is not None:
                    if int(round(self.found_short_stability)) != int(round(short_stability)):
                        message = data_flag.get_message("short_term_stability",
                                                        reported_value=self.found_short_stability,
                                                        user_value=short_stability)
                        result.append(red(message, None, self.table_title, self.table_type))
                else:
                    message = data_flag.get_message("no_short_term_stability")
                    result.append(yellow(message, None, self.table_title, self.table_type))
        except:
            pass

        try:
            if LLOQ > 0 and LLOQ is not None:
                if self.found_LLOQ > 0 and self.found_LLOQ is not None:
                    if self.found_LLOQ != LLOQ:
                        message = data_flag.get_message("lloq", reported_value=self.found_LLOQ, user_value=LLOQ)
                        result.append(red(message, None, self.table_title, self.table_type))
                else:
                    message = data_flag.get_message("no_lloq")
                    result.append(yellow(message, None, self.table_title, self.table_type))
        except:
            pass

        try:
            if "lm" in self.analysis_type:
                if ULOQ > 0 and ULOQ is not None:
                    if self.found_ULOQ > 0 and self.found_ULOQ is not None:
                        if self.found_ULOQ != ULOQ:
                            message = data_flag.get_message("uloq", reported_value=self.found_ULOQ, user_value=ULOQ)
                            result.append(red(message, None, self.table_title, self.table_type))
                    else:
                        message = data_flag.get_message("no_uloq")
                        result.append(yellow(message, None, self.table_title, self.table_type))
        except:
            pass

        try:
            if len(cs_header_values) > 0 and cs_header_values is not None:
                cs_header_values = [i for i in cs_header_values if i]
                if len(self.found_cs_concentrations) > 0 and self.found_cs_concentrations is not None:
                    found_cs_concentrations = [i for i in self.found_cs_concentrations if i]
                    if len(found_cs_concentrations) == 2:
                        if max(cs_header_values) not in found_cs_concentrations:
                            message = data_flag.get_message("cs_header", value=max(cs_header_values), unit=units)
                            result.append(red(message, None, self.table_title, self.table_type))

                        if min(cs_header_values) not in found_cs_concentrations:
                            message = data_flag.get_message("cs_header", value=min(cs_header_values), unit=units)
                            result.append(red(message, None, self.table_title, self.table_type))
                    else:
                        for item in cs_header_values:
                            if item not in found_cs_concentrations:
                                message = data_flag.get_message("cs_header", value=item, unit=units)
                                result.append(red(message, None, self.table_title, self.table_type))
                else:
                    message = data_flag.get_message("no_cs_header")
                    result.append(yellow(message, None, self.table_title, self.table_type))
        except:
            pass

        try:
            if len(qc_header_values) > 0 and qc_header_values is not None:
                qc_header_values = [i for i in qc_header_values if i]
                if len(self.found_qc_concentrations) > 0 and self.found_qc_concentrations is not None:
                    found_qc_concentrations = [i for i in self.found_qc_concentrations if i]
                    if len(found_qc_concentrations) == 2:
                        if max(qc_header_values) not in found_qc_concentrations:
                            message = data_flag.get_message("qc_header", value=max(qc_header_values), unit=units)
                            result.append(red(message, None, self.table_title, self.table_type))
                        if min(qc_header_values) not in found_qc_concentrations:
                            message = data_flag.get_message("qc_header", value=min(qc_header_values), unit=units)
                            result.append(red(message, None, self.table_title, self.table_type))
                    else:
                        for item in qc_header_values:
                            if item not in found_qc_concentrations:
                                message = data_flag.get_message("qc_header", value=item, unit=units)
                                result.append(red(message, None, self.table_title, self.table_type))
                else:
                    message = data_flag.get_message("no_qc_header")
                    result.append(yellow(message, None, self.table_title, self.table_type))
        except:
            pass

        try:
            if ap_intra_accuracies is not None:
                ap_lloq_intra_accuracy = ap_intra_accuracies.get("lloq_accuracy", None)
                if ap_lloq_intra_accuracy is not None:
                    if self.found_intra_accuracies is not None and len(self.found_intra_accuracies) > 0:
                        lloq_intra_accuracy = self.found_intra_accuracies.get("lloq_accuracy", None)
                        if lloq_intra_accuracy is not None and len(lloq_intra_accuracy) == 2:
                            if (round(lloq_intra_accuracy[0], 1) != round(ap_lloq_intra_accuracy[0], 1)) or (
                                    round(lloq_intra_accuracy[1], 1) != round(ap_lloq_intra_accuracy[1], 1)):
                                message = data_flag.get_message("lloq_intra_accuracy_mismatch",
                                                                calc_value1=ap_lloq_intra_accuracy[0],
                                                                calc_value2=ap_lloq_intra_accuracy[1],
                                                                reported_value1=lloq_intra_accuracy[0],
                                                                reported_value2=lloq_intra_accuracy[1])
                                result.append(red(message, None, self.table_title, self.table_type))
                        else:
                            message = data_flag.get_message("lloq_intra_accuracy_missing",
                                                            calc_value1=ap_lloq_intra_accuracy[0],
                                                            calc_value2=ap_lloq_intra_accuracy[1])

                            result.append(red(message, None, self.table_title, self.table_type))
                    else:
                        message = data_flag.get_message("lloq_intra_accuracy_missing",
                                                        calc_value1=ap_lloq_intra_accuracy[0],
                                                        calc_value2=ap_lloq_intra_accuracy[1])
                        result.append(red(message, None, self.table_title, self.table_type))

                ap_uloq_intra_accuracy = ap_intra_accuracies.get("uloq_accuracy", None)
                if ap_uloq_intra_accuracy is not None:
                    if self.found_intra_accuracies is not None and len(self.found_intra_accuracies) > 0:
                        uloq_intra_accuracy = self.found_intra_accuracies.get("uloq_accuracy", None)
                        if uloq_intra_accuracy is not None and len(uloq_intra_accuracy) == 2:
                            if (round(uloq_intra_accuracy[0], 1) != round(ap_uloq_intra_accuracy[0], 1)) or (
                                    round(uloq_intra_accuracy[1], 1) != round(ap_uloq_intra_accuracy[1], 1)):
                                message = data_flag.get_message("uloq_intra_accuracy_mismatch",
                                                                calc_value1=ap_uloq_intra_accuracy[0],
                                                                calc_value2=ap_uloq_intra_accuracy[1],
                                                                reported_value1=uloq_intra_accuracy[0],
                                                                reported_value2=uloq_intra_accuracy[1])
                                result.append(red(message, None, self.table_title, self.table_type))
                        else:
                            message = data_flag.get_message("uloq_intra_accuracy_missing",
                                                            calc_value1=ap_uloq_intra_accuracy[0],
                                                            calc_value2=ap_uloq_intra_accuracy[1])
                            result.append(red(message, None, self.table_title, self.table_type))
                    else:
                        message = data_flag.get_message("uloq_intra_accuracy_missing",
                                                        calc_value1=ap_uloq_intra_accuracy[0],
                                                        calc_value2=ap_uloq_intra_accuracy[1])
                        result.append(red(message, None, self.table_title, self.table_type))

                ap_lmh_intra_accuracy = ap_intra_accuracies.get("lmh_accuracy", None)
                if ap_lmh_intra_accuracy is not None:
                    if self.found_intra_accuracies is not None and len(self.found_intra_accuracies) > 0:
                        lmh_intra_accuracy = self.found_intra_accuracies.get("lmh_accuracy", None)
                        if lmh_intra_accuracy is not None and len(lmh_intra_accuracy) == 2:
                            if (round(lmh_intra_accuracy[0], 1) != round(ap_lmh_intra_accuracy[0], 1)) or (
                                    round(lmh_intra_accuracy[1], 1) != round(ap_lmh_intra_accuracy[1], 1)):
                                message = data_flag.get_message("lmh_intra_accuracy_mismatch",
                                                                calc_value1=ap_lmh_intra_accuracy[0],
                                                                calc_value2=ap_lmh_intra_accuracy[1],
                                                                reported_value1=lmh_intra_accuracy[0],
                                                                reported_value2=lmh_intra_accuracy[1])
                                result.append(red(message, None, self.table_title, self.table_type))
                        else:
                            message = data_flag.get_message("lmh_intra_accuracy_missing",
                                                            calc_value1=ap_lmh_intra_accuracy[0],
                                                            calc_value2=ap_lmh_intra_accuracy[1])
                            result.append(red(message, None, self.table_title, self.table_type))
                    else:
                        message = data_flag.get_message("lmh_intra_accuracy_missing",
                                                        calc_value1=ap_lmh_intra_accuracy[0],
                                                        calc_value2=ap_lmh_intra_accuracy[1])
                        result.append(red(message, None, self.table_title, self.table_type))
        except:
            pass
        try:
            if ap_inter_accuracies is not None:
                ap_lloq_inter_accuracy = ap_inter_accuracies.get("lloq_accuracy", None)
                if ap_lloq_inter_accuracy is not None:
                    if self.found_inter_accuracies is not None and len(self.found_inter_accuracies) > 0:
                        lloq_inter_accuracy = self.found_inter_accuracies.get("lloq_accuracy", None)
                        if lloq_inter_accuracy is not None and len(lloq_inter_accuracy) == 1:
                            if round(lloq_inter_accuracy[0], 1) != round(ap_lloq_inter_accuracy, 1):
                                message = data_flag.get_message("lloq_inter_accuracy_mismatch",
                                                                calc_value=ap_lloq_inter_accuracy,
                                                                reported_value=lloq_inter_accuracy[0])
                                result.append(red(message, None, self.table_title, self.table_type))
                        else:
                            message = data_flag.get_message("lloq_inter_accuracy_missing",
                                                            calc_value=ap_lloq_inter_accuracy)
                            result.append(red(message, None, self.table_title, self.table_type))
                    else:
                        message = data_flag.get_message("lloq_inter_accuracy_missing",
                                                        calc_value=ap_lloq_inter_accuracy)
                        result.append(red(message, None, self.table_title, self.table_type))

                ap_uloq_inter_accuracy = ap_inter_accuracies.get("uloq_accuracy", None)
                if ap_uloq_inter_accuracy is not None:
                    if self.found_inter_accuracies is not None and len(self.found_inter_accuracies) > 0:
                        uloq_inter_accuracy = self.found_inter_accuracies.get("uloq_accuracy", None)
                        if uloq_inter_accuracy is not None and len(uloq_inter_accuracy) == 1:
                            if round(uloq_inter_accuracy[0], 1) != round(ap_uloq_inter_accuracy, 1):
                                message = data_flag.get_message("uloq_inter_accuracy_mismatch",
                                                                calc_value=ap_uloq_inter_accuracy,
                                                                reported_value=uloq_inter_accuracy[0])
                                result.append(red(message, None, self.table_title, self.table_type))
                        else:
                            message = data_flag.get_message("uloq_inter_accuracy_missing",
                                                            calc_value=ap_uloq_inter_accuracy)
                            result.append(red(message, None, self.table_title, self.table_type))
                    else:
                        message = data_flag.get_message("uloq_inter_accuracy_missing",
                                                        calc_value=ap_uloq_inter_accuracy)
                        result.append(red(message, None, self.table_title, self.table_type))

                ap_lmh_inter_accuracy = ap_inter_accuracies.get("lmh_accuracy", None)
                if ap_lmh_inter_accuracy is not None:
                    if self.found_inter_accuracies is not None and len(self.found_inter_accuracies) > 0:
                        lmh_inter_accuracy = self.found_inter_accuracies.get("lmh_accuracy", None)
                        if lmh_inter_accuracy is not None and len(lmh_inter_accuracy) == 2:
                            if (round(lmh_inter_accuracy[0], 1) != round(ap_lmh_inter_accuracy[0], 1)) or (
                                    round(lmh_inter_accuracy[1], 1) != round(ap_lmh_inter_accuracy[1], 1)):
                                message = data_flag.get_message("lmh_inter_accuracy_mismatch",
                                                                calc_value1=ap_lmh_inter_accuracy[0],
                                                                calc_value2=ap_lmh_inter_accuracy[1],
                                                                reported_value1=lmh_inter_accuracy[0],
                                                                reported_value2=lmh_inter_accuracy[1])
                                result.append(red(message, None, self.table_title, self.table_type))
                        else:
                            message = data_flag.get_message("lmh_inter_accuracy_missing",
                                                            calc_value1=ap_lmh_inter_accuracy[0],
                                                            calc_value2=ap_lmh_inter_accuracy[1])
                            result.append(red(message, None, self.table_title, self.table_type))
                    else:
                        message = data_flag.get_message("lmh_inter_accuracy_missing",
                                                        calc_value1=ap_lmh_inter_accuracy[0],
                                                        calc_value2=ap_lmh_inter_accuracy[1])
                        result.append(red(message, None, self.table_title, self.table_type))
        except:
            pass

        try:
            if ap_inter_precisions is not None:
                ap_lloq_inter_precision = ap_inter_precisions.get("lloq_precision", None)
                if ap_lloq_inter_precision is not None:
                    if self.found_inter_precisions is not None and len(self.found_inter_precisions) > 0:
                        lloq_inter_precision = self.found_inter_precisions.get("lloq_precision", None)
                        if lloq_inter_precision is not None and len(lloq_inter_precision) == 1:
                            if round(lloq_inter_precision[0], 1) != round(ap_lloq_inter_precision, 1):
                                message = data_flag.get_message("lloq_inter_precision_mismatch",
                                                                calc_value=ap_lloq_inter_precision,
                                                                reported_value=lloq_inter_precision[0])
                                result.append(red(message, None, self.table_title, self.table_type))
                        else:
                            message = data_flag.get_message("lloq_inter_precision_missing",
                                                            calc_value=ap_lloq_inter_precision)
                            result.append(red(message, None, self.table_title, self.table_type))
                    else:
                        message = data_flag.get_message("lloq_inter_precision_missing",
                                                        calc_value=ap_lloq_inter_precision)
                        result.append(red(message, None, self.table_title, self.table_type))

                ap_uloq_inter_precision = ap_inter_precisions.get("uloq_precision", None)
                if ap_uloq_inter_precision is not None:
                    if self.found_inter_precisions is not None and len(self.found_inter_precisions) > 0:
                        uloq_inter_precision = self.found_inter_precisions.get("uloq_precision", None)
                        if uloq_inter_precision is not None and len(uloq_inter_precision) == 1:
                            if round(uloq_inter_precision[0], 1) != round(ap_uloq_inter_precision, 1):
                                message = data_flag.get_message("uloq_inter_precision_mismatch",
                                                                calc_value=ap_uloq_inter_precision,
                                                                reported_value=uloq_inter_precision[0])
                                result.append(red(message, None, self.table_title, self.table_type))
                        else:
                            message = data_flag.get_message("uloq_inter_precision_missing",
                                                            calc_value=ap_uloq_inter_precision)
                            result.append(red(message, None, self.table_title, self.table_type))
                    else:
                        message = data_flag.get_message("uloq_inter_precision_missing",
                                                        calc_value=ap_uloq_inter_precision)
                        result.append(red(message, None, self.table_title, self.table_type))

                ap_lmh_inter_precision = ap_inter_precisions.get("lmh_precision", None)
                if ap_lmh_inter_precision is not None:
                    if self.found_inter_precisions is not None and len(self.found_inter_precisions) > 0:
                        lmh_inter_precision = self.found_inter_precisions.get("lmh_precision", None)
                        if lmh_inter_precision is not None and len(lmh_inter_precision) == 2:
                            if (round(lmh_inter_precision[0], 1) != round(ap_lmh_inter_precision[0], 1)) or (
                                    round(lmh_inter_precision[1], 1) != round(ap_lmh_inter_precision[1], 1)):
                                message = data_flag.get_message("lmh_inter_precision_mismatch",
                                                                calc_value1=ap_lmh_inter_precision[0],
                                                                calc_value2=ap_lmh_inter_precision[1],
                                                                reported_value1=lmh_inter_precision[0],
                                                                reported_value2=lmh_inter_precision[1])
                                result.append(red(message, None, self.table_title, self.table_type))
                        else:
                            message = data_flag.get_message("lmh_inter_precision_missing",
                                                            calc_value1=ap_lmh_inter_precision[0],
                                                            calc_value2=ap_lmh_inter_precision[1])
                            result.append(red(message, None, self.table_title, self.table_type))
                    else:
                        message = data_flag.get_message("lmh_inter_precision_missing",
                                                        calc_value1=ap_lmh_inter_precision[0],
                                                        calc_value2=ap_lmh_inter_precision[1])
                        result.append(red(message, None, self.table_title, self.table_type))
        except:
            pass

        try:
            if ap_intra_precisions is not None:
                ap_lloq_intra_precision = ap_intra_precisions.get("lloq_precision", None)
                if ap_lloq_intra_precision is not None:
                    if self.found_intra_precisions is not None and len(self.found_intra_precisions) > 0:
                        lloq_intra_precision = self.found_intra_precisions.get("lloq_precision", None)
                        if lloq_intra_precision is not None and len(lloq_intra_precision) == 1:
                            if round(lloq_intra_precision[0], 1) != round(ap_lloq_intra_precision[1], 1):
                                message = data_flag.get_message("lloq_intra_precision_mismatch",
                                                                calc_value=ap_lloq_intra_precision[1],
                                                                reported_value=lloq_intra_precision[0])
                                result.append(red(message, None, self.table_title, self.table_type))
                        else:
                            message = data_flag.get_message("lloq_intra_precision_missing",
                                                            calc_value=ap_lloq_intra_precision[1])
                            result.append(red(message, None, self.table_title, self.table_type))
                    else:
                        message = data_flag.get_message("lloq_intra_precision_missing",
                                                        calc_value=ap_lloq_intra_precision[1])
                        result.append(red(message, None, self.table_title, self.table_type))

                ap_uloq_intra_precision = ap_intra_precisions.get("uloq_precision", None)
                if ap_uloq_intra_precision is not None:
                    if self.found_intra_precisions is not None and len(self.found_intra_precisions) > 0:
                        uloq_intra_precision = self.found_intra_precisions.get("uloq_precision", None)
                        if uloq_intra_precision is not None and len(uloq_intra_precision) == 1:
                            if round(uloq_intra_precision[0], 1) != round(ap_uloq_intra_precision[1], 1):
                                message = data_flag.get_message("uloq_intra_precision_mismatch",
                                                                calc_value=ap_uloq_intra_precision[1],
                                                                reported_value=uloq_intra_precision[0])
                                result.append(red(message, None, self.table_title, self.table_type))
                        else:
                            message = data_flag.get_message("uloq_intra_precision_missing",
                                                            calc_value=ap_uloq_intra_precision[1])
                            result.append(red(message, None, self.table_title, self.table_type))
                    else:
                        message = data_flag.get_message("uloq_intra_precision_missing",
                                                        calc_value=ap_uloq_intra_precision[1])
                        result.append(red(message, None, self.table_title, self.table_type))

                ap_lmh_intra_precision = ap_intra_precisions.get("lmh_precision", None)
                if ap_lmh_intra_precision is not None:
                    if self.found_intra_precisions is not None and len(self.found_intra_precisions) > 0:
                        lmh_intra_precision = self.found_intra_precisions.get("lmh_precision", None)
                        if lmh_intra_precision is not None and len(lmh_intra_precision) == 2:
                            if (round(lmh_intra_precision[0], 1) != round(ap_lmh_intra_precision[0], 1)) or (
                                    round(lmh_intra_precision[1], 1) != round(ap_lmh_intra_precision[1], 1)):
                                message = data_flag.get_message("lmh_intra_precision_mismatch",
                                                                calc_value1=ap_lmh_intra_precision[0],
                                                                calc_value2=ap_lmh_intra_precision[1],
                                                                reported_value1=lmh_intra_precision[0],
                                                                reported_value2=lmh_intra_precision[1])
                                result.append(red(message, None, self.table_title, self.table_type))
                        else:
                            message = data_flag.get_message("lmh_intra_precision_missing",
                                                            calc_value1=ap_lmh_intra_precision[0],
                                                            calc_value2=ap_lmh_intra_precision[1])
                            result.append(red(message, None, self.table_title, self.table_type))
                    else:
                        message = data_flag.get_message("lmh_intra_precision_missing",
                                                        calc_value1=ap_lmh_intra_precision[0],
                                                        calc_value2=ap_lmh_intra_precision[1])
                        result.append(red(message, None, self.table_title, self.table_type))
        except:
            pass

        try:
            if ap_inter_tes is not None:
                ap_lloq_inter_te = ap_inter_tes.get("lloq_te", None)
                if ap_lloq_inter_te is not None:
                    if self.found_inter_tes is not None and len(self.found_inter_tes) > 0:
                        lloq_inter_te = self.found_inter_tes.get("lloq_te", None)
                        if lloq_inter_te is not None and len(lloq_inter_te) == 1:
                            if round(lloq_inter_te[0], 1) != round(ap_lloq_inter_te, 1):
                                message = data_flag.get_message("lloq_inter_te_mismatch", calc_value=ap_lloq_inter_te,
                                                                reported_value=lloq_inter_te[0])
                                result.append(red(message, None, self.table_title, self.table_type))
                        else:
                            message = data_flag.get_message("lloq_inter_te_missing", calc_value=ap_lloq_inter_te)
                            result.append(red(message, None, self.table_title, self.table_type))
                    else:
                        message = data_flag.get_message("lloq_inter_te_missing", calc_value=ap_lloq_inter_te)
                        result.append(red(message, None, self.table_title, self.table_type))

                ap_uloq_inter_te = ap_inter_tes.get("uloq_te", None)
                if ap_uloq_inter_te is not None:
                    if self.found_inter_tes is not None and len(self.found_inter_tes) > 0:
                        uloq_inter_te = self.found_inter_tes.get("uloq_te", None)
                        if uloq_inter_te is not None and len(uloq_inter_te) == 1:
                            if round(uloq_inter_te[0], 1) != round(ap_uloq_inter_te, 1):
                                message = data_flag.get_message("uloq_inter_te_mismatch", calc_value=ap_uloq_inter_te,
                                                                reported_value=uloq_inter_te[0])
                                result.append(red(message, None, self.table_title, self.table_type))
                        else:
                            message = data_flag.get_message("uloq_inter_te_missing", calc_value=ap_uloq_inter_te)
                            result.append(red(message, None, self.table_title, self.table_type))
                    else:
                        message = data_flag.get_message("uloq_inter_te_missing", calc_value=ap_uloq_inter_te)
                        result.append(red(message, None, self.table_title, self.table_type))

                ap_lmh_inter_te = ap_inter_tes.get("lmh_te", None)
                if ap_lmh_inter_te is not None:
                    if self.found_inter_tes is not None and len(self.found_inter_tes) > 0:
                        lmh_inter_te = self.found_inter_tes.get("lmh_te", None)
                        if lmh_inter_te is not None and len(lmh_inter_te) == 2:
                            if (round(lmh_inter_te[0], 1) != round(ap_lmh_inter_te[0], 1)) or (
                                    round(lmh_inter_te[1], 1) != round(ap_lmh_inter_te[1], 1)):
                                message = data_flag.get_message("lmh_inter_te_mismatch", calc_value1=ap_lmh_inter_te[0],
                                                                calc_value2=ap_lmh_inter_te[1],
                                                                reported_value1=lmh_inter_te[0],
                                                                reported_value2=lmh_inter_te[1])
                                result.append(red(message, None, self.table_title, self.table_type))
                        else:
                            message = data_flag.get_message("lmh_inter_te_missing", calc_value1=ap_lmh_inter_te[0],
                                                            calc_value2=ap_lmh_inter_te[1])
                            result.append(red(message, None, self.table_title, self.table_type))
                    else:
                        message = data_flag.get_message("lmh_inter_te_missing", calc_value1=ap_lmh_inter_te[0],
                                                        calc_value2=ap_lmh_inter_te[1])
                        result.append(red(message, None, self.table_title, self.table_type))
        except KeyError:
            pass

        if len(result) == 0:
            result.append(
                green('Values reported in Validation Summary confirmed', None, self.table_title, self.table_type))

        self.result += result


class Parallelism(Table, Flags):
    def __init__(self, parsed_table, analytes, template_type):
        Table.__init__(self, parsed_table, template_type)
        self.analyte = get_analyte(analytes, self.tb_title)

    def validate_table_format(self):
        error_messages = self.error_messages
        required_col = {"Run ID", "Sample ID", "Nominal Concentration (ng/mL)", "Measured Concentration (ng/mL)", "%RE"}
        missing_col = required_col.intersection(self.missing_col)
        if missing_col:
            message = error_messages.get_message("missing_col", col_names=", ".join(missing_col))
            self.result.append(red(message, None, self.table_title, self.table_type))
            self.valid_format = False

    def process_table(self):
        table = self.table_df
        try:
            table["nominal"] = table["nominal"].apply(parse_decimal)
            table["conc"] = table["conc"].apply(parse_decimal)
            table["re"] = table["re"].apply(parse_decimal)

            table["calc_re"] = ((table["conc"] - table["nominal"]) / table["nominal"]) * 100
            table["per_diff_re"] = calculate_per_diff(table["re"], table["calc_re"])
        except:
            message = self.error_messages.get_message("data_error")
            self.result.append(red(message, None, self.table_title, self.table_type))
            self.valid_data = False
        table = table.fillna("")
        self.data_df = fill_val(table)

    @staticmethod
    def get_error_dict():
        data_flag = FlagProperties("parallelism")
        return data_flag

    def validate_table_data(self):
        valid_count_diff = Decimal(self.valid_difference_values.get_message("difference_n"))
        valid_re_cv_diff = Decimal(self.valid_difference_values.get_message("difference_re_cv"))
        valid_mean_diff = Decimal(self.valid_difference_values.get_message("difference_mean"))

        result = []
        data_flag = self.get_error_dict()
        table = self.data_df
        samples = table["samples"].to_list()
        batches = table["run_id"].to_list()
        reported_res = table["re"]
        calc_res = table["calc_re"]
        per_diff_res = table["per_diff_re"]
        test_threshold = Decimal(self.threshold_values.get_message("re_cv_threshold"))
        if "lm" in self.analysis_type:
            test_threshold += Decimal(5)
        for index, sample in enumerate(samples):
            batch = batches[index]

            reported_re = reported_res[index]
            calc_re = calc_res[index]
            per_diff_re = per_diff_res[index]
            re_error_para = {"sample": sample, "batch": batch, "reported_value": reported_re,
                             "calc_value": utils.format_value(calc_re)}
            if str(per_diff_re) != "":
                if per_diff_re > valid_re_cv_diff:
                    if abs(calc_re) <= Decimal(1) and abs(reported_re) <= Decimal(1):
                        message = data_flag.get_message("re_rounding", **re_error_para)
                        result.append(yellow(message, None, self.table_title, self.table_type))
                    else:
                        message = data_flag.get_message("re_error", **re_error_para)
                        result.append(red(message, None, self.table_title, self.table_type))

            red_threshold_para = {"sample": sample, "batch": batch, "threshold": test_threshold}
            yellow_threshold_para = {"sample": sample, "batch": batch, "threshold": test_threshold - Decimal(5)}
            if str(calc_re).strip() != "":
                if abs(calc_re) > test_threshold:
                    message = data_flag.get_message("re", **red_threshold_para)
                    result.append(red(message, None, self.table_title, self.table_type))

                elif (test_threshold - Decimal(5)) <= abs(calc_re) <= test_threshold:
                    message = data_flag.get_message("re", **yellow_threshold_para)
                    result.append(yellow(message, None, self.table_title, self.table_type))

        if len(result) == 0:
            result.append(green(f"All values within acceptable range", None, self.table_title, self.table_type))
        self.result += result
