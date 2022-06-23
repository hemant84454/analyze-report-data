from collections import defaultdict
from checks.check_result import red, yellow, green
from checks.utils import get_analyte


def take_table_type(table_dict):
    return table_dict["table_type"], table_dict.get("table_subtype", "")


def check_multi_analyte_across_title(tables: list, analytes: list):
    result = []
    table_count = 0
    unnecessary_table_list = ["Overall and Individual Batch Performance", "Samples Analyzed"]
    black_list_tables = ["Expiration Information"]
    message_group = "Analytical Note"
    tables.sort(key=take_table_type)
    analyte_dict = defaultdict(int)
    previous_table_type = tables[0]["table_type"]
    previous_table_subtype = tables[0].get("table_subtype", "")

    for table in tables:
        table_type = table['table_type']
        if table_type in black_list_tables:
            continue
        table_subtype = table.get("table_subtype", "")
        table_title = table["table_title"] + table['tb_title']
        if previous_table_subtype != table_subtype:
            title_count = sum(analyte_dict.values())
            if title_count == 0:
                message = f"Cannot find analyte name {', '.join(analytes)} in table title"
                if previous_table_type in unnecessary_table_list:
                    result.append(yellow(message, None, previous_table_type, message_group))
                else:
                    result.append(red(message, None, previous_table_type, message_group))
            else:
                for analyte in analytes:
                    count = analyte_dict.get(analyte, 0)
                    if count > 1:
                        message = f"Multiple occurrences detected, more than one table of this type found with " \
                                  f"analyte {analyte} in table title"
                        result.append(yellow(message, None, previous_table_type, message_group))
                    elif count == 0 and table_count != title_count:
                        message = f"Cannot find analyte name {analyte} in table title"
                        result.append(red(message, None, previous_table_type, message_group))
            analyte_dict = defaultdict(int)
            previous_table_type = table_type
            previous_table_subtype = table_subtype
            table_count = 0

        analyte = get_analyte(analytes, table_title)
        table_count += 1
        if analyte is not None:
            analyte_dict[analyte] += 1

    if not result:
        result.append(green("No multiple occurrences of tables detected for any given analyte", None, message_group,
                            message_group))

    return result
