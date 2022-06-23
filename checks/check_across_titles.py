from checks.check_result import green, yellow, red
from checks.utils import stability_name
from decimal import Decimal


def check_across_titles(tables, possible_anticoagulants, possible_sample_matrices, possible_sample_preps,
                        possible_sample_species, options):
    result = []
    match_all = True
    yellow_flags_table_list = ['Expiration Information', 'Overall and Individual Batch Performance', 'Samples Analyzed']
    analyte_array = set()
    analyte_name = options.get("analyte_name", None)
    if analyte_name is not None:
        analyte_array = set(analyte_name.lower().split(", "))
    user_sample_matrix = options.get('sample_matrix', None)
    if user_sample_matrix is not None:
        user_sample_matrix = user_sample_matrix.strip()
    user_anticoagulant = options.get('anticoagulant', None)
    if user_anticoagulant is not None:
        user_anticoagulant = user_anticoagulant.strip()
    user_species = options.get('species', None)
    if user_species is not None:
        user_species = user_species.strip()
    user_sample_prep = options.get('samplePrep', None)
    if user_sample_prep is not None:
        user_sample_prep = user_sample_prep.strip()
    user_benchtop_hours = options.get('benchtop_stability', None)
    if user_benchtop_hours is not None:
        user_benchtop_hours = str(user_benchtop_hours).strip()
    user_extraction_hours = options.get('extraction_stability', None)
    if user_extraction_hours is not None:
        user_extraction_hours = str(user_extraction_hours).strip()
    lts_days = options.get('lts80', None)
    if lts_days is not None:
        lts_days = str(lts_days).strip()
    for table in tables:
        found_analyte = False
        try:
            tb_title = f"{str(table['tb_title'])} {str(table['table_title'])}"
            table_title = str(table['table_title'])
            table_type = str(table["table_type"])

            if table_type != "Blood Stability":
                title_array = tb_title.lower().split()
                title_array = [str(x).replace("(", "").replace(")", "") for x in title_array]

                stability = stability_name(title_array)

                if stability is not None:
                    if "extract" in str(stability).lower() or "short" in str(stability).lower() or "bench" in \
                            str(stability).lower():
                        try:
                            hr_index = title_array.index("hours")
                            hour = title_array[(hr_index - 1)]
                            if "bench" in str(stability).lower() or "short" in str(stability).lower():
                                if Decimal(user_benchtop_hours) == Decimal(hour):
                                    pass
                                else:
                                    result.append(
                                        red(
                                            f'Contains the {stability} (hours) {hour} hours instead of '
                                            f'{user_benchtop_hours} hours', None, table_title, table_type))
                            else:
                                if Decimal(user_extraction_hours) == Decimal(hour):
                                    pass
                                else:
                                    result.append(
                                        red(
                                            f'Contains the {stability} (hours) {hour} hours instead of {user_extraction_hours} hours',
                                            None, table_title,
                                            table_type))
                        except IndexError:
                            pass

                    elif "long" in str(stability).lower():
                        try:
                            try:
                                day_index = title_array.index("day")
                            except ValueError:
                                day_index = title_array.index("days")
                            day = title_array[(day_index - 1)]
                            if Decimal(lts_days) == Decimal(day):
                                pass
                            else:
                                result.append(
                                    red(f'Contains the {stability} (days) {day} days instead of {lts_days} days', None,
                                        table_title,
                                        table_type))
                        except IndexError:
                            pass

                for species in possible_sample_species:
                    if species.lower() in user_species.lower():
                        pass
                    elif species.lower() in title_array:
                        match_all = False
                        result.append(
                            red('Contains the species ' + species + ' in table title instead of ' + user_species, None, table_title,
                                table_type))

                for sample_matrix in possible_sample_matrices:
                    if sample_matrix.lower() in user_sample_matrix.lower():
                        pass
                    elif sample_matrix.lower() in title_array:
                        match_all = False
                        result.append(
                            red('Contains the sample matrix ' + sample_matrix + ' in table title instead of ' +
                                user_sample_matrix, None, table_title, table_type))

                common_words = analyte_array.intersection(set(title_array))
                if len(common_words):
                    found_analyte = True
                elif analyte_name.lower().strip() in tb_title.lower():
                    found_analyte = True

                if not found_analyte:
                    if table_type in yellow_flags_table_list:
                        result.append(yellow(f"Can not find Analyte name {analyte_name} in table title", None,
                                             table_title, table_type))
                    else:
                        result.append(red(f"Can not find Analyte name {analyte_name} in table title", None, table_title,
                                          table_type))

                for anticoagulant in possible_anticoagulants:
                    if anticoagulant.lower() in user_anticoagulant.lower():
                        pass
                    elif anticoagulant.lower() in tb_title.lower():
                        match_all = False
                        result.append(
                            red('Contains the anticoagulant ' + str(anticoagulant) + ' in table title instead of ' +
                                str(user_anticoagulant), None,
                                table_title, table_type))
                        break

                for sample_preps in possible_sample_preps:
                    if sample_preps.lower() in user_sample_prep.lower():
                        pass
                    elif sample_preps.lower() in title_array:
                        result.append(
                            red('Contains the sample preparation ' + sample_preps + ' in table title instead of '
                                + user_sample_prep, None,
                                table_title, table_type))
        except:
            continue

    if match_all:
        result.append(green("No table titles contain an improper species, sample matrix, or anticoagulant", None,
                            "Entire Report", "Entire Report"))

    return result
