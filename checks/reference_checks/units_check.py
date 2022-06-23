from checks.check_result import red, green, yellow


def cross_check_unit(conc_units, all_conc_units, multiple_analyte, analytes):
    result = []
    for analyte in analytes:
        if multiple_analyte:
            units = conc_units.get(analyte, "")
            units = str(units).replace("(", "").replace(")", "").strip()
            all_units = all_conc_units.get(analyte, [])
            unit_mismatch_msg = f"Unit of (ng/mL) found in calibration curve table for analyte {analyte}. Possible unit reporting error in other tables for the same analyte"
            unit_confirm_msg = f"Unit of ({units}) confirmed for analyte {analyte} in all tables"
            unit_missing_msg = f"Calibration curve table not found for analyte {analyte}, unable to confirm if concentration unit was consistent in all tables for the same analyte"
        else:
            units = conc_units
            units = str(units).replace("(", "").replace(")", "").strip()
            all_units = all_conc_units
            unit_mismatch_msg = f"Unit of ({units}) found in calibration curve table. Possible unit reporting error in other tables"
            unit_confirm_msg = f"Unit of ({units}) confirmed in all tables"
            unit_missing_msg = "Calibration curve table not found, unable to confirm if concentration unit was consistent in all tables"

        all_units = [str(x).replace("(", "").replace(")", "").strip() for x in all_units if str(x).strip() != ""]
        if str(units).strip() == "":
            result.append(yellow(unit_missing_msg, None, "Entire Report", "Entire Report"))
        elif all_units:
            if len(set(all_units)) == 1 and str(all_units[0]).lower().strip() == str(units).lower().strip():
                result.append(green(unit_confirm_msg, None, "Entire Report", "Entire Report"))
            else:
                result.append(red(unit_mismatch_msg, None, "Entire Report", "Entire Report"))

    return result

