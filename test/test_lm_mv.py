import json

import boto3
import pytest

from report_analyzer import analyze_report

ddb_c = boto3.client("dynamodb")
s3_c = boto3.client("s3")


@pytest.fixture
def test_dataset_1():
    event = {
        "user_id": "mmcguigan",
        "analysis_id": "1636732566240",
        "analysis_det": {
            "analysisType": "LMR",
            "analysisSubtype": "MV0",
            "projectCodeSA": "LM 2in1 Val Report Test 1",
            "projectCodeMV": "null",
            "regressionModel": "5PL",
            "lloq": "null",
            "uloq": "null",
            "re_type": "true",
            "re_value": "null",
            "lts80_dialog": "99",
            "lts20_dialog": "null",
            "extraction_dialog": "null",
            "benchtop_dialog": "24",
            "freeze_thaw_dialog": "4",
            "multipleAnalytes": "True",
            "first_patient_dialog": "null",
            "sampleMatrix": "Serum",
            "analysisSubtype_ADA": "null",
            "antiCoagulant": "",
            "species": "Human",
            "samplePrep": "null",
            "sampleDate": "2021-02-05T05:00:00.000Z",
            "analyteNames": [
                "COMPOUND-X"
            ],
            "userId": "mmcguigan",
            "analysisId": "1636732566240",
            "analytes": {
                "xlsx": [],
                "docx": [
                    [
                        "mmcguigan/us-east-2:2a03920e-7333-461c-87b2-47f55ef222f9/1636732566240/COMPOUND-X/Report/LM PK Val two-Analyte Rpt_2_27Oct2021_MM.docx",
                        "/tmp/LM PK Val two-Analyte Rpt_2_27Oct2021_MM.docx",
                        "Method Validation Report"
                    ]
                ],
                "pdf": []
            },
            "analysisStatus": "TableSelected"
        },
        "file_type": "docx",
        "AnalysisType": "LMR",
        "bucket_name": "kcas-dev-upload",
        "table_name": "unit-test",
        "dash_table": "kcas-dev-Dash"
    }
    return event


def test_case_3(test_dataset_1):
    user_id = test_dataset_1["user_id"]
    analysis_id = test_dataset_1["analysis_id"]
    table_name = test_dataset_1["table_name"]
    location_dict = test_dataset_1['analysis_det']['analytes']

    reference_item = ddb_c.get_item(TableName=table_name, Key={
        'UserId': {'S': user_id},
        'AnalysisId': {'S': analysis_id}
    })['Item']
    reference_result = json.loads(reference_item["AnalysisResult"]["S"])

    overall_result = analyze_report({
        'analytes': test_dataset_1['analysis_det']['analyteNames'],
        'file_name': location_dict,
        'bucket_name': test_dataset_1['bucket_name'],
        'AnalysisDetails': json.loads(reference_item['AnalysisDetails']['S']),
        'SelectedTables': reference_item['SelectedTables']['S'],  # Request Modifications
        "species": test_dataset_1['analysis_det']['species'],
        "matrix": test_dataset_1['analysis_det']['sampleMatrix'],
        "antiCoagulant": test_dataset_1['analysis_det']['antiCoagulant']
    })
    process_result = overall_result[0]
    for r_result, p_result in zip(reference_result, process_result):
        for r_result_item, p_result_item in zip(r_result["resultItem"], p_result["resultItem"]):
            for r_item, p_item in zip(r_result_item['item'], p_result_item['item']):
                if len(set(r_item["message"]).difference(set(p_item["message"]))) != 0:
                    assert r_item["message"] == p_item["message"], f"Error in test case 1 LM MV for table type " \
                                                                   f"{r_result['heading']}"
                assert True

