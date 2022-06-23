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
        "analysis_id": "1636724540461",
        "analysis_det": {
            "analysisType": "SMR",
            "analysisSubtype": "SA0",
            "projectCodeSA": "SM PK 3in1 report test 4",
            "projectCodeMV": "null",
            "regressionModel": "LNR",
            "lloq": "null",
            "uloq": "null",
            "re_type": "true",
            "re_value": "null",
            "lts80_dialog": "90",
            "lts20_dialog": "null",
            "extraction_dialog": "98",
            "benchtop_dialog": "null",
            "freeze_thaw_dialog": "3",
            "multipleAnalytes": "True",
            "first_patient_dialog": "null",
            "sampleMatrix": "Serum",
            "analysisSubtype_ADA": "null",
            "antiCoagulant": "",
            "species": "Rat",
            "samplePrep": "null",
            "sampleDate": "2020-06-18T04:00:00.000Z",
            "analyteNames": [
                "Topiramate"
            ],
            "userId": "mmcguigan",
            "analysisId": "1636724540461",
            "analytes": {
                "xlsx": [],
                "docx": [
                    [
                        "mmcguigan/us-east-2:2a03920e-7333-461c-87b2-47f55ef222f9/1636724540461/Topiramate/Report/SM Mock PK GLP Assay Report_3in1_17Jun2021_C.docx",
                        "/tmp/SM Mock PK GLP Assay Report_3in1_17Jun2021_C.docx",
                        "Sample Analysis Report"
                    ]
                ],
                "pdf": []
            },
            "analysisStatus": "TableSelected"
        },
        "file_type": "docx",
        "AnalysisType": "SMR",
        "bucket_name": "kcas-dev-upload",
        "table_name": "unit-test",
        "dash_table": "kcas-dev-Dash"
    }
    return event


def test_case_1(test_dataset_1):
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
                    assert r_item["message"] == p_item["message"], f"Error in test case 1 SM SA for table type " \
                                                                   f"{r_result['heading']}"
                assert True
                

@pytest.fixture
def test_dataset_2():
    event = {
        "user_id": "rziegler",
        "analysis_id": "1637163028127",
        "analysis_det": {
            "analysisType": "SMR",
            "analysisSubtype": "SA0",
            "projectCodeSA": "Test",
            "projectCodeMV": "null",
            "regressionModel": "LNR",
            "lloq": "null",
            "uloq": "null",
            "re_type": "true",
            "re_value": "null",
            "lts80_dialog": "90",
            "lts20_dialog": "null",
            "extraction_dialog": "98",
            "benchtop_dialog": "null",
            "freeze_thaw_dialog": "3",
            "multipleAnalytes": "True",
            "first_patient_dialog": "null",
            "sampleMatrix": "Serum",
            "analysisSubtype_ADA": "null",
            "antiCoagulant": "",
            "species": "Rat",
            "samplePrep": "null",
            "sampleDate": "2020-06-18T04:00:00.000Z",
            "analyteNames": [
                "Topiramate"
            ],
            "userId": "rziegler",
            "analysisId": "1637163028127",
            "analytes": {
                "xlsx": [],
                "docx": [
                    [
                        "rziegler/us-east-2:6aa087cd-f94d-46a7-b3b5-ee64a3b2cc7d/1637163028127/Topiramate/Report/SM Mock PK GLP Assay Report_3in1_17Jun2021_C.docx",
                        "/tmp/SM Mock PK GLP Assay Report_3in1_17Jun2021_C.docx",
                        "Sample Analysis Report"
                    ]
                ],
                "pdf": []
            },
            "analysisStatus": "TableSelected"
        },
        "file_type": "docx",
        "AnalysisType": "SMR",
        "bucket_name": "kcas-dev-upload",
        "table_name": "unit-test",
        "dash_table": "kcas-dev-Dash"
    }
    return event


def test_case_6(test_dataset_2):
    user_id = test_dataset_2["user_id"]
    analysis_id = test_dataset_2["analysis_id"]
    table_name = test_dataset_2["table_name"]
    location_dict = test_dataset_2['analysis_det']['analytes']

    reference_item = ddb_c.get_item(TableName=table_name, Key={
        'UserId': {'S': user_id},
        'AnalysisId': {'S': analysis_id}
    })['Item']
    reference_result = json.loads(reference_item["AnalysisResult"]["S"])

    overall_result = analyze_report({
        'analytes': test_dataset_2['analysis_det']['analyteNames'],
        'file_name': location_dict,
        'bucket_name': test_dataset_2['bucket_name'],
        'AnalysisDetails': json.loads(reference_item['AnalysisDetails']['S']),
        'SelectedTables': reference_item['SelectedTables']['S'],  # Request Modifications
        "species": test_dataset_2['analysis_det']['species'],
        "matrix": test_dataset_2['analysis_det']['sampleMatrix'],
        "antiCoagulant": test_dataset_2['analysis_det']['antiCoagulant']
    })
    process_result = overall_result[0]
    for r_result, p_result in zip(reference_result, process_result):
        for r_result_item, p_result_item in zip(r_result["resultItem"], p_result["resultItem"]):
            for r_item, p_item in zip(r_result_item['item'], p_result_item['item']):
                if len(set(r_item["message"]).difference(set(p_item["message"]))) != 0:
                    assert r_item["message"] == p_item["message"], f"Error in test case 2 SM SA for table type " \
                                                                   f"{r_result['heading']}"
                assert True