from jinja2 import Environment, FileSystemLoader
from botocore.exceptions import ClientError
from horizontal_chart import horizontalChart
from report_analyzer import analyze_report
from datetime import datetime
import weasyprint
import traceback
import logging
import boto3
import json
import os

UTC_OFFSET = 4


env = Environment(loader=FileSystemLoader('.'))
template = env.get_template("export.html")

null = None
logger = logging.getLogger()
logger.setLevel(logging.INFO)
s3_c = boto3.client('s3')
s3_r = boto3.resource('s3')
ddb_c = boto3.client('dynamodb', region_name='us-east-2')
dynamodb = boto3.resource('dynamodb', region_name='us-east-2')
snsclient = boto3.client('sns', region_name='us-east-2')
snsArn = os.environ.get("snsARN", None)
domain = os.environ.get("domain", None)
company = os.environ.get("company", "")


def send_notification(sns_arn, analysis_id, message):
    message = f"Backend Runtime Error in Domain {domain} for AnalysisId: {analysis_id} \n\n" \
              f"###############################################################\n"\
              f"Error: \n"\
              f" {message}\n" \
              f"############################################"
    try:
        snsclient.publish(
            TargetArn=sns_arn,
            Subject=f'Execution error for Analysis Id {analysis_id}',
            Message=message
        )
    except ClientError as e:
        logger.error("An error occurred: %s" % e)


def lambda_handler(event, context):
    """Lambda Function to call Step functions depending on AnalysisType attribute in DynamoDBTable

        Parameters
        ----------
        event: dict, required
            DynamoDB NewImage type Stream Input Format
            Event doc: https://docs.aws.amazon.com/lambda/latest/dg/with-ddb.html

        context: object, required
            Lambda Context runtime methods and attributes
            Context doc: https://docs.aws.amazon.com/lambda/latest/dg/python-context-object.html

        Returns
        ------
        API Gateway Lambda Proxy Output Format: dict
            Return doc: https://docs.aws.amazon.com/apigateway/latest/developerguide/set-up-lambda-proxy-integrations
            .html
        """
    logger.info('EVENT:\n' + json.dumps(event))
    table_name = event['dash_table']
    table = dynamodb.Table(table_name)
    user_id = event['user_id']
    analysis_id = event['analysis_id']
    table_name = event['table_name']
    location_dict = event['analysis_det']['analytes']
    [s3_c.download_file(event['bucket_name'], file[0], file[1]) for file in location_dict['xlsx']]

    try:
        s3_c.download_file(event['bucket_name'], 'dropdowns.json', '/tmp/dropdowns.json')
        s3_c.download_file(event['bucket_name'], 'ObjectNotation.json', '/tmp/ObjectNotation.json')

    except Exception as e:
        tb = traceback.format_exc()
        logger.error('Core Fallacy %s', e)
        if snsArn:
            send_notification(snsArn, analysis_id, tb)
        raise e

    item = ddb_c.get_item(TableName=table_name, Key={
        'UserId': {'S': user_id},
        'AnalysisId': {'S': analysis_id}
    })['Item']

    try:
        if location_dict['pdf']:
            raw_path = '/'.join(str(location_dict['pdf'][0][0]).split('/')[:-1])
            selected_csv_tables = json.loads(item['SelectedTables']['S'])
            for csv_table in selected_csv_tables:
                table_index = csv_table['table_index']
                table_file = str(table_index.split('/')[-1])
                s3_csv_location = raw_path + table_index
                s3_c.download_file(event['bucket_name'], s3_csv_location, '/tmp/' + table_file)
            logger.info(f"Selected tables {selected_csv_tables}")

        overall_result = analyze_report({
            'analytes': event['analysis_det']['analyteNames'],
            'file_name': location_dict,
            'bucket_name': event['bucket_name'],
            'AnalysisDetails': json.loads(item['AnalysisDetails']['S']),
            'SelectedTables': item['SelectedTables']['S'],  # Request Modifications
            "species": event['analysis_det']['species'],
            "matrix": event['analysis_det']['sampleMatrix'],
            "antiCoagulant": event['analysis_det']['antiCoagulant']
        })

        # Generate Docx visualization summary
        bucket_name = event['bucket_name']
        analyte_list = event['analysis_det']['analytes']
        file_list = []
        if analyte_list['docx']:
            file_list = analyte_list['docx']
        elif analyte_list['xlsx']:
            file_list = analyte_list['xlsx']
        elif analyte_list['pdf']:
            file_list = analyte_list['pdf']

        s3_output_pdf = ''
        s3_output_xlsx = ''
        file_name = ''

        if file_list:
            file = file_list[0]
            file_name = str(file[0].split('/')[-1])
            s3_output_pdf = '/'.join(file[0].split('/')[:-2]) + '/output/' + file_name.split(".")[
                0] + "_output.pdf"
            s3_output_xlsx = '/'.join(file[0].split('/')[:-2]) + '/output/' + file_name.split(".")[
                0] + "_output.xlsx"

        output_file = '/tmp/' + file_name.split(".")[0] + "_output.pdf"
        analysis_type = str(event['analysis_det']['analysisType']) + '_' + str(event['analysis_det']['analysisSubtype'])
        date = datetime.now().strftime("%a %m/%d/%Y, %H:%M:%S")
        doc_date = datetime.now().strftime("%a %m/%d/%Y")
        try:
            if analysis_type is not None:
                analysis_type = analysis_type.lower()
        except Exception as e:
            analysis_type = ''
            print(e)

        if 'lm' in analysis_type:
            first_part = 'Large Molecule'
        elif 'sm' in analysis_type:
            first_part = 'Small Molecule'
        else:
            first_part = ''

        if 'mv' in analysis_type:
            second_part = ' Method Validation Findings'
        elif 'sa' in analysis_type:
            second_part = ' Sample Analysis Findings'
        else:
            second_part = ''
        doc_analysis_type = first_part + second_part
        template_vars = {
            "date": doc_date,
            "analysis_type": doc_analysis_type,
            "data": overall_result[0],
            "file_name": file_name
        }
        excel_file = overall_result[2]
        horizontalbar_chart = horizontalChart(overall_result[0], analysis_id, file_name, first_part, user_id)
        html_out = template.render(template_vars)
        weasyprint.HTML(string=html_out).write_pdf(output_file, stylesheets=["style.css"])

        if s3_output_pdf:
            s3_c.upload_file(output_file, bucket_name, s3_output_pdf)
            s3_r.Object(bucket_name, s3_output_xlsx).put(Body=excel_file, Tagging=f"company={company}")
            excel_file.close()

        analysis_details = event['analysis_det']
        summary_statistics = {
            'analysisId': analysis_id,
            'userId': user_id,
            'overall_g': overall_result[1]['green'],
            'overall_y': overall_result[1]['yellow'],
            'overall_r': overall_result[1]['red'],
            'date': date,
            'analytes': analysis_details['analyteNames'],
            'projectCodeSA': analysis_details['projectCodeSA'],
            'projectCodeMV': analysis_details['projectCodeMV'],
            'analysisType': analysis_details['analysisType'],
            'analysisSubtype': analysis_details['analysisSubtype']
        }
        item['SummaryStatistics'] = {'S': json.dumps(summary_statistics)}
        item['AnalysisResult'] = {'S': json.dumps(overall_result[0])}
        item['AnalysisStatus'] = {'S': 'Complete'}
        item['AnalysisResultFile'] = {'S': s3_output_pdf}
        item['AnalysisResultXlsxFile'] = {'S': s3_output_xlsx}

        ddb_c.put_item(TableName=table_name, Item=item)
        table.put_item(Item=horizontalbar_chart)
        return event
    except Exception as e:
        tb = traceback.format_exc()
        item['AnalysisStatus'] = {'S': 'Error'}
        item['AnalysisResult'] = {'S': 'Error occurred while processing data '}
        item['Error'] = {'S': tb}
        if snsArn:
            send_notification(snsArn, analysis_id, tb)
        ddb_c.put_item(TableName=table_name, Item=item)
        logger.error(f'Result wrote in DynamoDB, Error is {tb}')
        raise e
    finally:
        # add code to delete file
        pass


if __name__ == "__main__":
    with open("event.json") as f:
        data = json.load(f)
    lambda_handler(data, "context")
