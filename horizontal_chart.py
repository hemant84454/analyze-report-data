from datetime import datetime


def horizontalChart(overall_result, analysis_id, file_name, analysis_type, userid):
    horizontalbarchart = dict()
    red = []
    yellow = []
    green = []
    horizontalbarchart["ChartLabels"] = []
    overall_result = overall_result[1:]
    date = datetime.now().strftime("%b %d, %Y")
    for result in overall_result:
        pre_green = []
        pre_red = []
        pre_yellow = []
        horizontalbarchart["ChartLabels"].append(result["heading"])
        for resultitem in result["resultItem"]:
            for items in resultitem["item"]:
                if items["color"].lower() == "green":
                    pre_green.append(len(items["message"]))
                if items["color"].lower() == "red":
                    pre_red.append(len(items["message"]))
                if items["color"].lower() == "yellow":
                    pre_yellow.append(len(items["message"]))

        green.append(sum(pre_green))
        yellow.append(sum(pre_yellow))
        red.append(sum(pre_red))
    horizontalbarchart["ChartData"] = [{
        "data": red,
        "label": "Error",
        "stack": "a"
    }, {
        "data": yellow,
        "label": "Warning",
        "stack": "a"
    }, {
        "data": green,
        "label": "Information",
        "stack": "a"
    }]
    horizontalbarchart["AnalysisId"] = analysis_id
    horizontalbarchart["ChartType"] = "horizontalBar"
    horizontalbarchart["ChartId"] = file_name.split(".")[0]
    horizontalbarchart["AnalysisDate"] = date
    horizontalbarchart["AnalysisType"] = analysis_type
    horizontalbarchart["UserId"] = userid

    return horizontalbarchart