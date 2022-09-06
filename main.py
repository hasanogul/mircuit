"""
This script runs the application using a development server.
It contains the definition of routes and views for the application.
"""

from flask import Flask, render_template, request, redirect, url_for
from flask import send_file
from flask import current_app
from flask import send_from_directory
from werkzeug.utils import secure_filename
import pandas as pd 
import math as mt
app = Flask(__name__)



@app.route('/')
def index():
    return render_template('mainPage.html')

def difexp1(file_name_1, file_name_2, output_1, output_2):
    import pandas as pd 
    import math as mt
    folder_1 = pd.read_csv(file_name_1, sep = "\t")
    folder_2 = pd.read_csv(file_name_2, sep = "\t")
    for i, j in folder_1.iterrows():
        folder_1.loc[i, "log(FC)"] = mt.log((j["Tumor"]/j["Normal"]), 2)
    for k, l in folder_2.iterrows():
        folder_2.loc[k, "log(FC)"] = mt.log((l["Tumor"]/l["Normal"]), 2)

    thrs = 1.5                                                             
    folder_DE_1 = folder_1[abs(folder_1["log(FC)"]) > thrs] 

    thrs = 1.5 
    folder_DE_2 = folder_2[abs(folder_2["log(FC)"]) > thrs]

    folder_DE_1_sorted = folder_DE_1.sort_values("log(FC)", axis = 0, ascending = False)
    for a, b in folder_DE_1_sorted.iterrows():
        if b["log(FC)"] >= 0:
            folder_DE_1_sorted.loc[a, "Expression"] = "Up"
        else: 
            folder_DE_1_sorted.loc[a, "Expression"] = "Down"

    folder_DE_2_sorted = folder_DE_2.sort_values("log(FC)", axis = 0, ascending = False)
    for c, d in folder_DE_2_sorted.iterrows():
        if d["log(FC)"] >= 0:
            folder_DE_2_sorted.loc[c, "Expression"] = "Up"
        else: 
            folder_DE_2_sorted.loc[c, "Expression"] = "Down"

    folder_DE_1_sorted.to_csv(output_1, sep = "\t")
    folder_DE_2_sorted.to_csv(output_2, sep = "\t")

    with open("Other_output/merge.txt", "w") as f1:
        f1.write("mRNA\tExpression\tExpression\tmiRNA\n")
        for m, n in folder_DE_1_sorted.iterrows():
            for o, p in folder_DE_2_sorted.iterrows():
                if (n["Expression"] == "Up") & (p["Expression"] == "Down") | (n["Expression"] == "Down") & (p["Expression"] == "Up"):
                    f1.write("{0}\t{1}\t{2}\t{3}\n".format(n["mRNA"], n["Expression"], p["Expression"], p["miRNA"]))
        f1.close()
    
    merge = ((pd.read_csv("Other_output/merge.txt", sep = "\t").drop_duplicates(keep = "first")).loc[:, ["mRNA", "miRNA"]]).to_csv("Other_output/merge_edited.csv", sep = "\t")
    
    TargetScan = pd.read_csv("Edited_data/TargetScan.csv", sep = "\t")
    merge = pd.read_csv("Other_output/merge_edited.csv", sep = "\t")
    miRNA_target = (pd.merge(TargetScan, merge, how = "inner").drop_duplicates(keep = "first").loc[:, ["miRNA", "mRNA"]]).to_csv("Other_output/miRNA_target(mRNA).csv", sep = "\t")

    miRNA_target = pd.read_csv("Other_output/miRNA_target(mRNA).csv", sep = "\t").loc[:, ["miRNA", "mRNA"]]
    miRNA_target.columns = ["miRNA", "TF"]
    miRNA_target.to_csv("Other_output/miRNA_target(TF).csv", sep = "\t")

    hTFtarget = pd.read_csv("Edited_data/hTFtarget.csv", sep = "\t").loc[:, ["TF"]]
    miRNA_target_TF = pd.read_csv("Other_output/miRNA_target(TF).csv", sep = "\t").loc[:, ["miRNA", "TF"]]
    miRNA_TF = ((pd.merge(hTFtarget, miRNA_target_TF, how = "inner").loc[:, ["miRNA", "TF"]]).drop_duplicates(keep = "first")).to_csv("Other_output/miRNA_TF.csv", sep = "\t")

    miRNA_TF = pd.read_csv("Other_output/miRNA_TF.csv", sep = "\t").loc[:, ["miRNA", "TF"]]
    miRNA_target_TF =  pd.read_csv("Other_output/miRNA_target(TF).csv", sep = "\t").loc[:, ["miRNA", "TF"]]
    merge_2 = pd.merge(miRNA_TF, miRNA_target_TF, how = "right").drop_duplicates(keep = "first")
    miRNA_mRNA = pd.concat([merge_2, miRNA_TF]).drop_duplicates(keep = False).loc[:, ["miRNA", "TF"]]
    miRNA_mRNA.columns = ["miRNA", "mRNA"]
    miRNA_mRNA.to_csv("Other_output/miRNA_mRNA.csv", sep = "\t")

    miRNA_TF = pd.read_csv("Other_output/miRNA_TF.csv", sep = "\t").loc[:, ["miRNA", "TF"]]
    miRNA_mRNA = pd.read_csv("Other_output/miRNA_mRNA.csv", sep = "\t").loc[:, ["miRNA", "mRNA"]]
    miRNA_TF_mRNA = pd.merge(miRNA_TF, miRNA_mRNA, how = "inner").drop_duplicates(keep = "first")

    with open("Output/miRNA_TF_mRNA_1.txt", "w") as f2:
        f2.write("Open Circuits miRNA_TF_mRNA 1\n\n\nTF\t\tmiRNA\t\tmRNA\n\n\n")
        for i, j in miRNA_TF_mRNA.iterrows():
            f2.write(j["TF"] + "\t"  + "|------" + "\t" + j["miRNA"] + "\t" + "------|" + "\t" + j["mRNA"] + "\n\n")
        f2.close()

    merge_edited = pd.read_csv("Other_output/merge_edited.csv", sep = "\t").loc[:, ["mRNA", "miRNA"]]
    merge_edited.columns = ["TF", "miRNA"]
    TransmiR = pd.read_csv("Edited_data/TransmiR_new_2.csv", sep = "\t")
    TF_miRNA = (pd.merge(merge_edited, TransmiR, how = "inner").drop_duplicates(keep = "first")).to_csv("Other_output/TF_miRNA.csv", sep = "\t")

    TF_miRNA = pd.read_csv("Other_output/TF_miRNA.csv", sep = "\t").loc[:, ["TF", "miRNA"]]
    miRNA_TF_mRNA_2 = pd.merge(TF_miRNA, miRNA_TF_mRNA, how = "inner").drop_duplicates(keep = "first")

    with open("Output/miRNA_TF_mRNA_2.txt", "w") as f3:
        f3.write("Open Circuits miRNA_TF_mRNA 2\n\n\nTF\t\tmiRNA\t\tmRNA\n\n\n")
        for i, j in miRNA_TF_mRNA_2.iterrows():
            f3.write("         ------>\n" + j["TF"] + "\t"  + "|------" + "\t" + j["miRNA"] + "\t" + "------|" + "\t" + j["mRNA"] + "\n\n\n")
        f3.close()

    hTFtarget = pd.read_csv("Edited_data/hTFtarget_2.csv", sep = "\t").loc[:, ["mRNA", "TF"]]     
    miRNA_TF_mRNA_3 = pd.merge(hTFtarget, miRNA_TF_mRNA_2).drop_duplicates(keep = "first")       
    miRNA_TF_mRNA_3.to_csv("Other_output/mRNA_TF_miRNA.csv", sep = "\t")

    with open("Output/Circuits.txt", "w") as f4:
        f4.write("CLOSED CIRCUITS\n\n\nmRNA\t\tTF\t\tmiRNA\n\n\n")
        miRNA_TF_mRNA_3 = pd.read_csv("Other_output/mRNA_TF_miRNA.csv", sep = "\t")
        for i, j in miRNA_TF_mRNA_3.iterrows():
            f4.write("\n\n" + "\t\t  |----------\n" + j["mRNA"] + "\t"  + "<---" + "\t" + j["TF"] + "\t" + "--->" + "\t" + j["miRNA"] + "\n"  + "\t" + "|-------------------" + "\n\n\n")
        f4.close()   


    
    with open("Output/Circuits_list.txt", "w") as f5:
        f5.write("CLOSED CIRCUITS\n\nmRNA\tTF\tmiRNA\n")
        miRNA_TF_mRNA_3 = pd.read_csv("Other_output/mRNA_TF_miRNA.csv", sep = "\t")
        for i, j in miRNA_TF_mRNA_3.iterrows():
            f5.write(j["mRNA"] + "\t" + j["TF"] + "\t" + j["miRNA"] + "\n")
        f5.close()  


def difexp2_yedek(file_name_1, file_name_2, output_1, output_2):
    import pandas as pd 
    import math as mt
    folder_1 = pd.read_csv(file_name_1, sep = "\t")
    folder_2 = pd.read_csv(file_name_2, sep = "\t")
    for i, j in folder_1.iterrows():
        folder_1.loc[i, "log(FC)"] = mt.log((j["Tumor"]/j["Normal"]), 2)
    for k, l in folder_2.iterrows():
        folder_2.loc[k, "log(FC)"] = mt.log((l["Tumor"]/l["Normal"]), 2)

    thrs = 1.5                                                             
    folder_DE_1 = folder_1[abs(folder_1["log(FC)"]) > thrs] 

    thrs = 2 
    folder_DE_2 = folder_2[abs(folder_2["log(FC)"]) > thrs]

    folder_DE_1_sorted = folder_DE_1.sort_values("log(FC)", axis = 0, ascending = False)
    for a, b in folder_DE_1_sorted.iterrows():
        if b["log(FC)"] >= 0:
            folder_DE_1_sorted.loc[a, "Expression"] = "Up"
        else: 
            folder_DE_1_sorted.loc[a, "Expression"] = "Down"

    folder_DE_2_sorted = folder_DE_2.sort_values("log(FC)", axis = 0, ascending = False)
    for c, d in folder_DE_2_sorted.iterrows():
        if d["log(FC)"] >= 0:
            folder_DE_2_sorted.loc[c, "Expression"] = "Up"
        else: 
            folder_DE_2_sorted.loc[c, "Expression"] = "Down"

    folder_DE_1_sorted.to_csv(output_1, sep = "\t")
    folder_DE_2_sorted.to_csv(output_2, sep = "\t")


    with open("Other_output/merge.txt", "w") as f1:
        f1.write("mRNA\tExpression\tExpression\tmiRNA\n")
        for m, n in folder_DE_1_sorted.iterrows():
            for o, p in folder_DE_2_sorted.iterrows():
                if (n["Expression"] == "Up") & (p["Expression"] == "Down") | (n["Expression"] == "Down") & (p["Expression"] == "Up"):
                    f1.write("{0}\t{1}\t{2}\t{3}\n".format(n["mRNA"], n["Expression"], p["Expression"], p["miRNA"]))
        f1.close()
    
    merge = ((pd.read_csv("Other_output/merge.txt", sep = "\t").drop_duplicates(keep = "first")).loc[:, ["mRNA", "miRNA"]]).to_csv("Other_output/merge_edited.csv", sep = "\t")
    
    TargetScan = pd.read_csv("Edited_data/TargetScan.csv", sep = "\t")
    merge = pd.read_csv("Other_output/merge_edited.csv", sep = "\t")
    miRNA_target = (pd.merge(TargetScan, merge, how = "inner").drop_duplicates(keep = "first").loc[:, ["miRNA", "mRNA"]]).to_csv("Other_output/miRNA_target(mRNA).csv", sep = "\t")

    miRNA_target = pd.read_csv("Other_output/miRNA_target(mRNA).csv", sep = "\t").loc[:, ["miRNA", "mRNA"]]
    miRNA_target.columns = ["miRNA", "TF"]
    miRNA_target.to_csv("Other_output/miRNA_target(TF).csv", sep = "\t")

    hTFtarget = pd.read_csv("Edited_data/hTFtarget.csv", sep = "\t").loc[:, ["TF"]]
    miRNA_target_TF = pd.read_csv("Other_output/miRNA_target(TF).csv", sep = "\t").loc[:, ["miRNA", "TF"]]
    miRNA_TF = ((pd.merge(hTFtarget, miRNA_target_TF, how = "inner").loc[:, ["miRNA", "TF"]]).drop_duplicates(keep = "first")).to_csv("Other_output/miRNA_TF.csv", sep = "\t")

    miRNA_TF = pd.read_csv("Other_output/miRNA_TF.csv", sep = "\t").loc[:, ["miRNA", "TF"]]
    miRNA_target_TF =  pd.read_csv("Other_output/miRNA_target(TF).csv", sep = "\t").loc[:, ["miRNA", "TF"]]
    merge_2 = pd.merge(miRNA_TF, miRNA_target_TF, how = "right").drop_duplicates(keep = "first")
    miRNA_mRNA = pd.concat([merge_2, miRNA_TF]).drop_duplicates(keep = False).loc[:, ["miRNA", "TF"]]
    miRNA_mRNA.columns = ["miRNA", "mRNA"]
    miRNA_mRNA.to_csv("Other_output/miRNA_mRNA.csv", sep = "\t")

    miRNA_TF = pd.read_csv("Other_output/miRNA_TF.csv", sep = "\t").loc[:, ["miRNA", "TF"]]
    miRNA_mRNA = pd.read_csv("Other_output/miRNA_mRNA.csv", sep = "\t").loc[:, ["miRNA", "mRNA"]]
    miRNA_TF_mRNA = pd.merge(miRNA_TF, miRNA_mRNA, how = "inner").drop_duplicates(keep = "first")

    with open("Output/miRNA_TF_mRNA_1.txt", "w") as f2:
        f2.write("Open Circuits miRNA_TF_mRNA 1\n\n\nTF\t\tmiRNA\t\tmRNA\n\n\n")
        for i, j in miRNA_TF_mRNA.iterrows():
            f2.write(j["TF"] + "\t"  + "|------" + "\t" + j["miRNA"] + "\t" + "------|" + "\t" + j["mRNA"] + "\n\n")
        f2.close()

    merge_edited = pd.read_csv("Other_output/merge_edited.csv", sep = "\t").loc[:, ["mRNA", "miRNA"]]
    merge_edited.columns = ["TF", "miRNA"]
    TransmiR = pd.read_csv("Edited_data/TransmiR_new_2.csv", sep = "\t")
    TF_miRNA = (pd.merge(merge_edited, TransmiR, how = "inner").drop_duplicates(keep = "first")).to_csv("Other_output/TF_miRNA.csv", sep = "\t")

    TF_miRNA = pd.read_csv("Other_output/TF_miRNA.csv", sep = "\t").loc[:, ["TF", "miRNA"]]
    miRNA_TF_mRNA_2 = pd.merge(TF_miRNA, miRNA_TF_mRNA, how = "inner").drop_duplicates(keep = "first")

    with open("Output/miRNA_TF_mRNA_2.txt", "w") as f3:
        f3.write("Open Circuits miRNA_TF_mRNA 2\n\n\nTF\t\tmiRNA\t\tmRNA\n\n\n")
        for i, j in miRNA_TF_mRNA_2.iterrows():
            f3.write("         ------>\n" + j["TF"] + "\t"  + "|------" + "\t" + j["miRNA"] + "\t" + "------|" + "\t" + j["mRNA"] + "\n\n\n")
        f3.close()

    hTFtarget = pd.read_csv("Edited_data/hTFtarget_2.csv", sep = "\t").loc[:, ["mRNA", "TF"]]     
    miRNA_TF_mRNA_3 = pd.merge(hTFtarget, miRNA_TF_mRNA_2).drop_duplicates(keep = "first")       
    miRNA_TF_mRNA_3.to_csv("Other_output/mRNA_TF_miRNA.csv", sep = "\t")

    with open("Output/Circuits.txt", "w") as f4:
        f4.write("CLOSED CIRCUITS\n\n\nmRNA\t\tTF\t\tmiRNA\n\n\n")
        miRNA_TF_mRNA_3 = pd.read_csv("Other_output/mRNA_TF_miRNA.csv", sep = "\t")
        for i, j in miRNA_TF_mRNA_3.iterrows():
            f4.write("\n\n" + "\t\t  |----------\n" + j["mRNA"] + "\t"  + "<---" + "\t" + j["TF"] + "\t" + "--->" + "\t" + j["miRNA"] + "\n"  + "\t" + "|-------------------" + "\n\n\n")
        f4.close()   


    
    with open("Output/Circuits_list.txt", "w") as f5:
        f5.write("CLOSED CIRCUITS\n\nmRNA\tTF\tmiRNA\n")
        miRNA_TF_mRNA_3 = pd.read_csv("Other_output/mRNA_TF_miRNA.csv", sep = "\t")
        for i, j in miRNA_TF_mRNA_3.iterrows():
            f5.write(j["mRNA"] + "\t" + j["TF"] + "\t" + j["miRNA"] + "\n")
        f5.close()  

def difexp2_yedek(file_name_1, file_name_2, output_1, output_2):
    import pandas as pd 
    import math as mt
    folder_1 = pd.read_csv(file_name_1, sep = "\t")
    folder_2 = pd.read_csv(file_name_2, sep = "\t")
    for i, j in folder_1.iterrows():
        folder_1.loc[i, "log(FC)"] = mt.log((j["Tumor"]/j["Normal"]), 2)
    for k, l in folder_2.iterrows():
        folder_2.loc[k, "log(FC)"] = mt.log((l["Tumor"]/l["Normal"]), 2)

    thrs = 1.5                                                             
    folder_DE_1 = folder_1[abs(folder_1["log(FC)"]) > thrs] 

    thrs = 2 
    folder_DE_2 = folder_2[abs(folder_2["log(FC)"]) > thrs]

    folder_DE_1_sorted = folder_DE_1.sort_values("log(FC)", axis = 0, ascending = False)
    for a, b in folder_DE_1_sorted.iterrows():
        if b["log(FC)"] >= 0:
            folder_DE_1_sorted.loc[a, "Expression"] = "Up"
        else: 
            folder_DE_1_sorted.loc[a, "Expression"] = "Down"

    folder_DE_2_sorted = folder_DE_2.sort_values("log(FC)", axis = 0, ascending = False)
    for c, d in folder_DE_2_sorted.iterrows():
        if d["log(FC)"] >= 0:
            folder_DE_2_sorted.loc[c, "Expression"] = "Up"
        else: 
            folder_DE_2_sorted.loc[c, "Expression"] = "Down"

    folder_DE_1_sorted.to_csv(output_1, sep = "\t")
    folder_DE_2_sorted.to_csv(output_2, sep = "\t")

    f = open("Other_output/merge.txt", "w")
    f.write("mRNA\tExpression\tExpression\tmiRNA\n")
    for m, n in folder_DE_1_sorted.iterrows():
        for o, p in folder_DE_2_sorted.iterrows():
            if (n["Expression"] == "Up") & (p["Expression"] == "Down") | (n["Expression"] == "Down") & (p["Expression"] == "Up"):
                f.write("{0}\t{1}\t{2}\t{3}\n".format(n["mRNA"], n["Expression"], p["Expression"], p["miRNA"]))
    f.close()
    
    merge = ((pd.read_csv("Other_output/merge.txt", sep = "\t").drop_duplicates(keep = "first")).loc[:, ["mRNA", "miRNA"]]).to_csv("Other_output/merge_edited.csv", sep = "\t")
    
    TargetScan = pd.read_csv("Edited_data/TargetScan.csv", sep = "\t")
    merge = pd.read_csv("Other_output/merge_edited.csv", sep = "\t")
    miRNA_target = (pd.merge(TargetScan, merge, how = "inner").drop_duplicates(keep = "first").loc[:, ["miRNA", "mRNA"]]).to_csv("Other_output/miRNA_target(mRNA).csv", sep = "\t")

    miRNA_target = pd.read_csv("Other_output/miRNA_target(mRNA).csv", sep = "\t").loc[:, ["miRNA", "mRNA"]]
    miRNA_target.columns = ["miRNA", "TF"]
    miRNA_target.to_csv("Other_output/miRNA_target(TF).csv", sep = "\t")

    hTFtarget = pd.read_csv("Edited_data/hTFtarget.csv", sep = "\t").loc[:, ["TF"]]
    miRNA_target_TF = pd.read_csv("Other_output/miRNA_target(TF).csv", sep = "\t").loc[:, ["miRNA", "TF"]]
    miRNA_TF = ((pd.merge(hTFtarget, miRNA_target_TF, how = "inner").loc[:, ["miRNA", "TF"]]).drop_duplicates(keep = "first")).to_csv("Other_output/miRNA_TF.csv", sep = "\t")

    miRNA_TF = pd.read_csv("Other_output/miRNA_TF.csv", sep = "\t").loc[:, ["miRNA", "TF"]]
    miRNA_target_TF =  pd.read_csv("Other_output/miRNA_target(TF).csv", sep = "\t").loc[:, ["miRNA", "TF"]]
    merge_2 = pd.merge(miRNA_TF, miRNA_target_TF, how = "right").drop_duplicates(keep = "first")
    miRNA_mRNA = pd.concat([merge_2, miRNA_TF]).drop_duplicates(keep = False).loc[:, ["miRNA", "TF"]]
    miRNA_mRNA.columns = ["miRNA", "mRNA"]
    miRNA_mRNA.to_csv("Other_output/miRNA_mRNA.csv", sep = "\t")

    miRNA_TF = pd.read_csv("Other_output/miRNA_TF.csv", sep = "\t").loc[:, ["miRNA", "TF"]]
    miRNA_mRNA = pd.read_csv("Other_output/miRNA_mRNA.csv", sep = "\t").loc[:, ["miRNA", "mRNA"]]
    miRNA_TF_mRNA = pd.merge(miRNA_TF, miRNA_mRNA, how = "inner").drop_duplicates(keep = "first")

    f = open("Output/miRNA_TF_mRNA_1.txt", "w")
    f.write("Open Circuits miRNA_TF_mRNA 1\n\n\nTF\t\tmiRNA\t\tmRNA\n\n\n")
    for i, j in miRNA_TF_mRNA.iterrows():
        f.write(j["TF"] + "\t"  + "|------" + "\t" + j["miRNA"] + "\t" + "------|" + "\t" + j["mRNA"] + "\n\n")
    f.close()

    merge_edited = pd.read_csv("Other_output/merge_edited.csv", sep = "\t").loc[:, ["mRNA", "miRNA"]]
    merge_edited.columns = ["TF", "miRNA"]
    TransmiR = pd.read_csv("Edited_data/TransmiR_new_2.csv", sep = "\t")
    TF_miRNA = (pd.merge(merge_edited, TransmiR, how = "inner").drop_duplicates(keep = "first")).to_csv("Other_output/TF_miRNA.csv", sep = "\t")

    TF_miRNA = pd.read_csv("Other_output/TF_miRNA.csv", sep = "\t").loc[:, ["TF", "miRNA"]]
    miRNA_TF_mRNA_2 = pd.merge(TF_miRNA, miRNA_TF_mRNA, how = "inner").drop_duplicates(keep = "first")

    f = open("Output/miRNA_TF_mRNA_2.txt", "w")
    f.write("Open Circuits miRNA_TF_mRNA 2\n\n\nTF\t\tmiRNA\t\tmRNA\n\n\n")
    for i, j in miRNA_TF_mRNA_2.iterrows():
        f.write("         ------>\n" + j["TF"] + "\t"  + "|------" + "\t" + j["miRNA"] + "\t" + "------|" + "\t" + j["mRNA"] + "\n\n\n")
    f.close()

    hTFtarget = pd.read_csv("Edited_data/hTFtarget_2.csv", sep = "\t").loc[:, ["mRNA", "TF"]]     
    miRNA_TF_mRNA_3 = pd.merge(hTFtarget, miRNA_TF_mRNA_2).drop_duplicates(keep = "first")       
    miRNA_TF_mRNA_3.to_csv("Other_output/mRNA_TF_miRNA.csv", sep = "\t")

    f = open("Output/Circuits.txt", "w")
    f.write("CLOSED CIRCUITS\n\n\nmRNA\t\tTF\t\tmiRNA\n\n\n")
    miRNA_TF_mRNA_3 = pd.read_csv("Other_output/mRNA_TF_miRNA.csv", sep = "\t")
    for i, j in miRNA_TF_mRNA_3.iterrows():
        f.write("\n\n" + "\t\t  |----------\n" + j["mRNA"] + "\t"  + "<---" + "\t" + j["TF"] + "\t" + "--->" + "\t" + j["miRNA"] + "\n"  + "\t" + "|-------------------" + "\n\n\n")
    f.close()   


    f = open("Output/Circuits_list.txt", "w")
    f.write("CLOSED CIRCUITS\n\nmRNA\tTF\tmiRNA\n")
    miRNA_TF_mRNA_3 = pd.read_csv("Other_output/mRNA_TF_miRNA.csv", sep = "\t")
    for i, j in miRNA_TF_mRNA_3.iterrows():
        f.write(j["mRNA"] + "\t" + j["TF"] + "\t" + j["miRNA"] + "\n")
    f.close()  


def difexp3(file_name_1, file_name_2, output_1, output_2):
    import pandas as pd 
    import math as mt
    folder_1 = pd.read_csv(file_name_1, sep = "\t")
    folder_2 = pd.read_csv(file_name_2, sep = "\t")
    for i, j in folder_1.iterrows():
        folder_1.loc[i, "log(FC)"] = mt.log((j["Tumor"]/j["Normal"]), 2)
    for k, l in folder_2.iterrows():
        folder_2.loc[k, "log(FC)"] = mt.log((l["Tumor"]/l["Normal"]), 2)

    thrs = 2                                                             # Gerçek kod için değerler değiştirilecek !
    folder_DE_1 = folder_1[abs(folder_1["log(FC)"]) > thrs] 

    thrs = 1.5 
    folder_DE_2 = folder_2[abs(folder_2["log(FC)"]) > thrs]

    folder_DE_1_sorted = folder_DE_1.sort_values("log(FC)", axis = 0, ascending = False)
    for a, b in folder_DE_1_sorted.iterrows():
        if b["log(FC)"] >= 0:
            folder_DE_1_sorted.loc[a, "Expression"] = "Up"
        else: 
            folder_DE_1_sorted.loc[a, "Expression"] = "Down"

    folder_DE_2_sorted = folder_DE_2.sort_values("log(FC)", axis = 0, ascending = False)
    for c, d in folder_DE_2_sorted.iterrows():
        if d["log(FC)"] >= 0:
            folder_DE_2_sorted.loc[c, "Expression"] = "Up"
        else: 
            folder_DE_2_sorted.loc[c, "Expression"] = "Down"

    folder_DE_1_sorted.to_csv(output_1, sep = "\t")
    folder_DE_2_sorted.to_csv(output_2, sep = "\t")

    f = open("Other_output/merge.txt", "w")
    f.write("mRNA\tExpression\tExpression\tmiRNA\n")
    for m, n in folder_DE_1_sorted.iterrows():
        for o, p in folder_DE_2_sorted.iterrows():
            if (n["Expression"] == "Up") & (p["Expression"] == "Down") | (n["Expression"] == "Down") & (p["Expression"] == "Up"):
                f.write("{0}\t{1}\t{2}\t{3}\n".format(n["mRNA"], n["Expression"], p["Expression"], p["miRNA"]))
    f.close()
    
    merge = ((pd.read_csv("Other_output/merge.txt", sep = "\t").drop_duplicates(keep = "first")).loc[:, ["mRNA", "miRNA"]]).to_csv("Other_output/merge_edited.csv", sep = "\t")
    
    TargetScan = pd.read_csv("Edited_data/TargetScan.csv", sep = "\t")
    merge = pd.read_csv("Other_output/merge_edited.csv", sep = "\t")
    miRNA_target = (pd.merge(TargetScan, merge, how = "inner").drop_duplicates(keep = "first").loc[:, ["miRNA", "mRNA"]]).to_csv("Other_output/miRNA_target(mRNA).csv", sep = "\t")

    miRNA_target = pd.read_csv("Other_output/miRNA_target(mRNA).csv", sep = "\t").loc[:, ["miRNA", "mRNA"]]
    miRNA_target.columns = ["miRNA", "TF"]
    miRNA_target.to_csv("Other_output/miRNA_target(TF).csv", sep = "\t")

    hTFtarget = pd.read_csv("Edited_data/hTFtarget.csv", sep = "\t").loc[:, ["TF"]]
    miRNA_target_TF = pd.read_csv("Other_output/miRNA_target(TF).csv", sep = "\t").loc[:, ["miRNA", "TF"]]
    miRNA_TF = ((pd.merge(hTFtarget, miRNA_target_TF, how = "inner").loc[:, ["miRNA", "TF"]]).drop_duplicates(keep = "first")).to_csv("Other_output/miRNA_TF.csv", sep = "\t")

    miRNA_TF = pd.read_csv("Other_output/miRNA_TF.csv", sep = "\t").loc[:, ["miRNA", "TF"]]
    miRNA_target_TF =  pd.read_csv("Other_output/miRNA_target(TF).csv", sep = "\t").loc[:, ["miRNA", "TF"]]
    merge_2 = pd.merge(miRNA_TF, miRNA_target_TF, how = "right").drop_duplicates(keep = "first")
    miRNA_mRNA = pd.concat([merge_2, miRNA_TF]).drop_duplicates(keep = False).loc[:, ["miRNA", "TF"]]
    miRNA_mRNA.columns = ["miRNA", "mRNA"]
    miRNA_mRNA.to_csv("Other_output/miRNA_mRNA.csv", sep = "\t")

    miRNA_TF = pd.read_csv("Other_output/miRNA_TF.csv", sep = "\t").loc[:, ["miRNA", "TF"]]
    miRNA_mRNA = pd.read_csv("Other_output/miRNA_mRNA.csv", sep = "\t").loc[:, ["miRNA", "mRNA"]]
    miRNA_TF_mRNA = pd.merge(miRNA_TF, miRNA_mRNA, how = "inner").drop_duplicates(keep = "first")

    f = open("Output/miRNA_TF_mRNA_1.txt", "w")
    f.write("Open Circuits miRNA_TF_mRNA 1\n\n\nTF\t\tmiRNA\t\tmRNA\n\n\n")
    for i, j in miRNA_TF_mRNA.iterrows():
        f.write(j["TF"] + "\t"  + "|------" + "\t" + j["miRNA"] + "\t" + "------|" + "\t" + j["mRNA"] + "\n\n")
    f.close()

    merge_edited = pd.read_csv("Other_output/merge_edited.csv", sep = "\t").loc[:, ["mRNA", "miRNA"]]
    merge_edited.columns = ["TF", "miRNA"]
    TransmiR = pd.read_csv("Edited_data/TransmiR_new_2.csv", sep = "\t")
    TF_miRNA = (pd.merge(merge_edited, TransmiR, how = "inner").drop_duplicates(keep = "first")).to_csv("Other_output/TF_miRNA.csv", sep = "\t")

    TF_miRNA = pd.read_csv("Other_output/TF_miRNA.csv", sep = "\t").loc[:, ["TF", "miRNA"]]
    miRNA_TF_mRNA_2 = pd.merge(TF_miRNA, miRNA_TF_mRNA, how = "inner").drop_duplicates(keep = "first")

    f = open("Output/miRNA_TF_mRNA_2.txt", "w")
    f.write("Open Circuits miRNA_TF_mRNA 2\n\n\nTF\t\tmiRNA\t\tmRNA\n\n\n")
    for i, j in miRNA_TF_mRNA_2.iterrows():
        f.write("         ------>\n" + j["TF"] + "\t"  + "|------" + "\t" + j["miRNA"] + "\t" + "------|" + "\t" + j["mRNA"] + "\n\n\n")
    f.close()

    hTFtarget = pd.read_csv("Edited_data/hTFtarget_2.csv", sep = "\t").loc[:, ["mRNA", "TF"]]    
    miRNA_TF_mRNA_3 = pd.merge(hTFtarget, miRNA_TF_mRNA_2).drop_duplicates(keep = "first")       
    miRNA_TF_mRNA_3.to_csv("Other_output/mRNA_TF_miRNA.csv", sep = "\t")

    f = open("Output/Circuits.txt", "w")
    f.write("CLOSED CIRCUITS\n\n\nmRNA\t\tTF\t\tmiRNA\n\n\n")
    miRNA_TF_mRNA_3 = pd.read_csv("Other_output/mRNA_TF_miRNA.csv", sep = "\t")
    for i, j in miRNA_TF_mRNA_3.iterrows():
        f.write("\n\n" + "\t\t  |----------\n" + j["mRNA"] + "\t"  + "<---" + "\t" + j["TF"] + "\t" + "--->" + "\t" + j["miRNA"] + "\n"  + "\t" + "|-------------------" + "\n\n\n")
    f.close()   


    
    f = open("Output/Circuits_list.txt", "w")
    f.write("CLOSED CIRCUITS\n\nmRNA\tTF\tmiRNA\n")
    miRNA_TF_mRNA_3 = pd.read_csv("Other_output/mRNA_TF_miRNA.csv", sep = "\t")
    for i, j in miRNA_TF_mRNA_3.iterrows():
        f.write(j["mRNA"] + "\t" + j["TF"] + "\t" + j["miRNA"] + "\n")
    f.close()  

def difexp4(file_name_1, file_name_2, output_1, output_2):
    import pandas as pd 
    import math as mt
    folder_1 = pd.read_csv(file_name_1, sep = "\t")
    folder_2 = pd.read_csv(file_name_2, sep = "\t")
    for i, j in folder_1.iterrows():
        folder_1.loc[i, "log(FC)"] = mt.log((j["Tumor"]/j["Normal"]), 2)
    for k, l in folder_2.iterrows():
        folder_2.loc[k, "log(FC)"] = mt.log((l["Tumor"]/l["Normal"]), 2)

    thrs = 2                                                             
    folder_DE_1 = folder_1[abs(folder_1["log(FC)"]) > thrs] 

    thrs = 2
    folder_DE_2 = folder_2[abs(folder_2["log(FC)"]) > thrs]

    folder_DE_1_sorted = folder_DE_1.sort_values("log(FC)", axis = 0, ascending = False)
    for a, b in folder_DE_1_sorted.iterrows():
        if b["log(FC)"] >= 0:
            folder_DE_1_sorted.loc[a, "Expression"] = "Up"
        else: 
            folder_DE_1_sorted.loc[a, "Expression"] = "Down"

    folder_DE_2_sorted = folder_DE_2.sort_values("log(FC)", axis = 0, ascending = False)
    for c, d in folder_DE_2_sorted.iterrows():
        if d["log(FC)"] >= 0:
            folder_DE_2_sorted.loc[c, "Expression"] = "Up"
        else: 
            folder_DE_2_sorted.loc[c, "Expression"] = "Down"

    folder_DE_1_sorted.to_csv(output_1, sep = "\t")
    folder_DE_2_sorted.to_csv(output_2, sep = "\t")

    f = open("Other_output/merge.txt", "w")
    f.write("mRNA\tExpression\tExpression\tmiRNA\n")
    for m, n in folder_DE_1_sorted.iterrows():
        for o, p in folder_DE_2_sorted.iterrows():
            if (n["Expression"] == "Up") & (p["Expression"] == "Down") | (n["Expression"] == "Down") & (p["Expression"] == "Up"):
                f.write("{0}\t{1}\t{2}\t{3}\n".format(n["mRNA"], n["Expression"], p["Expression"], p["miRNA"]))
    f.close()
    
    merge = ((pd.read_csv("Other_output/merge.txt", sep = "\t").drop_duplicates(keep = "first")).loc[:, ["mRNA", "miRNA"]]).to_csv("Other_output/merge_edited.csv", sep = "\t")
    
    TargetScan = pd.read_csv("Edited_data/TargetScan.csv", sep = "\t")
    merge = pd.read_csv("Other_output/merge_edited.csv", sep = "\t")
    miRNA_target = (pd.merge(TargetScan, merge, how = "inner").drop_duplicates(keep = "first").loc[:, ["miRNA", "mRNA"]]).to_csv("Other_output/miRNA_target(mRNA).csv", sep = "\t")

    miRNA_target = pd.read_csv("Other_output/miRNA_target(mRNA).csv", sep = "\t").loc[:, ["miRNA", "mRNA"]]
    miRNA_target.columns = ["miRNA", "TF"]
    miRNA_target.to_csv("Other_output/miRNA_target(TF).csv", sep = "\t")

    hTFtarget = pd.read_csv("Edited_data/hTFtarget.csv", sep = "\t").loc[:, ["TF"]]
    miRNA_target_TF = pd.read_csv("Other_output/miRNA_target(TF).csv", sep = "\t").loc[:, ["miRNA", "TF"]]
    miRNA_TF = ((pd.merge(hTFtarget, miRNA_target_TF, how = "inner").loc[:, ["miRNA", "TF"]]).drop_duplicates(keep = "first")).to_csv("Other_output/miRNA_TF.csv", sep = "\t")

    miRNA_TF = pd.read_csv("Other_output/miRNA_TF.csv", sep = "\t").loc[:, ["miRNA", "TF"]]
    miRNA_target_TF =  pd.read_csv("Other_output/miRNA_target(TF).csv", sep = "\t").loc[:, ["miRNA", "TF"]]
    merge_2 = pd.merge(miRNA_TF, miRNA_target_TF, how = "right").drop_duplicates(keep = "first")
    miRNA_mRNA = pd.concat([merge_2, miRNA_TF]).drop_duplicates(keep = False).loc[:, ["miRNA", "TF"]]
    miRNA_mRNA.columns = ["miRNA", "mRNA"]
    miRNA_mRNA.to_csv("Other_output/miRNA_mRNA.csv", sep = "\t")

    miRNA_TF = pd.read_csv("Other_output/miRNA_TF.csv", sep = "\t").loc[:, ["miRNA", "TF"]]
    miRNA_mRNA = pd.read_csv("Other_output/miRNA_mRNA.csv", sep = "\t").loc[:, ["miRNA", "mRNA"]]
    miRNA_TF_mRNA = pd.merge(miRNA_TF, miRNA_mRNA, how = "inner").drop_duplicates(keep = "first")

    f = open("Output/miRNA_TF_mRNA_1.txt", "w")
    f.write("Open Circuits miRNA_TF_mRNA 1\n\n\nTF\t\tmiRNA\t\tmRNA\n\n\n")
    for i, j in miRNA_TF_mRNA.iterrows():
        f.write(j["TF"] + "\t"  + "|------" + "\t" + j["miRNA"] + "\t" + "------|" + "\t" + j["mRNA"] + "\n\n")
    f.close()

    merge_edited = pd.read_csv("Other_output/merge_edited.csv", sep = "\t").loc[:, ["mRNA", "miRNA"]]
    merge_edited.columns = ["TF", "miRNA"]
    TransmiR = pd.read_csv("Edited_data/TransmiR_new_2.csv", sep = "\t")
    TF_miRNA = (pd.merge(merge_edited, TransmiR, how = "inner").drop_duplicates(keep = "first")).to_csv("Other_output/TF_miRNA.csv", sep = "\t")

    TF_miRNA = pd.read_csv("Other_output/TF_miRNA.csv", sep = "\t").loc[:, ["TF", "miRNA"]]
    miRNA_TF_mRNA_2 = pd.merge(TF_miRNA, miRNA_TF_mRNA, how = "inner").drop_duplicates(keep = "first")

    f = open("Output/miRNA_TF_mRNA_2.txt", "w")
    f.write("Open Circuits miRNA_TF_mRNA 2\n\n\nTF\t\tmiRNA\t\tmRNA\n\n\n")
    for i, j in miRNA_TF_mRNA_2.iterrows():
        f.write("         ------>\n" + j["TF"] + "\t"  + "|------" + "\t" + j["miRNA"] + "\t" + "------|" + "\t" + j["mRNA"] + "\n\n\n")
    f.close()

    hTFtarget = pd.read_csv("Edited_data/hTFtarget_2.csv", sep = "\t").loc[:, ["mRNA", "TF"]]     
    miRNA_TF_mRNA_3 = pd.merge(hTFtarget, miRNA_TF_mRNA_2).drop_duplicates(keep = "first")       
    miRNA_TF_mRNA_3.to_csv("Other_output/mRNA_TF_miRNA.csv", sep = "\t")

    f = open("Output/Circuits.txt", "w")
    f.write("CLOSED CIRCUITS\n\n\nmRNA\t\tTF\t\tmiRNA\n\n\n")
    miRNA_TF_mRNA_3 = pd.read_csv("Other_output/mRNA_TF_miRNA.csv", sep = "\t")
    for i, j in miRNA_TF_mRNA_3.iterrows():
        f.write("\n\n" + "\t\t  |----------\n" + j["mRNA"] + "\t"  + "<---" + "\t" + j["TF"] + "\t" + "--->" + "\t" + j["miRNA"] + "\n"  + "\t" + "|-------------------" + "\n\n\n")
    f.close()   


    
    f = open("Output/Circuits_list.txt", "w")
    f.write("CLOSED CIRCUITS\n\nmRNA\tTF\tmiRNA\n")
    miRNA_TF_mRNA_3 = pd.read_csv("Other_output/mRNA_TF_miRNA.csv", sep = "\t")
    for i, j in miRNA_TF_mRNA_3.iterrows():
        f.write(j["mRNA"] + "\t" + j["TF"] + "\t" + j["miRNA"] + "\n")
    f.close()  


def difexp5(file_name_1, file_name_2, output_1, output_2):
    import pandas as pd 
    import math as mt
    folder_DE_1 = pd.read_csv(file_name_1, sep = "\t")
    folder_DE_2 = pd.read_csv(file_name_2, sep = "\t")

    folder_DE_1_sorted = folder_DE_1.sort_values("log(FC)", axis = 0, ascending = False)
    for a, b in folder_DE_1_sorted.iterrows():
        if b["log(FC)"] >= 0:
            folder_DE_1_sorted.loc[a, "Expression"] = "Up"
        else: 
            folder_DE_1_sorted.loc[a, "Expression"] = "Down"

    folder_DE_2_sorted = folder_DE_2.sort_values("log(FC)", axis = 0, ascending = False)
    for c, d in folder_DE_2_sorted.iterrows():
        if d["log(FC)"] >= 0:
            folder_DE_2_sorted.loc[c, "Expression"] = "Up"
        else: 
            folder_DE_2_sorted.loc[c, "Expression"] = "Down"

    folder_DE_1_sorted.to_csv(output_1, sep = "\t")
    folder_DE_2_sorted.to_csv(output_2, sep = "\t")


    with open("Other_output/merge.txt", "w") as f1:
        f1.write("mRNA\tExpression\tExpression\tmiRNA\n")
        for m, n in folder_DE_1_sorted.iterrows():
            for o, p in folder_DE_2_sorted.iterrows():
                if (n["Expression"] == "Up") & (p["Expression"] == "Down") | (n["Expression"] == "Down") & (p["Expression"] == "Up"):
                    f1.write("{0}\t{1}\t{2}\t{3}\n".format(n["mRNA"], n["Expression"], p["Expression"], p["miRNA"]))
        f1.close()
    
    merge = ((pd.read_csv("Other_output/merge.txt", sep = "\t").drop_duplicates(keep = "first")).loc[:, ["mRNA", "miRNA"]]).to_csv("Other_output/merge_edited.csv", sep = "\t")
    
    TargetScan = pd.read_csv("Edited_data/TargetScan.csv", sep = "\t")
    merge = pd.read_csv("Other_output/merge_edited.csv", sep = "\t")
    miRNA_target = (pd.merge(TargetScan, merge, how = "inner").drop_duplicates(keep = "first").loc[:, ["miRNA", "mRNA"]]).to_csv("Other_output/miRNA_target(mRNA).csv", sep = "\t")

    miRNA_target = pd.read_csv("Other_output/miRNA_target(mRNA).csv", sep = "\t").loc[:, ["miRNA", "mRNA"]]
    miRNA_target.columns = ["miRNA", "TF"]
    miRNA_target.to_csv("Other_output/miRNA_target(TF).csv", sep = "\t")

    hTFtarget = pd.read_csv("Edited_data/hTFtarget.csv", sep = "\t").loc[:, ["TF"]]
    miRNA_target_TF = pd.read_csv("Other_output/miRNA_target(TF).csv", sep = "\t").loc[:, ["miRNA", "TF"]]
    miRNA_TF = ((pd.merge(hTFtarget, miRNA_target_TF, how = "inner").loc[:, ["miRNA", "TF"]]).drop_duplicates(keep = "first")).to_csv("Other_output/miRNA_TF.csv", sep = "\t")

    miRNA_TF = pd.read_csv("Other_output/miRNA_TF.csv", sep = "\t").loc[:, ["miRNA", "TF"]]
    miRNA_target_TF =  pd.read_csv("Other_output/miRNA_target(TF).csv", sep = "\t").loc[:, ["miRNA", "TF"]]
    merge_2 = pd.merge(miRNA_TF, miRNA_target_TF, how = "right").drop_duplicates(keep = "first")
    miRNA_mRNA = pd.concat([merge_2, miRNA_TF]).drop_duplicates(keep = False).loc[:, ["miRNA", "TF"]]
    miRNA_mRNA.columns = ["miRNA", "mRNA"]
    miRNA_mRNA.to_csv("Other_output/miRNA_mRNA.csv", sep = "\t")

    miRNA_TF = pd.read_csv("Other_output/miRNA_TF.csv", sep = "\t").loc[:, ["miRNA", "TF"]]
    miRNA_mRNA = pd.read_csv("Other_output/miRNA_mRNA.csv", sep = "\t").loc[:, ["miRNA", "mRNA"]]
    miRNA_TF_mRNA = pd.merge(miRNA_TF, miRNA_mRNA, how = "inner").drop_duplicates(keep = "first")

    with open("Output/miRNA_TF_mRNA_1.txt", "w") as f2:
        f2.write("Open Circuits miRNA_TF_mRNA 1\n\n\nTF\t\tmiRNA\t\tmRNA\n\n\n")
        for i, j in miRNA_TF_mRNA.iterrows():
            f2.write(j["TF"] + "\t"  + "|------" + "\t" + j["miRNA"] + "\t" + "------|" + "\t" + j["mRNA"] + "\n\n")
        f2.close()

    merge_edited = pd.read_csv("Other_output/merge_edited.csv", sep = "\t").loc[:, ["mRNA", "miRNA"]]
    merge_edited.columns = ["TF", "miRNA"]
    TransmiR = pd.read_csv("Edited_data/TransmiR_new_2.csv", sep = "\t")
    TF_miRNA = (pd.merge(merge_edited, TransmiR, how = "inner").drop_duplicates(keep = "first")).to_csv("Other_output/TF_miRNA.csv", sep = "\t")

    TF_miRNA = pd.read_csv("Other_output/TF_miRNA.csv", sep = "\t").loc[:, ["TF", "miRNA"]]
    miRNA_TF_mRNA_2 = pd.merge(TF_miRNA, miRNA_TF_mRNA, how = "inner").drop_duplicates(keep = "first")

    with open("Output/miRNA_TF_mRNA_2.txt", "w") as f3:
        f3.write("Open Circuits miRNA_TF_mRNA 2\n\n\nTF\t\tmiRNA\t\tmRNA\n\n\n")
        for i, j in miRNA_TF_mRNA_2.iterrows():
            f3.write("         ------>\n" + j["TF"] + "\t"  + "|------" + "\t" + j["miRNA"] + "\t" + "------|" + "\t" + j["mRNA"] + "\n\n\n")
        f3.close()

    hTFtarget = pd.read_csv("Edited_data/hTFtarget_2.csv", sep = "\t").loc[:, ["mRNA", "TF"]]     
    miRNA_TF_mRNA_3 = pd.merge(hTFtarget, miRNA_TF_mRNA_2).drop_duplicates(keep = "first")       
    miRNA_TF_mRNA_3.to_csv("Other_output/mRNA_TF_miRNA.csv", sep = "\t")

    with open("Output/Circuits.txt", "w") as f4:
        f4.write("CLOSED CIRCUITS\n\n\nmRNA\t\tTF\t\tmiRNA\n\n\n")
        miRNA_TF_mRNA_3 = pd.read_csv("Other_output/mRNA_TF_miRNA.csv", sep = "\t")
        for i, j in miRNA_TF_mRNA_3.iterrows():
            f4.write("\n\n" + "\t\t  |----------\n" + j["mRNA"] + "\t"  + "<---" + "\t" + j["TF"] + "\t" + "--->" + "\t" + j["miRNA"] + "\n"  + "\t" + "|-------------------" + "\n\n\n")
        f4.close()   


    
    with open("Output/Circuits_list.txt", "w") as f5:
        f5.write("CLOSED CIRCUITS\n\nmRNA\tTF\tmiRNA\n")
        miRNA_TF_mRNA_3 = pd.read_csv("Other_output/mRNA_TF_miRNA.csv", sep = "\t")
        for i, j in miRNA_TF_mRNA_3.iterrows():
            f5.write(j["mRNA"] + "\t" + j["TF"] + "\t" + j["miRNA"] + "\n")
        f5.close()  



@app.route('/', methods=['POST'])   #file post edildiginde burasi cagirilsin diye
def upload_file():
    
    mrna_expression_file = request.files['file1']
    miRNA_expression_file = request.files['file2']
    

    threshold_mrna = request.form.get('options-1', type=float)
    threshold_mirna = request.form.get('options-2', type=float)
 

    mrna_expression_file.save(mrna_expression_file.filename)
    miRNA_expression_file.save(miRNA_expression_file.filename)

    
    if request.form.get('ForFileType') == 'Raw Data' and request.form.get('analysis1') == 'without GSEA and miSEA'  :
        if(threshold_mrna == 1.5 and threshold_mirna == 1.5):
          difexp1(mrna_expression_file.filename, miRNA_expression_file.filename, "Output/DE_mRNA.csv", "Output/DE_miRNA.csv")
        
        elif(threshold_mrna == 1.5 and threshold_mirna == 2):
          difexp2(mrna_expression_file.filename, miRNA_expression_file.filename, "Output/DE_mRNA.csv", "Output/DE_miRNA.csv")
        
        elif(threshold_mrna == 2 and threshold_mirna == 1.5):
          difexp3(mrna_expression_file.filename, miRNA_expression_file.filename, "Output/DE_mRNA.csv", "Output/DE_miRNA.csv")
        
        elif(threshold_mrna == 2 and threshold_mirna == 2 ):
          difexp4(mrna_expression_file.filename, miRNA_expression_file.filename, "Output/DE_mRNA.csv", "Output/DE_miRNA.csv")
    
    elif  request.form.get('ForFileType') == 'Expression Values'  and request.form.get('analysis1') == 'without GSEA and miSEA' :
        difexp5(mrna_expression_file.filename, miRNA_expression_file.filename, "Output/DE_mRNA.csv", "Output/DE_miRNA.csv")
    
                
    app.config['SEND_FILE_MAX_AGE_DEFAULT'] = 0   #to change file content in every file changing       
    return redirect(url_for('index'))  #kendi uzerine donsun diye index yazilir

@app.route('/upload1') #to see content of file without downloading
def upload1():
    import csv
    with open('result.txt', "w") as my_output_file:
        with open('Raw_data_mRNA_Files.csv', "r") as my_input_file:
             [ my_output_file.write(" ".join(row)+'\n') for row in csv.reader(my_input_file)]
        my_output_file.close()           
    return send_file("result.txt")

@app.route('/upload2') #to see content of file without downloading
def upload2():
    import csv
    with open('result.txt', "w") as my_output_file:
        with open('Raw_data_miRNA_Files.csv', "r") as my_input_file:
             [ my_output_file.write(" ".join(row)+'\n') for row in csv.reader(my_input_file)]
        my_output_file.close()           
    return send_file("result.txt")

@app.route('/upload3') #to see content of file without downloading
def upload3():
    import csv
    with open('result.txt', "w") as my_output_file:
        with open('Expression_Values_mRNA_File.csv', "r") as my_input_file:
             [ my_output_file.write(" ".join(row)+'\n') for row in csv.reader(my_input_file)]
        my_output_file.close()           
    return send_file("result.txt")

@app.route('/upload4') #to see content of file without downloading
def upload4():
    import csv
    with open('result.txt', "w") as my_output_file:
        with open('Expression_Values_miRNA_File.csv', "r") as my_input_file:
             [ my_output_file.write(" ".join(row)+'\n') for row in csv.reader(my_input_file)]
        my_output_file.close()           
    return send_file("result.txt")

@app.route('/upload5') #to see content of file without downloading
def upload5():
    import csv
    with open('result.txt', "w") as my_output_file:
        with open('Output/DE_mRNA.csv', "r") as my_input_file:
             [ my_output_file.write(" ".join(row)+'\n') for row in csv.reader(my_input_file)]
        my_output_file.close()           
    return send_file("result.txt")

@app.route('/upload6') #to see content of file without downloading
def upload6():
    import csv
    with open('result.txt', "w") as my_output_file:
        with open('Output/DE_miRNA.csv', "r") as my_input_file:
             [ my_output_file.write(" ".join(row)+'\n') for row in csv.reader(my_input_file)]
        my_output_file.close()           
    return send_file("result.txt")

@app.route('/upload7') #to see content of file without downloading
def upload7():
    import csv
    with open('result.txt', "w") as my_output_file:
        with open('Output/miRNA_TF_mRNA_1.txt', "r") as my_input_file:
             [ my_output_file.write(" ".join(row)+'\n') for row in csv.reader(my_input_file)]
        my_output_file.close()           
    return send_file("result.txt")


@app.route('/upload8') #to see content of file without downloading
def upload8():
    import csv
    with open('result.txt', "w") as my_output_file:
        with open('Output/miRNA_TF_mRNA_2.txt', "r") as my_input_file:
             [ my_output_file.write(" ".join(row)+'\n') for row in csv.reader(my_input_file)]
        my_output_file.close()           
    return send_file("result.txt")


@app.route('/upload9') #to see content of file without downloading
def upload9():
    import csv
    with open('result.txt', "w") as my_output_file:
        with open('Output/Circuits.txt', "r") as my_input_file:
             [ my_output_file.write(" ".join(row)+'\n') for row in csv.reader(my_input_file)]
        my_output_file.close()           
    return send_file("result.txt")

@app.route('/upload10') #to see content of file without downloading
def upload10():
    import csv
    with open('result.txt', "w") as my_output_file:
        with open('Output/Circuits_list.txt', "r") as my_input_file:
             [ my_output_file.write(" ".join(row)+'\n') for row in csv.reader(my_input_file)]
        my_output_file.close()           
    return send_file("result.txt")







@app.route('/download')
def download_file():        #to download TF file
    try:
        path = "Output/DE_mRNA.csv"
        return send_file(path , as_attachment=True)
    except:
        return redirect(url_for('index')) 
@app.route('/download2')
def download_file2():        #to download TF file
    try:
        path = "Output/DE_miRNA.csv"
        return send_file(path , as_attachment=True)
    except:
        return redirect(url_for('index')) 
@app.route('/download3')
def download_file3():        #to download TF file
    try:
        path = "Output/miRNA_TF_mRNA_1.txt"
        return send_file(path , as_attachment=True)
    except:
        return redirect(url_for('index')) 
@app.route('/download4')
def download_file4():        #to download TF file
    try:
        path = "Output/miRNA_TF_mRNA_2.txt"
        return send_file(path , as_attachment=True)
    except:
        return redirect(url_for('index')) 

@app.route('/download5')
def download_file5():        #to download TF file
    try:
        path = "Output/Circuits.txt"
        return send_file(path , as_attachment=True)
    except:
        return redirect(url_for('index')) 

@app.route('/download6')
def download_file6():        #to download TF file
    try:
        path = "Output/Circuits_list.txt"
        return send_file(path , as_attachment=True)
    except:
        return redirect(url_for('index')) 

@app.route('/download7')
def download_file7():        #to download mRNA example file
    try:
        path = "example_download_mRNA.csv"
        return send_file(path , as_attachment=True)
    except:
        return redirect(url_for('index')) 

@app.route('/download8')
def download_file8():        #to download miRNA example file
    try:
        path = "example_download_miRNA.csv"
        return send_file(path , as_attachment=True)
    except:
        return redirect(url_for('index'))

if __name__ == '__main__':  # projenin baslatilmasi, portun acilmz
    import os
    HOST = os.environ.get('SERVER_HOST', 'localhost')
    try:
        PORT = int(os.environ.get('SERVER_PORT', '5555'))
    except ValueError:
        PORT = 5555
    app.run(HOST, PORT)
 

