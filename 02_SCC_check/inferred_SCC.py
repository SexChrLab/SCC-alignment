import pandas as pd

summary_file = "indexcov-indexcov.ped"
summary_info = pd.read_csv(summary_file, sep="\t")

summary_info[["sample","coverage"]] = summary_info["sample_id"].str.split("_",expand=True)

summary_chrY = summary_info.loc[:,["sample","coverage","CNchrY"]]
summary_chrY["has_chrY"] = "no"
summary_chrY.loc[summary_chrY["CNchrY"] > 0.25, "has_chrY"] = "yes"

output_file = "inferred_SCC.csv"
summary_chrY.to_csv(output_file,index=False)
