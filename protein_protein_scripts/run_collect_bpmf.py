"""
collect BPMF data
"""
from __future__ import print_function

import os
import glob
import argparse
import pickle

from _affinity_data import AffinityData
from _plots import scater_plot, plot_histogram

parser = argparse.ArgumentParser()
parser.add_argument("--postprocess_dir", type=str, default="postprocess")
parser.add_argument("--affinity_data_dir", type=str, default="affinity")
args = parser.parse_args()

BPMF_PKL = "bpmf.pkl"

AFFINITY_FILE_NAMES = ["affinity_v1.tsv",  "affinity_v2.tsv"]
affinity_data_files = [os.path.join(args.affinity_data_dir, file) for file in AFFINITY_FILE_NAMES]

experimental_dG = AffinityData(affinity_data_files).get_dG()

bpmf_files = glob.glob(os.path.join(args.postprocess_dir, "*", BPMF_PKL))
quantities_to_output = ["bpmf", "mean_Psi", "min_Psi"]

scores = {}
for file in bpmf_files:
    complex_id = file.split("/")[-2]
    results = pickle.load(open(file, "r"))

    for quantity in quantities_to_output:
        for phase in results[quantity].keys():
            score_id = quantity + "_" + phase
            if score_id not in scores.keys():
                scores[score_id] = {}
            scores[score_id][complex_id] = results[quantity][phase]


for score_id in scores.keys():
    out_string = "# complex    score\n"
    for complex_id in scores[score_id].keys():
        out_string += "%10s %15.5f\n"%(complex_id, scores[score_id][complex_id])
    open(score_id+".score", "w").write(out_string)

# histogram
plot_histogram(experimental_dG.values(), "experimental $\Delta G$ (kcal/mol)", "hist_experimental.pdf")

for score in scores.keys():
    data_to_plot = [v for v in scores[score].values() if str(v) not in ["nan", "inf", "-inf"]]
    plot_histogram(data_to_plot, score, "hist_"+score+".pdf")

# scatter plots

for score in scores.keys():
    complex_names = set(experimental_dG.keys()).intersection(scores[score].keys())
    complex_names = list(complex_names)

    exp_data = []
    cal_data = []
    for name in complex_names:
        if str(scores[score][name]) not in ["nan", "inf", "-inf"]:
            exp_data.append(experimental_dG[name])
            cal_data.append(scores[score][name])

    x_label = "experimental $\Delta G$ (kcal/mol)"
    y_label = score
    out = "scatter_"+score+".pdf"
    scater_plot(exp_data, cal_data, x_label, y_label, out, text_to_title=True)

