from shiny import App, ui, render, reactive

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import json

sns.set_theme("notebook")

GLOBAL_GROUPS = json.load(open("global_groups.json"))
GENE_ACC = pd.read_table("./GENE_ACC.txt", sep="\t")
GENE_ACC["GeneID"] = GENE_ACC["GeneID"].astype(str).values
GSM2GRP = pd.read_table("./GSM2GRP.txt", sep="\t")

DATA = json.load(open("DATA.json"))

app_ui = ui.page_sidebar(
    ui.sidebar(
        ui.input_selectize(
            id="gse",
            label="Select a GEO dataset from the list below",
            choices=list(GLOBAL_GROUPS.keys()),
        ),
        ui.card(
            ui.card_header("Title"),
            ui.output_ui("title_contents")
        ),
        ui.card(
            ui.card_header("Summary"),
            ui.output_ui("summary_contents")
        ),
        ui.card(
            ui.card_header("Overall Design"),
            ui.output_ui("design_contents")
        ),
        width=500
    ),
    #ui.panel_title("GEO Browser"),
    ui.card(
        ui.card_header("Expression Levels (ELs) of Selected Probes"),
        ui.layout_columns(
        ui.output_ui("select_probe2plot"),
        ui.output_plot("box"),
        col_widths=[3, 9],
        height=400,
        fullscreen=True
        )
    ),
    ui.card(
        ui.card_header("Differential Expression Calculated by LIMMA"),
        ui.layout_columns(
            ui.output_ui("select_probe2tab"),
            ui.output_data_frame("limma_tt"),
            col_widths=[3, 9],
        )
    ),
    title="GEO Browser",
    fillable=True
)

def server(input, output, session):

    @render.text
    @reactive.event(input.gse, ignore_init=False)
    def title_contents():
        return DATA[input.gse()]["title"]
    
    @render.text
    @reactive.event(input.gse, ignore_init=False)
    def summary_contents():
        return DATA[input.gse()]["summary"]
    
    @render.text
    @reactive.event(input.gse, ignore_init=False)
    def design_contents():
        return DATA[input.gse()]["design"]

    @render.ui
    @reactive.event(input.gse, ignore_init=False)
    def select_probe2plot():
        genes = DATA[input.gse()]["f_tab"]["GeneName"]
        genes = sorted(genes.items(), key=lambda x: x[1])
        genes = dict(map(lambda x: (f"{input.gse()}_{x[0]}", f"{x[0]} ({x[1]})"), genes))
        widget = ui.input_selectize(id="probe2plot", label="Select probe(s) to display in the plot",
            choices=genes, multiple=True, selected=list(genes.keys())[0])
        return widget
    
    @render.plot
    @reactive.event(input.probe2plot, ignore_init=False, ignore_none=False)
    def box():
        p2p = ["_".join(p.split("_")[1:]) for p in input.probe2plot()]
        val = pd.DataFrame(DATA[input.gse()]["f_val"]).T.sort_index(axis=0)
        grp = GSM2GRP.loc[GSM2GRP["GEO.sample"].isin(val.index),:].sort_values("GEO.sample")["group.ID"].values
        val["grp"] = grp
        val = val.melt(id_vars="grp", value_vars=p2p)
        fig, ax = plt.subplots(layout="tight")
        sns.boxplot(data=val, ax=ax, x="grp", y="value", hue="variable", legend="brief", fill=True)
        ax.set_xlabel("")
        ax.set_xticks(sorted(list(set(grp))))
        ax.set_xticklabels([f"{g} (#{GLOBAL_GROUPS[input.gse()].index(g)})" for g in GLOBAL_GROUPS[input.gse()]],
            fontdict={"rotation": 45, "ha": "right", "fontsize": 10})
        ymin, ymax = ax.get_ylim()
        ydiff = (ymax - ymin) / 3
        ax.set_yticks(np.arange(ymin, ymax+ydiff, ydiff))
        ax.set_ylim((ymin, ymax))
        ax.set_yticklabels([f"{float(t.get_text()):.2f}" for t in ax.get_yticklabels()], fontdict={"fontsize": 10})
        ax.set_ylabel("Log2-Transformed EL", fontdict={"fontsize": 10})
        ax.get_legend().set_title("")
        return fig, ax
    
    @render.ui
    @reactive.event(input.gse, ignore_init=False)
    def select_probe2tab():
        genes = DATA[input.gse()]["f_tab"]["GeneName"]
        genes = sorted(genes.items(), key=lambda x: x[1])
        genes = dict(map(lambda x: (f"{input.gse()}_{x[0]}", f"{x[0]} ({x[1]})"), genes))
        widget = ui.input_selectize(id="probe2tab", label="Select probe(s) to display in the plot",
            choices=genes, multiple=True, selected=list(genes.keys()))
        return widget

    @render.data_frame
    @reactive.event(input.probe2tab, ignore_init=False)
    def limma_tt():
        p2t = ["_".join(p.split("_")[1:]) for p in input.probe2tab()]
        tt = pd.DataFrame(DATA[input.gse()]["f_tt"])
        tt = pd.concat((tt.index.to_series(), tt), axis=1)
        tt = tt.rename(columns={0: "probe"})
        tt = tt.loc[tt["probe"].isin(p2t),:].reindex(p2t)
        cnames = tt.columns.map(lambda x: x.replace(".grp", " vs #").replace("grp", "#"))
        tt = tt.rename(dict(zip(tt.columns, cnames)), axis="columns")
        for ii, cc in enumerate(tt.columns):
            if 0 < ii < 4:
                tt[cc] = tt[cc].map(lambda x: f"{x:.2g}")
            elif ii >= 4:
                tt[cc] = tt[cc].map(lambda x: f"{x:.2f}")
        return render.DataGrid(tt)
    
app = App(app_ui, server)