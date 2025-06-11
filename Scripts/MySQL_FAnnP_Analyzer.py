import marimo

__generated_with = "0.13.6"
app = marimo.App()


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    # MySQL FAnnP Analyzer
    <center>To make sure you can use the tools to analyze the data. Follow the three steps in this document.<br>
    You can change the analysis to whatever you like.<br>
    The tools are used for examples of analysis and can be modified however you see fit.</center>
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ## Table of contents <a id="1"></a>
    * [Installation & Environment Activation](#2)
        * [Initialize Database Functions](#2.2)
        * [Initialize Table Functions](#2.3)
    * [Marimo Notebook Guide](#guide)
    * [Configuration](#3)
        * [Global Variables](#3.1)
        * [User Run List](#list_of_runs)
    * [Analysis Tools](#tools)
        * [Contig and Gene counts per Bin](#contig_gene_counts)
        * [Gene Searcher](#gene_searcher)
        * [Contig Viewer](#contig_viewer)
        * [Gene Viewer](#gene_viewer)
        * [NCBI NR](#diamond)
            * [NCBI NR taxonomy bar chart](#diamond_bar_chart)
            * [NCBI NR taxonomy counts](#diamond_taxonomy)
            * [NCBI NR taxonomy pie chart](#diamond_taxonomy_pie_chart)
        * [Cog](#cog)
            * [Cog counts per Cog Family](#cog_per_cog_family)
            * [Cog counts for Cog Family per Bin](#cog_for_cog_per_bin)
            * [Cog Family counts per Bin heatmap](#cog_family_bin_heatmap)
        * [Pfam](#pfam)
            * [Pfam Clans with the highest Pfam count](#highest_pfam_count_per_clan)
            * [Pfam Clan counts per Bin heatmap](#pfam_clan_bin_heatmap)
        * [KO](#ko)
            * [KO Jaccard Similarity Index for Bins](#ko_jaccard_bins)
            * [Pathway Modules coverage for Bins](#pathway_module_cov_bins)
            * [Pathway map KO coverage for Bins](#pathway_ko_map_bins)
        * [Benchmark Run Information](#benchmark_table)
        * [Benchmark Pie Chart](#benchmark_pie_chart)
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Installation & Environment Activation <a id="2"></a>""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ### Step 1: Initiates all the necessary functions to connect to the database. <a id="2.2"></a>
    #### Stores MySQL credentials, imports used libraries, and initiates database connection functions.
    """
    )
    return


@app.cell(hide_code=True)
def _():
    # MySQL credentials.
    mysql_username = "MySQL_Name_Placeholder" #"MySQL_Name_Placeholder"
    mysql_password = "MySQL_Password_Placeholder" #"MySQL_Password_Placeholder"
    mysql_database = "MySQL_Database_Placeholder" #"MySQL_Database_Placeholder"

    # Make sure the useful libraries we want to use are activated by running this code.
    import os
    import math
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    import mysql.connector
    from mysql.connector import errorcode
    import seaborn as sns
    import marimo as mo
    from pygenomeviz import GenomeViz
    from dna_features_viewer import GraphicFeature, GraphicRecord
    from Bio.KEGG import REST
    from Bio.KEGG.KGML import KGML_parser
    from Bio.Graphics.KGML_vis import KGMLCanvas
    from pdf2image import convert_from_path
    from PIL import Image as PIL_Image
    from io import BytesIO
    from io import StringIO
    from bokeh.plotting import figure, show
    from bokeh.io import output_notebook
    from bokeh.palettes import TolPRGn4, RdYlGn11, Spectral11, Set1_9, Set3_12, interp_palette
    from bokeh.transform import cumsum
    from bokeh.models import HoverTool

    # Checks if we can establish a connection to the FAnnP MySQL database.
    def connect_to_mysql():
        try:
            return mysql.connector.connect(user=mysql_username, password=mysql_password, host='localhost', database=mysql_database)
        except mysql.connector.Error as err:
            if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
                print("Something is wrong with your user name or password")
                return None
            elif err.errno == errorcode.ER_BAD_DB_ERROR:
                print("Database does not exist")
                return None
            else:
                print(err)
                return None

    # Tries to execute the given SQL query using the FAnnP MySQL database.
    def run_sql_query(query):
        cnx = connect_to_mysql()
        if cnx and cnx.is_connected():
            print("Connection with the database established.")
            with cnx.cursor() as cursor:
                cursor.execute(query)
                return [cursor.column_names, cursor.fetchall()]
            cnx.close()
            print("Connection with the database closed.")
        elif cnx:
            print("Error, no connection with the server established.")
            cnx.close()
    return (
        BytesIO,
        GenomeViz,
        GraphicFeature,
        GraphicRecord,
        KGMLCanvas,
        KGML_parser,
        PIL_Image,
        REST,
        Set1_9,
        Set3_12,
        Spectral11,
        StringIO,
        convert_from_path,
        cumsum,
        figure,
        interp_palette,
        math,
        mo,
        np,
        os,
        pd,
        plt,
        run_sql_query,
        sns,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ### Step 2: Initiates all the necessary functions to fetch specific tables from the database. <a id="2.3"></a>
    #### Initiates MySQL table functions, which will give all information for the designated table in a pandas dataframe.
    """
    )
    return


@app.cell(hide_code=True)
def _(pd, run_sql_query):
    # More useful functions to fetch data from run specific data/tables.
    # Information about the run.
    def get_run(run_id):
        columns, data = run_sql_query(f"SELECT * FROM Run;")
        annotation_data = pd.DataFrame(data, columns=columns)
        return annotation_data.sort_values("run_id").query(f"run_id == '{run_id}'")

    # Information about the bins.
    def get_bin(run_id):
        columns, data = run_sql_query(f"SELECT * FROM Bin;")
        annotation_data = pd.DataFrame(data, columns=columns)
        return annotation_data.sort_values("bin_id").query(f"run_id == '{run_id}'")

    # Information about metadata from the bins.
    def get_metadata(run_id):
        columns, data = run_sql_query(f"SELECT * FROM Metadata;")
        annotation_data = pd.DataFrame(data, columns=columns)
        return annotation_data.sort_values("bin_id").query(f"run_id == '{run_id}'")

    # Information about the contigs.
    def get_contig(run_id):
        columns, data = run_sql_query(f"SELECT * FROM Contig_run_{run_id};")
        annotation_data = pd.DataFrame(data, columns=columns)
        return annotation_data

    # Information about the genes.
    def get_feature(run_id):
        columns, data = run_sql_query(f"SELECT * FROM Feature_run_{run_id};")
        annotation_data = pd.DataFrame(data, columns=columns)
        return annotation_data

    # Information about the diamond annotations.
    def get_diamond(run_id):
        columns, data = run_sql_query(f"SELECT * FROM Diamond_Annotation_run_{run_id};")
        annotation_data = pd.DataFrame(data, columns=columns)
        return annotation_data

    # Information about the other annotation databases (e.g. arCOG, TIGR, CAZy, SignalP, MEROPs, TransporterDB, and HydDB).
    def get_annotation(run_id):
        columns, data = run_sql_query(f"SELECT * FROM Annotation;")
        annotation_data = pd.DataFrame(data, columns=columns)
        return annotation_data.sort_values("annotation_id").query(f"run_id == '{run_id}'")

    # Information about the other annotations (e.g. arCOG, TIGR, CAZy, SignalP, MEROPs, TransporterDB, and HydDB).
    def get_annotation_feature(run_id):
        columns, data = run_sql_query(f"SELECT * FROM Annotation_has_Feature_run_{run_id};")
        annotation_data = pd.DataFrame(data, columns=columns)
        return annotation_data

    # Information about the pfam annotations.
    def get_feature_pfam(run_id):
        columns, data = run_sql_query(f"SELECT * FROM Feature_has_Pfam_run_{run_id};")
        annotation_data = pd.DataFrame(data, columns=columns)
        return annotation_data

    # Information about the cog annotations.
    def get_feature_cog(run_id):
        columns, data = run_sql_query(f"SELECT * FROM Feature_has_Cog_run_{run_id};")
        annotation_data = pd.DataFrame(data, columns=columns)
        return annotation_data

    # Information about the ko annotations.
    def get_feature_ko(run_id):
        columns, data = run_sql_query(f"SELECT * FROM Feature_has_KO_run_{run_id};")
        annotation_data = pd.DataFrame(data, columns=columns)
        return annotation_data

    # Functions to fetch data from static tables and preloading them.
    # Tables that aren't specific to a run.
    # Pfam information table.
    def get_pfam():
        columns, data = run_sql_query(f"SELECT * FROM Pfam;")
        annotation_data = pd.DataFrame(data, columns=columns)
        return annotation_data
    # Preload table.
    global_pfam_dataframe = get_pfam()
    def get_pfam_dataframe():
        return global_pfam_dataframe

    # Pfam Clan information table.
    def get_pfam_clan():
        columns, data = run_sql_query(f"SELECT * FROM Pfam_clan;")
        annotation_data = pd.DataFrame(data, columns=columns)
        return annotation_data
    # Preload table.
    global_pfam_clan_dataframe = get_pfam_clan()
    def get_pfam_clan_dataframe():
        return global_pfam_clan_dataframe

    # Cog information table.
    def get_cog():
        columns, data = run_sql_query(f"SELECT * FROM Cog;")
        annotation_data = pd.DataFrame(data, columns=columns)
        return annotation_data
    # Preload table.
    global_cog_dataframe = get_cog()
    def get_cog_dataframe():
        return global_cog_dataframe

    # Cog - Cog Family information table.
    def get_cog_cog_family():
        columns, data = run_sql_query(f"SELECT * FROM Cog_has_Cog_family;")
        annotation_data = pd.DataFrame(data, columns=columns)
        return annotation_data
    # Preload table.
    global_cog_cog_family_dataframe = get_cog_cog_family()
    def get_cog_cog_family_dataframe():
        return global_cog_cog_family_dataframe

    # Cog Family information table.
    def get_cog_family():
        columns, data = run_sql_query(f"SELECT * FROM Cog_family;")
        annotation_data = pd.DataFrame(data, columns=columns)
        return annotation_data
    # Preload table.
    global_cog_family_dataframe = get_cog_family()
    def get_cog_family_dataframe():
        return global_cog_family_dataframe

    # KO information table.
    def get_ko():
        columns, data = run_sql_query(f"SELECT * FROM KO;")
        annotation_data = pd.DataFrame(data, columns=columns)
        return annotation_data
    # Preload table.
    global_ko_dataframe = get_ko()
    def get_ko_dataframe():
        return global_ko_dataframe

    # Module - KO information table.
    def get_module_ko():
        columns, data = run_sql_query(f"SELECT * FROM Module_has_KO;")
        annotation_data = pd.DataFrame(data, columns=columns)
        return annotation_data
    # Preload table.
    global_module_ko_dataframe = get_module_ko()
    def get_module_ko_dataframe():
        return global_module_ko_dataframe

    # Module information table.
    def get_module():
        columns, data = run_sql_query(f"SELECT * FROM Module;")
        annotation_data = pd.DataFrame(data, columns=columns)
        return annotation_data
    # Preload table.
    global_module_dataframe = get_module()
    def get_module_dataframe():
        return global_module_dataframe

    # Pathway - KO information table.
    def get_pathway_ko():
        columns, data = run_sql_query(f"SELECT * FROM Pathway_has_KO;")
        annotation_data = pd.DataFrame(data, columns=columns)
        return annotation_data
    # Preload table.
    global_pathway_ko_dataframe = get_pathway_ko()
    def get_pathway_ko_dataframe():
        return global_pathway_ko_dataframe

    # Pathway information table.
    def get_pathway():
        columns, data = run_sql_query(f"SELECT * FROM Pathway;")
        annotation_data = pd.DataFrame(data, columns=columns)
        return annotation_data
    # Preload table.
    global_pathway_dataframe = get_pathway()
    def get_pathway_dataframe():
        return global_pathway_dataframe

    # Pathway - Module information table.
    def get_pathway_module():
        columns, data = run_sql_query(f"SELECT * FROM Pathway_has_Module;")
        annotation_data = pd.DataFrame(data, columns=columns)
        return annotation_data
    # Preload table.
    global_pathway_module_dataframe = get_pathway_module()
    def get_pathway_module_dataframe():
        return global_pathway_module_dataframe

    # KO - EC information table.
    def get_ko_ec():
        columns, data = run_sql_query(f"SELECT * FROM KO_has_EC;")
        annotation_data = pd.DataFrame(data, columns=columns)
        return annotation_data
    # Preload table.
    global_ko_ec_dataframe = get_ko_ec()
    def get_ko_ec_dataframe():
        return global_ko_ec_dataframe

    # EC information table.
    def get_ec():
        columns, data = run_sql_query(f"SELECT * FROM EC;")
        annotation_data = pd.DataFrame(data, columns=columns)
        return annotation_data
    # Preload table.
    global_ec_dataframe = get_ec()
    def get_ec_dataframe():
        return global_ec_dataframe

    # Non run speficic functions.
    # Get all runs to find your own run ID.
    def get_run_all():
        columns, data = run_sql_query(f"SELECT * FROM Run;")
        annotation_data = pd.DataFrame(data, columns=columns)
        return annotation_data.sort_values("run_id")

    # Get all the bins for the list of given run IDs.
    def get_bin_multiple(run_id_list):
        columns, data = run_sql_query(f"SELECT * FROM Bin;")
        annotation_data = pd.DataFrame(data, columns=columns)
        return annotation_data[annotation_data["run_id"].isin(run_id_list)]

    # Return table from custom SQL query.
    # Check table names before using this function.
    def get_sql_table(sql_query):
        columns, data = run_sql_query(sql_query)
        annotation_data = pd.DataFrame(data, columns=columns)
        return annotation_data
    return (
        get_annotation,
        get_annotation_feature,
        get_bin,
        get_bin_multiple,
        get_cog_cog_family_dataframe,
        get_cog_dataframe,
        get_cog_family_dataframe,
        get_contig,
        get_diamond,
        get_feature,
        get_feature_cog,
        get_feature_ko,
        get_feature_pfam,
        get_ko_dataframe,
        get_ko_ec_dataframe,
        get_metadata,
        get_module_dataframe,
        get_pathway_dataframe,
        get_pathway_module_dataframe,
        get_pfam_clan_dataframe,
        get_pfam_dataframe,
        get_run,
        get_run_all,
        get_sql_table,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md("""## Marimo Notebook Guide <a id="guide"></a>""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        """
    ### Marimo Explanation
    __[Marimo](https://marimo.io/)__ is a tool that can open Marimo Notebooks, creating an interactive coding environment based on the python language similar to Jupyter Notebooks. The Notebook environment is stored in a python file which can be opened, read, and edited by the user.<br><br> 
    Unlike with Jupyter, Marimo executes the code in every cell automatically when a Notebook is opened. Marimo also has a few more useful features which aren't included with Jupyter. Marimo can download libraries, manage files, let you expore variables, and lets you configure several other things in the settings.<br><br> 
    Marimo lets you look at and change code which, is usually hidden for readability. If Marimo displays a Pandas table you will be able to search through it with the search icon. There is also a download button which might not separate the data correctly. A better way of downloading Pandas tables is with the *to_csv* function.
    ### Libraries
    Libraries are activated in [step one](#2.2) and more can be downloaded through Marimo or the conda environment.
    """
    )
    return


@app.cell(hide_code=True)
def _():
    """ Marimo Tutorial '''
    In the code cell under Step 2, you can see several function that are made for extracting information out of the FAnnP MySQL database. The information is set into a Pandas dataframe. You can apply Pandas functions to this Dataframe and download the table for example:
    -
    get_bin_dataframe().to_csv("bin_table.csv")
    get_diamond_dataframe().to_csv("diamond_table.tsv", sep='\t')


    Pandas Dataframes =
    To manipulate the dataframes, more coding exmaples can be found in the other analysis tools. Tables can be merged for example. You are able to set and reset the index of a dataframe and print specific columns. Pandas Documentation can be found here:
    https://pandas.pydata.org/docs/index.html
    -
    get_feature_dataframe().merge(get_diamond_dataframe(), on="feature_id")
    get_bin_dataframe().set_index(bin_data.columns[0]).reset_index()
    get_feature_pfam_dataframe()[["feature_id", "pfam_id"]]
    get_sql_table(f"SELECT * FROM Feature_run_{run_id};")       # MySQL query


    Marimo Printing =
    Marimo will attempt to automatically print the last line of code in a cell. Marimo can automatically print Pandas dataframes as well. Some objects like plots or markdown need a special Marimo function to print:
    -
    mo.md("Markdown Text.")                       # Marimo Markdown.
    mo.hstack([object]) or mo.vstack([object])     # Usually for images/UI elements.
    mo.ui.data_explorer(get_feature_dataframe())    # Data explorer for Pandas dataframes.


    Loose coding examples =
    -
    # Remove rows with empty values:
    get_contig_dataframe().replace('', np.nan).dropna()
    -
    # Get entries based on their index:
    get_contig_dataframe().set_index("contig_id").iloc[0:3]
    get_contig_dataframe().set_index("contig_id").loc["contig_1"]
    -
    # Get specific entries based on columns:
    get_contig_dataframe().query("bin_id == 'bin_1'")
    bin_list = ["bin_1", "bin_2"]
    get_contig_dataframe().query("bin_id in @bin_list")
    get_contig_dataframe().query("bin_id not in @bin_list")
    get_contig_dataframe()[get_contig_dataframe()["bin_id"].isin(["bin_1", "bin_2"])]
    get_contig_dataframe()[~get_contig_dataframe()["bin_id"].isin(["bin_1", "bin_2"])]
    -
    # Aggregate functions:
    get_contig_dataframe()["length"].agg(["sum", "mean", "max", "min", "count"])
    get_contig_dataframe().groupby("bin_id")["contig_id"].agg(Count="count")
    -
    # Sort table:
    get_contig_dataframe().sort_values(["bin_id", "contig_id"], ascending=False)
    -
    # Table transformations:
    get_contig_dataframe().transpose()

    """


    "If you press on the eye icon below, you will be able to see the code tutorial."
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Configuration <a id="3"></a>""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ### Step 3: Change global variables specific to your run. <a id="3.1"></a>
    #### Pick which run you want to access by selecting the run ID. The run ID will automatically be picked from the execution of the functional annotation pipeline.
    """
    )
    return


@app.cell(hide_code=True)
def _(get_run_all, mo):
    # These values are supposed to be replaced through the pipeline.
    current_user_name = "User_Name_Placeholder"
    current_run_id = "Run_ID_Placeholder"
    # All your run IDs.
    user_run_id_list = []
    # Here you can pick your run ID.
    # Prints all if pipeline placeholder values are still present.
    if "Placeholder" in current_user_name and "Placeholder" in current_run_id:
        user_run_id_list = get_run_all()["run_id"]
        run_id_dropdown = mo.ui.dropdown(options=user_run_id_list, label="Choose run ID:", value=user_run_id_list.values[0])
    else:
        user_run_id_list = get_run_all().query(f"user == '{current_user_name}'")["run_id"]
        run_id_dropdown = mo.ui.dropdown(options=user_run_id_list, label="Choose run ID:", value=current_run_id)
    mo.hstack([run_id_dropdown])
    return current_run_id, current_user_name, run_id_dropdown, user_run_id_list


@app.cell(hide_code=True)
def _(
    get_annotation,
    get_annotation_feature,
    get_bin,
    get_contig,
    get_diamond,
    get_feature,
    get_feature_cog,
    get_feature_ko,
    get_feature_pfam,
    get_metadata,
    get_run,
    run_id,
):
    # Dataframe version stored into memory. This will improve speed on large datasets.
    # Information about the run.
    global_run_dataframe = get_run(run_id)
    def get_run_dataframe():
        return global_run_dataframe

    # Information about the bins. Preloaded version.
    global_bin_dataframe = get_bin(run_id)
    def get_bin_dataframe():
        return global_bin_dataframe

    # Information about metadata from the bins. Preloaded version.
    global_metadata_dataframe = get_metadata(run_id)
    def get_metadata_dataframe():
        return global_metadata_dataframe

    # Information about the contigs. Preloaded version.
    global_contig_dataframe = get_contig(run_id)
    def get_contig_dataframe():
        return global_contig_dataframe

    # Information about the genes. Preloaded version.
    global_feature_dataframe = get_feature(run_id)
    def get_feature_dataframe():
        return global_feature_dataframe

    # Information about the diamond annotations. Preloaded version.
    global_diamond_dataframe = get_diamond(run_id)
    def get_diamond_dataframe():
        return global_diamond_dataframe

    # Information about the other annotation databases (e.g. arCOG, TIGR, CAZy, SignalP, MEROPs, TransporterDB, and HydDB). Preloaded version.
    global_annotation_dataframe = get_annotation(run_id)
    def get_annotation_dataframe():
        return global_annotation_dataframe

    # Information about the other annotations (e.g. arCOG, TIGR, CAZy, SignalP, MEROPs, TransporterDB, and HydDB). Preloaded version.
    global_annotation_feature_dataframe = get_annotation_feature(run_id)
    def get_annotation_feature_dataframe():
        return global_annotation_feature_dataframe

    # Information about the pfam annotations. Preloaded version.
    global_feature_pfam_dataframe = get_feature_pfam(run_id)
    def get_feature_pfam_dataframe():
        return global_feature_pfam_dataframe

    # Information about the cog annotations. Preloaded version.
    global_feature_cog_dataframe = get_feature_cog(run_id)
    def get_feature_cog_dataframe():
        return global_feature_cog_dataframe

    # Information about the ko annotations. Preloaded version.
    global_feature_ko_dataframe = get_feature_ko(run_id)
    def get_feature_ko_dataframe():
        return global_feature_ko_dataframe
    return (
        get_annotation_dataframe,
        get_annotation_feature_dataframe,
        get_bin_dataframe,
        get_contig_dataframe,
        get_diamond_dataframe,
        get_feature_cog_dataframe,
        get_feature_dataframe,
        get_feature_ko_dataframe,
        get_feature_pfam_dataframe,
        get_run_dataframe,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        """
    ### User Run List <a id="list_of_runs"></a>
    #### Prints a table displaying all the runs the user stores inside the database.
    """
    )
    return


@app.cell(hide_code=True)
def _(current_run_id, current_user_name, get_run_all, run_id_dropdown):
    # Here you can see the runs. The run ID will be stored in this variable and is used for the rest of the notebook.
    run_id = run_id_dropdown.value
    # Print run table, all entries if a user name isn't provided.
    if "Placeholder" in current_user_name and "Placeholder" in current_run_id:
        runs_df = get_run_all()
    else:
        runs_df = get_run_all().query(f"user == '{current_user_name}'")
    runs_df
    return (run_id,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Analysis Tools <a id="tools"></a>""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ### Contig and Gene counts per Bin <a id=contig_gene_counts></a>
    #### Prints the bin table along with the contig and gene counts within the bins.
    """
    )
    return


@app.cell(hide_code=True)
def _(get_bin_dataframe, get_contig_dataframe, get_feature_dataframe):
    # Count the gene amount and contig amount for each bin and then merge them into the bin table.
    contig_data = get_contig_dataframe().merge(get_bin_dataframe(), on="bin_id").groupby('bin_id')
    gene_data = get_feature_dataframe().merge(get_bin_dataframe(), on="bin_id").groupby('bin_id')
    get_bin_dataframe().merge(contig_data["bin_id"].agg(contig_amount="count"), on="bin_id").merge(gene_data["bin_id"].agg(gene_amount="count"), on="bin_id")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        """
    ### Gene Searcher <a id=gene_searcher></a>
    #### Pick a database and search for protein features and their descriptions.
    """
    )
    return


@app.cell(hide_code=True)
def _(mo, run_id, user_run_id_list):
    # If we should select more run IDs.
    all_selected_run_ids = mo.ui.multiselect(options=user_run_id_list, label="Choose run IDs:", value=[run_id])
    mo.hstack([all_selected_run_ids])
    return (all_selected_run_ids,)


@app.cell(hide_code=True)
def _(all_selected_run_ids, get_annotation_dataframe, get_bin_multiple, mo):
    # Which database type to pick.
    all_optional_databases = get_annotation_dataframe().set_index("database")["annotation_id"].to_dict()
    database_type_options = ["NCBI NR", "Pfam", "Cog", "Kegg"] + list(all_optional_databases.keys())
    database_type_dropdown = mo.ui.dropdown(options=database_type_options, label="Choose database:", value=database_type_options[0])
    # Which bins you want to pick.
    multiselect_bin_search_options = get_bin_multiple(all_selected_run_ids.value)["bin_id"].unique()
    bin_search_id_dropdown = mo.ui.multiselect(options=multiselect_bin_search_options, label="Choose bin IDs:", value=get_bin_multiple(all_selected_run_ids.value)["bin_id"].unique())
    # Text field with what feature should be searched for.
    feature_text_box = mo.ui.text(placeholder="Search...", label="Search for a feature:")
    # If we should account for letters being lower of upper case.
    case_sensitive_search = mo.ui.checkbox(label="Case sensitive search:", value=False)
    # If we should only view genes with multiple entries from the selected database..
    should_only_keep_duplicates = mo.ui.checkbox(label="Multiple gene entries only search:", value=False)
    mo.vstack([mo.hstack([database_type_dropdown]), mo.hstack([bin_search_id_dropdown]), mo.hstack([feature_text_box]), mo.hstack([case_sensitive_search]), mo.hstack([should_only_keep_duplicates])])
    return (
        all_optional_databases,
        bin_search_id_dropdown,
        case_sensitive_search,
        database_type_dropdown,
        feature_text_box,
        should_only_keep_duplicates,
    )


@app.cell(hide_code=True)
def _(
    all_optional_databases,
    all_selected_run_ids,
    bin_search_id_dropdown,
    case_sensitive_search,
    database_type_dropdown,
    feature_text_box,
    get_annotation_feature_dataframe,
    get_cog_dataframe,
    get_diamond_dataframe,
    get_feature_cog_dataframe,
    get_feature_dataframe,
    get_feature_ko_dataframe,
    get_feature_pfam_dataframe,
    get_ko_dataframe,
    get_pfam_dataframe,
    get_sql_table,
    pd,
    run_id,
    should_only_keep_duplicates,
):
    # Looks through the accesion and discription of the selected database to search for an gene annotation entry.
    # If the var is not replaced by a table, it couldn't find an entry.
    feature_results_table = f"No features with the name '{feature_text_box.value}' found."
    a = 0
    for run_id_entry in all_selected_run_ids.value:
        match database_type_dropdown.value:
            # NCBI database.
            case "NCBI NR":
                if run_id_entry == run_id:
                    feature_table_information_dataframe = get_feature_dataframe().merge(get_diamond_dataframe(), on="feature_id")
                else:
                    feature_table_information_dataframe = get_sql_table(f"SELECT Feature_run_{run_id_entry}.feature_id, contig_id, bin_id, name, c_start, c_end, direction, f_type, definition, feature_annotation_id, acc, description, tax_id, taxonomy, query_start, query_end, subject_start, subject_end, evalue, bitscore, identity FROM Feature_run_{run_id_entry} INNER JOIN Diamond_Annotation_run_{run_id_entry} ON Feature_run_{run_id_entry}.feature_id = Diamond_Annotation_run_{run_id_entry}.feature_id;")
                feature_table_information_dataframe["regex_string"] = feature_table_information_dataframe["acc"] + " " + feature_table_information_dataframe["description"] + " " + feature_table_information_dataframe["taxonomy"]
                feature_table_information_dataframe = feature_table_information_dataframe[feature_table_information_dataframe["regex_string"].str.contains(feature_text_box.value, case=case_sensitive_search.value)]
                # Remove genes that only have a single feature entry.
                if should_only_keep_duplicates.value:
                    feature_table_information_dataframe = feature_table_information_dataframe[feature_table_information_dataframe["feature_id"].duplicated(keep=False)]
                # If we found results, make sure to add it to the table we are going to print.
                if not feature_table_information_dataframe.empty:
                    if isinstance(feature_results_table, pd.DataFrame):
                        feature_results_table = pd.concat([feature_results_table, feature_table_information_dataframe])
                    else:
                        feature_results_table = feature_table_information_dataframe
            # Pfam database.
            case "Pfam":
                if run_id_entry == run_id:
                    feature_table_information_dataframe = get_feature_dataframe().merge(get_feature_pfam_dataframe(), on="feature_id").merge(get_pfam_dataframe(), on="pfam_id")
                else:
                    feature_table_information_dataframe = get_sql_table(f"SELECT Feature_run_{run_id_entry}.feature_id, contig_id, bin_id, name, c_start, c_end, direction, f_type, definition, Feature_has_Pfam_run_{run_id_entry}.pfam_id, start, end, evalue, score, pfam_acc, pfam_desc FROM Feature_run_{run_id_entry} INNER JOIN Feature_has_Pfam_run_{run_id_entry} ON Feature_run_{run_id_entry}.feature_id = Feature_has_Pfam_run_{run_id_entry}.feature_id INNER JOIN Pfam ON Feature_has_Pfam_run_{run_id_entry}.pfam_id = Pfam.pfam_id;")
                feature_table_information_dataframe["regex_string"] = feature_table_information_dataframe["pfam_id"] + " " + feature_table_information_dataframe["pfam_acc"] + " " + feature_table_information_dataframe["pfam_desc"]
                feature_table_information_dataframe = feature_table_information_dataframe[feature_table_information_dataframe["regex_string"].str.contains(feature_text_box.value, case=case_sensitive_search.value)]
                # Remove genes that only have a single feature entry.
                if should_only_keep_duplicates.value:
                    feature_table_information_dataframe = feature_table_information_dataframe[feature_table_information_dataframe["feature_id"].duplicated(keep=False)]
                # If we found results, make sure to add it to the table we are going to print.
                if not feature_table_information_dataframe.empty:
                    if isinstance(feature_results_table, pd.DataFrame):
                        feature_results_table = pd.concat([feature_results_table, feature_table_information_dataframe])
                    else:
                        feature_results_table = feature_table_information_dataframe
            # Cog database.
            case "Cog":
                if run_id_entry == run_id:
                    feature_table_information_dataframe = get_feature_dataframe().merge(get_feature_cog_dataframe(), on="feature_id").merge(get_cog_dataframe(), on="cog_id")
                else:
                    feature_table_information_dataframe = get_sql_table(f"SELECT Feature_run_{run_id_entry}.feature_id, contig_id, bin_id, name, c_start, c_end, direction, f_type, definition, Feature_has_Cog_run_{run_id_entry}.cog_id, start, end, evalue, score, cog_desc FROM Feature_run_{run_id_entry} INNER JOIN Feature_has_Cog_run_{run_id_entry} ON Feature_run_{run_id_entry}.feature_id = Feature_has_Cog_run_{run_id_entry}.feature_id INNER JOIN Cog ON Feature_has_Cog_run_{run_id_entry}.cog_id = Cog.cog_id;")
                feature_table_information_dataframe["regex_string"] = feature_table_information_dataframe["cog_id"] + " " + feature_table_information_dataframe["cog_desc"]
                feature_table_information_dataframe = feature_table_information_dataframe[feature_table_information_dataframe["regex_string"].str.contains(feature_text_box.value, case=case_sensitive_search.value)]
                # Remove genes that only have a single feature entry.
                if should_only_keep_duplicates.value:
                    feature_table_information_dataframe = feature_table_information_dataframe[feature_table_information_dataframe["feature_id"].duplicated(keep=False)]
                # If we found results, make sure to add it to the table we are going to print.
                if not feature_table_information_dataframe.empty:
                    if isinstance(feature_results_table, pd.DataFrame):
                        feature_results_table = pd.concat([feature_results_table, feature_table_information_dataframe])
                    else:
                        feature_results_table = feature_table_information_dataframe
            # Kegg database.
            case "Kegg":
                if run_id_entry == run_id:
                    feature_table_information_dataframe = get_feature_dataframe().merge(get_feature_ko_dataframe(), on="feature_id").merge(get_ko_dataframe(), on="ko_id")
                else:
                    feature_table_information_dataframe = get_sql_table(f"SELECT Feature_run_{run_id_entry}.feature_id, contig_id, bin_id, name, c_start, c_end, direction, f_type, definition, Feature_has_KO_run_{run_id_entry}.ko_id, start, end, evalue, score, ko_symbol, ko_desc FROM Feature_run_{run_id_entry} INNER JOIN Feature_has_KO_run_{run_id_entry} ON Feature_run_{run_id_entry}.feature_id = Feature_has_KO_run_{run_id_entry}.feature_id INNER JOIN KO ON Feature_has_KO_run_{run_id_entry}.ko_id = KO.ko_id;")
                feature_table_information_dataframe["regex_string"] = feature_table_information_dataframe["ko_id"] + " " + feature_table_information_dataframe["ko_symbol"] + " " + feature_table_information_dataframe["ko_desc"]
                feature_table_information_dataframe = feature_table_information_dataframe[feature_table_information_dataframe["regex_string"].str.contains(feature_text_box.value, case=case_sensitive_search.value)]
                # Remove genes that only have a single feature entry.
                if should_only_keep_duplicates.value:
                    feature_table_information_dataframe = feature_table_information_dataframe[feature_table_information_dataframe["feature_id"].duplicated(keep=False)]
                # If we found results, make sure to add it to the table we are going to print.
                if not feature_table_information_dataframe.empty:
                    if isinstance(feature_results_table, pd.DataFrame):
                        feature_results_table = pd.concat([feature_results_table, feature_table_information_dataframe])
                    else:
                        feature_results_table = feature_table_information_dataframe
            # The rest of the databases.
            case _:
                if run_id_entry == run_id:
                    feature_table_information_dataframe = get_feature_dataframe().merge(get_annotation_feature_dataframe(), on="feature_id").query(f"annotation_id == {all_optional_databases[database_type_dropdown.value]}")
                else:
                    feature_table_information_dataframe = get_sql_table(f"SELECT Feature_run_{run_id_entry}.feature_id, contig_id, bin_id, name, c_start, c_end, direction, f_type, definition, annotation_id, start, end, acc, description, evalue, score, identity FROM Feature_run_{run_id_entry} INNER JOIN Annotation_has_Feature_run_{run_id_entry} ON Feature_run_{run_id_entry}.feature_id = Annotation_has_Feature_run_{run_id_entry}.feature_id WHERE annotation_id='{all_optional_databases[database_type_dropdown.value]}';")
                if database_type_dropdown.value == "SignalP":
                    feature_table_information_dataframe["regex_string"] = feature_table_information_dataframe["acc"]
                else:
                    feature_table_information_dataframe["regex_string"] = feature_table_information_dataframe["acc"] + " " + feature_table_information_dataframe["description"]
                feature_table_information_dataframe = feature_table_information_dataframe[feature_table_information_dataframe["regex_string"].str.contains(feature_text_box.value, case=case_sensitive_search.value, na=False)]
                # Remove genes that only have a single feature entry.
                if should_only_keep_duplicates.value:
                    feature_table_information_dataframe = feature_table_information_dataframe[feature_table_information_dataframe["feature_id"].duplicated(keep=False)]
                # If we found results, make sure to add it to the table we are going to print.
                if not feature_table_information_dataframe.empty:
                    if isinstance(feature_results_table, pd.DataFrame):
                        feature_results_table = pd.concat([feature_results_table, feature_table_information_dataframe])
                    else:
                        feature_results_table = feature_table_information_dataframe
    if isinstance(feature_results_table, pd.DataFrame):
        feature_results_table = feature_results_table[feature_results_table["bin_id"].isin(bin_search_id_dropdown.value)]
    feature_results_table
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        """
    ### Contig Viewer <a id=contig_viewer></a>
    #### Display all genes on the selected contigs.
    """
    )
    return


@app.cell(hide_code=True)
def _(get_bin_dataframe, mo):
    # Here you can pick your bin IDs.
    bin_id_multiselect = mo.ui.multiselect(options=get_bin_dataframe()["bin_id"], label="Choose bin ID:", value=[get_bin_dataframe()["bin_id"].iloc[0]])
    mo.hstack([bin_id_multiselect])
    return (bin_id_multiselect,)


@app.cell(hide_code=True)
def _(bin_id_multiselect, get_contig_dataframe, mo):
    # Here you can pick contig IDs.
    multiselect_contig_options = get_contig_dataframe()[get_contig_dataframe()["bin_id"].isin(bin_id_multiselect.value)]["contig_id"].values
    contig_id_dropdown = mo.ui.multiselect(options=multiselect_contig_options, label="Choose contig ID:", value=[get_contig_dataframe()[get_contig_dataframe()["bin_id"].isin(bin_id_multiselect.value)]["contig_id"].iloc[0]])
    # Pick the range of the contig to display.
    min_max_print_number_contig = mo.ui.range_slider(start=0, stop=100, value=[0, 100], label="Contig Range:")
    # If we should print the whole contigs (might take long depending on the contig sizes).
    should_print_contig_image = mo.ui.checkbox(label="Print first 100 genes only:", value=True)
    mo.vstack([mo.hstack([contig_id_dropdown]), mo.hstack([min_max_print_number_contig]), mo.hstack([should_print_contig_image])])
    return (
        contig_id_dropdown,
        min_max_print_number_contig,
        should_print_contig_image,
    )


@app.cell(hide_code=True)
def _(
    GenomeViz,
    contig_id_dropdown,
    get_contig_dataframe,
    get_feature_dataframe,
    min_max_print_number_contig,
    should_print_contig_image,
):
    # Print contig after gathering gene locations.
    gv_contig = GenomeViz()
    for gv_gene_target_entry in contig_id_dropdown.value:
        contig_length = get_contig_dataframe().query(f"contig_id == '{gv_gene_target_entry}'")["length"].iloc[0]
        contig_length_min = round(contig_length*(min_max_print_number_contig.value[0]/100))
        contig_length_max = round(contig_length*(min_max_print_number_contig.value[1]/100))
        contig_track = gv_contig.add_feature_track(gv_gene_target_entry, (contig_length_min, contig_length_max))
        contig_track.add_sublabel()
        for contig_num, single_contig in enumerate(get_feature_dataframe().query(f"contig_id == '{gv_gene_target_entry}'").to_numpy()):
            if should_print_contig_image.value:
                if contig_num > 99:
                    break
            if single_contig[4] > contig_length_min and single_contig[5] < contig_length_max:
                if single_contig[6] == "-":
                    contig_track.add_feature(single_contig[4], single_contig[5], 1, label=single_contig[0])
                else:
                    contig_track.add_feature(single_contig[4], single_contig[5], -1, label=single_contig[0])
    gv_contig.plotfig()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        """
    ### Gene Viewer <a id=gene_viewer></a>
    #### Display all features on the selected genes.
    """
    )
    return


@app.cell(hide_code=True)
def _(contig_id_dropdown, get_feature_dataframe, mo):
    # Here you can pick your contig ID.
    multiselect_gene_options = get_feature_dataframe()[get_feature_dataframe()["contig_id"].isin(contig_id_dropdown.value)]["feature_id"].values
    gene_id_dropdown = mo.ui.multiselect(options=multiselect_gene_options, label="Choose gene ID:", value=[get_feature_dataframe()[get_feature_dataframe()["contig_id"].isin(contig_id_dropdown.value)]["feature_id"].iloc[0]])
    # Pick the range of the gene to display.
    min_max_print_number_gene = mo.ui.range_slider(start=0, stop=100, value=[0, 100], label="Gene Range:")
    # If we should print the pathway next to the original pathway.
    should_print_gene_name = mo.ui.checkbox(label="Print the gene name:", value=True)
    mo.vstack([mo.hstack([gene_id_dropdown]), mo.hstack([min_max_print_number_gene]), mo.hstack([should_print_gene_name])])
    return gene_id_dropdown, min_max_print_number_gene, should_print_gene_name


@app.cell(hide_code=True)
def _(
    GraphicFeature,
    GraphicRecord,
    gene_id_dropdown,
    get_annotation_feature_dataframe,
    get_diamond_dataframe,
    get_feature_cog_dataframe,
    get_feature_dataframe,
    get_feature_ko_dataframe,
    get_feature_pfam_dataframe,
    min_max_print_number_gene,
    should_print_gene_name,
):
    # List of pictures containing the features for each gene.
    gene_features_record_crop_list = []
    for gene_entry_var in gene_id_dropdown.value:
        # Print gene after gathering feature locations.
        target_gene_name = get_feature_dataframe().query(f"feature_id == '{gene_entry_var}'")
        gene_features_list = []
        # Fetch gene length.
        gene_length = get_feature_dataframe().query(f"feature_id == '{gene_entry_var}'")["c_end"].iloc[0] - get_feature_dataframe().query(f"feature_id == '{gene_entry_var}'")["c_start"].iloc[0]
        # Feature annotation table.
        for single_gene in target_gene_name.merge(get_annotation_feature_dataframe(), on="feature_id").to_numpy():
            gene_features_list.append(GraphicFeature(start=single_gene[10]*3, end=single_gene[11]*3, color="#ffe7b0", label=single_gene[12]))
        # Feature Pfam table.
        for single_gene in target_gene_name.merge(get_feature_pfam_dataframe(), on="feature_id").to_numpy():
            gene_features_list.append(GraphicFeature(start=single_gene[10]*3, end=single_gene[11]*3, color="#b7ffb7", label=single_gene[9]))
        # Feature Cog table.
        for single_gene in target_gene_name.merge(get_feature_cog_dataframe(), on="feature_id").to_numpy():
            gene_features_list.append(GraphicFeature(start=single_gene[10]*3, end=single_gene[11]*3, color="#ffbfbf", label=single_gene[9]))
        # Feature KO table.
        for single_gene in target_gene_name.merge(get_feature_ko_dataframe(), on="feature_id").to_numpy():
            gene_features_list.append(GraphicFeature(start=single_gene[10]*3, end=single_gene[11]*3, color="#c7c7ff", label=single_gene[9]))
        # Feature Diamond table.
        for single_gene in target_gene_name.merge(get_diamond_dataframe(), on="feature_id").to_numpy():
            gene_features_list.append(GraphicFeature(start=single_gene[14]*3, end=single_gene[15]*3, color="#ffd3be", label=single_gene[10]))
        gene_features_record = GraphicRecord(sequence_length=gene_length, features=gene_features_list)
        gene_features_record_crop = gene_features_record.crop((gene_length*(min_max_print_number_gene.value[0]/100), gene_length*(min_max_print_number_gene.value[1]/100)))
        if should_print_gene_name.value:
            gene_features_record_crop_list.append(gene_features_record_crop.plot()[0].set_title(gene_entry_var, loc='left'))
        else:
            gene_features_record_crop_list.append(gene_features_record_crop.plot()[0])
    gene_features_record_crop_list
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        """
    ## NCBI NR <a id="diamond"></a>
    The taxonomy of the NCBI NR entries can have multiple taxonomies. Some genes will have multiple taxonomies since in the NCBI NR database, genes are shared when they are samilar enough. This can cause there to be a larger representation of the taxonomy than the genes.
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        """
    ### NCBI NR taxonomy bar chart <a id="diamond_bar_chart"></a>
    #### Prints a bar chart displaying the taxonomy of each bin based on the selected taxonomy level.
    """
    )
    return


@app.cell(hide_code=True)
def _(get_feature_dataframe, mo, multiselect_bin_diamond_options):
    # Which bins you want to display.
    bin_diamond_id_multiselect_table = mo.ui.multiselect(options=multiselect_bin_diamond_options, label="Choose bin IDs:", value=get_feature_dataframe()["bin_id"].unique())
    # If we should only view genes with multiple entries from the selected database..
    print_taxonomy_twenty = mo.ui.checkbox(label="Print in mutiple barplots for 40 or more entries:", value=True)
    # What taxonomy level you want to select for displaying the bar plot.
    diamond_tax_level_list_table_options = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
    diamond_tax_level_table_dropdown = mo.ui.dropdown(options=diamond_tax_level_list_table_options, label="Choose taxonomy level for legend:", value=diamond_tax_level_list_table_options[0])
    # How many unique taxonomy values we should print. Too many might take a while to print.
    index_bar_plot_size_number = mo.ui.number(start=5, stop=100, label="Amount of unique taxonomy entries to display:", value=10)
    mo.vstack([mo.hstack([bin_diamond_id_multiselect_table]), mo.hstack([print_taxonomy_twenty]), mo.hstack([diamond_tax_level_table_dropdown]), mo.hstack([index_bar_plot_size_number])])
    return (
        bin_diamond_id_multiselect_table,
        diamond_tax_level_table_dropdown,
        index_bar_plot_size_number,
        print_taxonomy_twenty,
    )


@app.cell(hide_code=True)
def _(
    Set1_9,
    Set3_12,
    bin_diamond_id_multiselect_table,
    diamond_tax_level_table_dropdown,
    get_diamond_dataframe,
    get_feature_dataframe,
    index_bar_plot_size_number,
    mo,
    np,
    pd,
    print_taxonomy_twenty,
    sns,
):
    # Create a stacked bar plot using a certain taxonomy level and displays it for each bin in percentages.
    # Split entries at this amount.
    split_at_tax_amount = 40
    unclassified_tax_string = "Unclassified;Unclassified;Unclassified;Unclassified;Unclassified;Unclassified;Unclassified;Unclassified"
    diamond_tax_level_list_bar = ["Superkingdom", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
    all_feature_tax_graph_df = get_diamond_dataframe().replace('', unclassified_tax_string).dropna().merge(get_feature_dataframe(), on="feature_id").set_index(["bin_id"])["taxonomy"].apply(lambda x: x.split(',')).explode().reset_index()
    all_feature_tax_graph_df["taxonomy"] = all_feature_tax_graph_df["taxonomy"].str.split(";")
    all_feature_tax_graph_df["taxonomy"] = pd.DataFrame(all_feature_tax_graph_df.taxonomy.tolist(), columns=diamond_tax_level_list_bar)[f"{diamond_tax_level_table_dropdown.value}"]
    all_feature_tax_graph_df = all_feature_tax_graph_df = all_feature_tax_graph_df.replace('none', "Unclassified").dropna()
    all_feature_tax_graph_df = all_feature_tax_graph_df[all_feature_tax_graph_df["bin_id"].isin(bin_diamond_id_multiselect_table.value)]
    # Check which taxonomic values to keep and which to turn into Other Classified values.
    top_taxonomy_listed = all_feature_tax_graph_df.groupby(["taxonomy"])["taxonomy"].agg(Count="count").sort_values(["Count"], ascending=False).index[:index_bar_plot_size_number.value-1]
    all_feature_tax_graph_df.loc[~all_feature_tax_graph_df["taxonomy"].isin(top_taxonomy_listed), "taxonomy"] = "Other Classified"
    # Make count table for the plot.
    all_feature_tax_graph_df = all_feature_tax_graph_df.groupby(["taxonomy", "bin_id"])["taxonomy"].agg(Count="count").sort_values(["Count"], ascending=False).reset_index()
    all_feature_tax_graph_df = all_feature_tax_graph_df.pivot_table(values='Count', index='bin_id', columns='taxonomy').replace(np.nan, 0)
    all_feature_tax_graph_df = all_feature_tax_graph_df[all_feature_tax_graph_df.sum().sort_values(ascending=False).index].transpose()
    all_feature_tax_graph_df = all_feature_tax_graph_df/all_feature_tax_graph_df.sum()
    all_feature_tax_graph_df = all_feature_tax_graph_df.transpose()
    # Make plot and print it.
    plot_list = []
    if print_taxonomy_twenty.value and len(all_feature_tax_graph_df) > split_at_tax_amount:
        for split_graph_i in range(int(len(all_feature_tax_graph_df)/split_at_tax_amount)+1):
            split_graph_i *= split_at_tax_amount
            sns.set(style='white')
            bar_plot = all_feature_tax_graph_df.iloc[split_graph_i:split_graph_i+split_at_tax_amount].plot(kind='bar', stacked=True, width=0.9, color=Set1_9+Set3_12)
            sns.move_legend(bar_plot, "upper left", bbox_to_anchor=(1, 1))
            plot_list.append(bar_plot)
    else:
        sns.set(style='white')
        bar_plot = all_feature_tax_graph_df.plot(kind='bar', stacked=True, width=0.9, color=Set1_9+Set3_12)
        sns.move_legend(bar_plot, "upper left", bbox_to_anchor=(1, 1))
        plot_list = bar_plot
    mo.vstack([plot_list])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ### NCBI NR taxonomy counts <a id="diamond_taxonomy"></a>
    #### Prints a table displaying the highest taxonomic counts of the selected bins and taxonomy level.
    """
    )
    return


@app.cell(hide_code=True)
def _(get_bin_dataframe, get_feature_dataframe, mo):
    # What taxonomy level you want to select.
    diamond_tax_level_dict = {"Superkingdom": 0, "Kingdom": 1, "Phylum": 2, "Class": 3, "Order": 4, "Family": 5, "Genus": 6, "Species": 7}
    diamond_tax_level_dropdown = mo.ui.dropdown(options=diamond_tax_level_dict.keys(), label="Choose taxonomy level:", value=list(diamond_tax_level_dict.keys())[-1])
    # Which bins you want to pick.
    multiselect_bin_diamond_options = get_bin_dataframe()["bin_id"].unique()
    bin_diamond_id_multiselect = mo.ui.multiselect(options=multiselect_bin_diamond_options, label="Choose bin IDs:", value=get_feature_dataframe()["bin_id"].unique())
    # If we should print the highest taxonomy count only for each bin to make the table easier to read.
    print_highest_bins_tax = mo.ui.checkbox(label="Print highest taxonomy count per bin:", value=False)
    mo.vstack([mo.hstack([diamond_tax_level_dropdown]), mo.hstack([bin_diamond_id_multiselect]), mo.hstack([print_highest_bins_tax])])
    return (
        bin_diamond_id_multiselect,
        diamond_tax_level_dict,
        diamond_tax_level_dropdown,
        multiselect_bin_diamond_options,
        print_highest_bins_tax,
    )


@app.cell(hide_code=True)
def _(
    bin_diamond_id_multiselect,
    diamond_tax_level_dict,
    diamond_tax_level_dropdown,
    get_diamond_dataframe,
    get_feature_dataframe,
    np,
    pd,
    print_highest_bins_tax,
):
    # Only prints the first four entries to keep things conpact.
    all_feature_tax_df = get_diamond_dataframe().replace('', np.nan).dropna().merge(get_feature_dataframe(), on="feature_id").set_index(["bin_id"])["taxonomy"].apply(lambda x: x.split(',')).explode().reset_index()
    all_feature_tax_df["taxonomy"] = all_feature_tax_df["taxonomy"].str.split(";")
    all_feature_tax_df["taxonomy"] = pd.DataFrame(all_feature_tax_df.taxonomy.tolist(), columns=diamond_tax_level_dict.keys())[f"{diamond_tax_level_dropdown.value}"]
    all_feature_tax_df = all_feature_tax_df.replace('none', np.nan).dropna().groupby(["taxonomy", "bin_id"])["taxonomy"].agg(Count="count").sort_values(["Count"], ascending=False).reset_index()
    if print_highest_bins_tax.value:
        all_feature_tax_df = all_feature_tax_df.sort_values("Count", ascending=False).drop_duplicates("bin_id")
    all_feature_tax_df[all_feature_tax_df["bin_id"].isin(bin_diamond_id_multiselect.value)]
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        """
    ### NCBI NR taxonomy pie chart <a id="diamond_taxonomy_pie_chart"></a>
    #### Prints a pie chart displaying eight different taxonomy levels for the selected bin.
    """
    )
    return


@app.cell(hide_code=True)
def _(get_bin_dataframe, get_feature_dataframe, mo):
    # Which bin you want to make a pie chart out of.
    dropdown_bin_diamond_options = get_bin_dataframe()["bin_id"].unique()
    bin_diamond_id_dropdown = mo.ui.dropdown(options=dropdown_bin_diamond_options, label="Choose bin ID:", value=get_feature_dataframe()["bin_id"].iloc[0])
    # What taxonomy level you want to select for displaying its legend.
    diamond_tax_level_list_options = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
    diamond_tax_level_legend_dropdown = mo.ui.dropdown(options=diamond_tax_level_list_options, label="Choose taxonomy level for legend:")
    # If we should display taxonomy values that have less x amount of entries.
    taxonomy_required_minimum = mo.ui.number(start=1, stop=1000, label="Amount of entries needed to display the taxonomy on the pie chart:", value=1)
    # The pie chart relation to the next layer isn't displayed accurately without using the default settings. This is because the same taxonomy ratios are being sorted but stops being acurate when none values merge or taxonomy entries are being hidden.
    should_hide_none_values = mo.ui.radio(options=["Show none values:", "Hide none values:", "Use parent values of none values:"], value="Use parent values of none values:")
    mo.vstack([mo.hstack([bin_diamond_id_dropdown]), mo.hstack([diamond_tax_level_legend_dropdown]), mo.hstack([taxonomy_required_minimum]), mo.hstack([should_hide_none_values])])
    return (
        bin_diamond_id_dropdown,
        diamond_tax_level_legend_dropdown,
        should_hide_none_values,
        taxonomy_required_minimum,
    )


@app.cell(hide_code=True)
def _(
    Spectral11,
    bin_diamond_id_dropdown,
    cumsum,
    diamond_tax_level_legend_dropdown,
    figure,
    get_diamond_dataframe,
    get_feature_dataframe,
    interp_palette,
    math,
    mo,
    np,
    pd,
    should_hide_none_values,
    taxonomy_required_minimum,
):
    # Makes a pie chart from taxonomy of the genes found in the NCBI database.
    diamond_tax_level_list = ["Superkingdom", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
    all_feature_tax_pie_df = get_diamond_dataframe().replace('', np.nan).dropna().merge(get_feature_dataframe(), on="feature_id").query(f"bin_id == '{bin_diamond_id_dropdown.value}'").set_index(["bin_id"])["taxonomy"].apply(lambda x: x.split(',')).explode().reset_index()
    all_feature_tax_pie_df["taxonomy"] = all_feature_tax_pie_df["taxonomy"].str.split(";")
    all_feature_tax_pie_df = pd.DataFrame(all_feature_tax_pie_df.taxonomy.tolist(), columns=diamond_tax_level_list).sort_values(diamond_tax_level_list, ascending=False)
    # If a taxonomy value is none, use the next available taxonomy from a higher taxonomy level.
    if should_hide_none_values.value == "Use parent values of none values:":
        for tax_index_val in range(len(diamond_tax_level_list) - 1):
            all_feature_tax_pie_df.loc[all_feature_tax_pie_df[diamond_tax_level_list[tax_index_val+1]] == 'none', diamond_tax_level_list[tax_index_val+1]] = all_feature_tax_pie_df[diamond_tax_level_list[tax_index_val]]
    # Uses a different layout if a taxonomy level got selected.
    if diamond_tax_level_legend_dropdown.value:
        tax_pie_chart = figure(width=600, height=400, title=f"{diamond_tax_level_legend_dropdown.value} Taxonomy Pie Chart - - - {bin_diamond_id_dropdown.value}", tooltips="@taxonomy: @count", x_range=(-2, 1))
    else:
        tax_pie_chart = figure(width=600, height=600, title=f"Taxonomy Pie Chart - - - {bin_diamond_id_dropdown.value}", tooltips="@taxonomy: @count")

    # Loop through each taxonomy level to print them in a pie chart.
    for diamond_level_num, diamond_tax_level in enumerate(diamond_tax_level_list):
        single_feature_tax_pie_df = all_feature_tax_pie_df.groupby([f"{diamond_tax_level}"], sort=False)[diamond_tax_level].agg(count="count").reset_index().rename(columns = {diamond_tax_level: 'taxonomy'})
        single_feature_tax_pie_df = single_feature_tax_pie_df.query(f"count >= {taxonomy_required_minimum.value}")
        # Hides the none values present in each taxonomy level.
        if should_hide_none_values.value == "Hide none values:":
            single_feature_tax_pie_df = single_feature_tax_pie_df.query(f"taxonomy != 'none'")
        single_feature_tax_pie_df["angle"] = single_feature_tax_pie_df["count"]/single_feature_tax_pie_df["count"].sum() * 2 * math.pi
        # If a taxonomy level is selected, paint the level and the rest grey.
        if diamond_tax_level == diamond_tax_level_legend_dropdown.value:
            single_feature_tax_pie_df["color"] = interp_palette(Spectral11, len(single_feature_tax_pie_df))
        # Use colours if no taxonomy level is selected.
        elif not diamond_tax_level_legend_dropdown.value:
            single_feature_tax_pie_df["color"] = interp_palette(Spectral11, len(single_feature_tax_pie_df))
        # Rest are painted grey here.
        else:
            single_feature_tax_pie_df["color"] = interp_palette(['#333333', '#666666', '#999999'], len(single_feature_tax_pie_df))
        # Skip printing taxonomy level if dataframe is empty.
        if single_feature_tax_pie_df.empty:
            pass
        # Use regular wedge for first taxonomy level and add it to the pie chart.
        elif diamond_tax_level == "Superkingdom":
            tax_pie_chart.wedge(x=0, y=1, radius=0.2, start_angle=cumsum('angle', include_zero=True), end_angle=cumsum('angle'), line_color="white", color='color', alpha=0.7, source=single_feature_tax_pie_df)
        # Add legend for selected taxonomy level and add it using annular wedge.
        elif diamond_tax_level == diamond_tax_level_legend_dropdown.value:
            tax_pie_chart.annular_wedge(x=0, y=1, inner_radius=(diamond_level_num+1)/10, outer_radius=(diamond_level_num+2)/10, start_angle=cumsum('angle', include_zero=True), end_angle=cumsum('angle'), line_color="white", fill_color='color', alpha=0.7, legend_field='taxonomy', source=single_feature_tax_pie_df)
        # Add taxonomy level using annular wedge.
        else:
            tax_pie_chart.annular_wedge(x=0, y=1, inner_radius=(diamond_level_num+1)/10, outer_radius=(diamond_level_num+2)/10, start_angle=cumsum('angle', include_zero=True), end_angle=cumsum('angle'), line_color="white", fill_color='color', alpha=0.7, source=single_feature_tax_pie_df)
    # Hide the grid from the pie chart.
    tax_pie_chart.axis.axis_label = None
    tax_pie_chart.axis.visible = False
    tax_pie_chart.grid.grid_line_color = None
    # Clean up legend if we selected a taxonomy level.
    if diamond_tax_level_legend_dropdown.value:
        tax_pie_chart.legend.location = "top_left"
        tax_pie_chart.legend.label_text_font_size = "8px"
        tax_pie_chart.legend.label_text_font_style = "bold"
        tax_pie_chart.legend.border_line_alpha = 0
        tax_pie_chart.legend.background_fill_alpha = 0
    mo.vstack([tax_pie_chart])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""## Cog <a id="cog"></a>""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ### Cog counts per Cog Family <a id="cog_per_cog_family"></a>
    #### Prints a bar chart displaying how many Cog entries are present in each the Cog family.
    """
    )
    return


@app.cell(hide_code=True)
def _(
    get_cog_cog_family_dataframe,
    get_cog_family_dataframe,
    get_feature_cog_dataframe,
    plt,
):
    # Finds how the Cog Families are distributed in our dataset.
    # Merges 'feature has cog' and 'cog has cog family' tables.
    cog_data = get_feature_cog_dataframe().merge(get_cog_cog_family_dataframe(), on="cog_id").merge(get_cog_family_dataframe(), on="cog_family_id").groupby('cog_family_name')
    cog_data["cog_family_name"].agg(Count="count").sort_values(["Count"], ascending=False).plot.barh()
    plt.xlabel("Amount of Cogs")
    plt.ylabel("Cog Families")
    plt.show()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ### Cog counts for Cog Family per Bin <a id="cog_for_cog_per_bin"></a>
    #### Prints a bar chart displaying how many Cog entries from the selected Cog family are present in each bin. Takes the Cog family with the most amount of Cogs by default.
    """
    )
    return


@app.cell(hide_code=True)
def _(
    get_bin_dataframe,
    get_cog_cog_family_dataframe,
    get_cog_family_dataframe,
    get_feature_cog_dataframe,
    mo,
):
    # Takes the most common Cog Family as value.
    cog_family_counts_df = get_feature_cog_dataframe().merge(get_cog_cog_family_dataframe(), on="cog_id").merge(get_cog_family_dataframe(), on="cog_family_id").groupby(['cog_family_name'])["cog_family_name"].agg(Count="count").sort_values(["Count"], ascending=False).reset_index()
    # Dropdown with all Cog Family types.
    cog_family_name_dropdown = mo.ui.dropdown(options=cog_family_counts_df["cog_family_name"], label="Choose Cog Family:", value=cog_family_counts_df["cog_family_name"].iloc[0])
    # Pick the range of bins to display for large amount of datasets only.
    # Requires dataset with more than 30 bins by default.
    max_value_cbp = 40
    bin_count_cbp = int(len(get_bin_dataframe())/max_value_cbp)*max_value_cbp
    if bin_count_cbp < max_value_cbp:
        max_value_cbp = len(get_bin_dataframe())
    range_slider_bin_cog_barplot = mo.ui.slider(start=0, stop=bin_count_cbp, value=0, step=max_value_cbp, label="Display bin range:")
    mo.vstack([mo.hstack([cog_family_name_dropdown]), mo.hstack([range_slider_bin_cog_barplot])])
    return (
        cog_family_name_dropdown,
        max_value_cbp,
        range_slider_bin_cog_barplot,
    )


@app.cell(hide_code=True)
def _(
    cog_family_name_dropdown,
    get_cog_cog_family_dataframe,
    get_cog_family_dataframe,
    get_feature_cog_dataframe,
    get_feature_dataframe,
    max_value_cbp,
    plt,
    range_slider_bin_cog_barplot,
):
    # Looks how many Cogs from the selected Cog Family appear in each of the bins.
    chosen_cog_family = get_cog_family_dataframe().query(f"cog_family_name == '{cog_family_name_dropdown.value}'")["cog_family_id"].iloc[0]
    cog_ids = get_cog_cog_family_dataframe().query("cog_family_id == @chosen_cog_family")["cog_id"]
    get_feature_cog_dataframe().query("cog_id in @cog_ids").merge(get_feature_dataframe(), on="feature_id").groupby('bin_id')["bin_id"].agg(Count="count").iloc[range_slider_bin_cog_barplot.value:range_slider_bin_cog_barplot.value+max_value_cbp].plot.barh()
    plt.xlabel(f"Amount of Cogs from Cog Family {chosen_cog_family}")
    plt.ylabel("Bin Names")
    plt.title(f"{cog_family_name_dropdown.value} - - - Distribution")
    plt.show()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ### Cog Family counts per Bin heatmap <a id="cog_family_bin_heatmap"></a>
    #### Prints a heatmap displaying how many Cog entries are present between the Cog families and the bins.
    """
    )
    return


@app.cell(hide_code=True)
def _(get_bin_dataframe, mo):
    # If we should sort the heatmap or not based on its count instead of its name.
    sort_cogfam_count_heatmap = mo.ui.checkbox(label="Sort by amount:", value=True)
    # If we should normalize the data or not.
    normalize_cogfam_count_heatmap = mo.ui.checkbox(label="Normalize data:", value=True)
    # Pick the range of bins to display for large amount of datasets only.
    # Requires dataset with more than 30 bins by default.
    max_value_chm = 30
    bin_count_chm = int(len(get_bin_dataframe())/max_value_chm)*max_value_chm
    if bin_count_chm < max_value_chm:
        max_value_chm = len(get_bin_dataframe())
    range_slider_bin_cog_heatmap = mo.ui.slider(start=0, stop=bin_count_chm, value=0, step=max_value_chm, label="Display bin range:")
    mo.vstack([mo.hstack([sort_cogfam_count_heatmap]), mo.hstack([normalize_cogfam_count_heatmap]), mo.hstack([range_slider_bin_cog_heatmap])])
    return (
        max_value_chm,
        normalize_cogfam_count_heatmap,
        range_slider_bin_cog_heatmap,
        sort_cogfam_count_heatmap,
    )


@app.cell(hide_code=True)
def _(
    get_cog_cog_family_dataframe,
    get_cog_family_dataframe,
    get_feature_cog_dataframe,
    get_feature_dataframe,
    max_value_chm,
    normalize_cogfam_count_heatmap,
    np,
    plt,
    range_slider_bin_cog_heatmap,
    sns,
    sort_cogfam_count_heatmap,
):
    # Creates a heatmap using the amount of cogs for each cog family per bin.
    # Sorts the heatmap by count is True and otherwise sorts it by name.
    cog_counts = get_feature_cog_dataframe().merge(get_cog_cog_family_dataframe(), on="cog_id").merge(get_cog_family_dataframe(), on="cog_family_id").merge(get_feature_dataframe(), on="feature_id").groupby(["bin_id", "cog_family_name"])["cog_family_name"].agg(count="count").reset_index().pivot_table(values='count', index='bin_id', columns='cog_family_name').iloc[range_slider_bin_cog_heatmap.value:range_slider_bin_cog_heatmap.value+max_value_chm].T.replace(np.nan, 0)
    sns.set(rc={"figure.figsize":(10, 10)})
    # Normalization here is ** 0.5 which is the same as 2log().
    normalize_cogfam_amount = 1
    if normalize_cogfam_count_heatmap.value:
        normalize_cogfam_amount = 0.5
    # Create heatmap sorted by counts.
    if sort_cogfam_count_heatmap.value:
        cog_counts["sum"] = cog_counts.sum(axis=1, numeric_only=True)
        sns.heatmap(cog_counts.sort_values(["sum"], ascending=False).drop(columns="sum") ** normalize_cogfam_amount, cmap="inferno", square=True, vmin=0).set(xlabel="Bin Names", ylabel="Cog Family Names")
    # Create heatmap sorted by name.
    else:
        sns.heatmap(cog_counts ** normalize_cogfam_amount, cmap="inferno", square=True, vmin=0).set(xlabel="Bin Names", ylabel="Cog Family Names")
    plt.show()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""## Pfam <a id="pfam"></a>""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ### Pfam Clans with the highest Pfam count <a id="highest_pfam_count_per_clan"></a>
    #### Prints a bar chart displaying how many Pfam entries are present in each the Pfam clan.
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    # Here you can use the slider to pick how many of the most common Pfam Clans you want to display.
    pfam_clan_amount = mo.ui.slider(start=4, stop=20, label="Print amount:", value=20)
    mo.hstack([pfam_clan_amount])
    return (pfam_clan_amount,)


@app.cell(hide_code=True)
def _(
    get_feature_pfam_dataframe,
    get_pfam_clan_dataframe,
    get_pfam_dataframe,
    np,
    pfam_clan_amount,
    plt,
):
    # Finds the most common Pfam Clans present in our dataset.
    # Merges 'feature has pfam' and 'pfam' tables while removing all pfam ids that don't have a pfam clan.
    pfam_data = get_feature_pfam_dataframe().merge(get_pfam_dataframe(), on="pfam_id").replace('', np.nan).dropna().merge(get_pfam_clan_dataframe(), on="pfam_clan_id").groupby('pfam_clan_name')
    pfam_data["pfam_clan_name"].agg(Count="count").sort_values(["Count"], ascending=False).head(pfam_clan_amount.value).plot.barh()
    plt.ylabel("Pfam Clans")
    plt.xlabel("Amount of Pfams")
    plt.show()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ### Pfam Clan counts per Bin heatmap <a id="pfam_clan_bin_heatmap"></a>
    #### Prints a heatmap displaying how many Pfam entries are present between the Pfam clans and the bins.
    """
    )
    return


@app.cell(hide_code=True)
def _(get_bin_dataframe, mo):
    # Here you can use the slider to pick how many of the most common Pfam clans you want to display.
    pfam_clan_amount_heatmap = mo.ui.slider(start=4, stop=30, label="Print clan amount:", value=30)
    # If we should normalize the data or not.
    normalize_pfam_clan_count_heatmap = mo.ui.checkbox(label="Normalize data:", value=True)
    # Pick the range of bins to display for large amount of datasets only.
    # Requires dataset with more than 30 bins by default.
    max_value_phm = 30
    bin_count_phm = int(len(get_bin_dataframe())/max_value_phm)*max_value_phm
    if bin_count_phm < max_value_phm:
        max_value_phm = len(get_bin_dataframe())
    range_slider_bin_pfam_heatmap = mo.ui.slider(start=0, stop=bin_count_phm, value=0, step=max_value_phm, label="Display bin range:")
    mo.vstack([mo.hstack([pfam_clan_amount_heatmap]), mo.hstack([normalize_pfam_clan_count_heatmap]), mo.hstack([range_slider_bin_pfam_heatmap])])
    return (
        max_value_phm,
        normalize_pfam_clan_count_heatmap,
        pfam_clan_amount_heatmap,
        range_slider_bin_pfam_heatmap,
    )


@app.cell(hide_code=True)
def _(
    get_feature_dataframe,
    get_feature_pfam_dataframe,
    get_pfam_clan_dataframe,
    get_pfam_dataframe,
    max_value_phm,
    normalize_pfam_clan_count_heatmap,
    np,
    pfam_clan_amount_heatmap,
    plt,
    range_slider_bin_pfam_heatmap,
    sns,
):
    # Creates a heatmap using the amount of pfams for the pfam clan with the most found pfams per bin.
    pfam_counts = get_feature_pfam_dataframe().merge(get_pfam_dataframe(), on="pfam_id").replace('', np.nan).dropna().merge(get_pfam_clan_dataframe(), on="pfam_clan_id").merge(get_feature_dataframe(), on="feature_id").groupby(["bin_id", "pfam_clan_name"])["pfam_clan_name"].agg(count="count").reset_index().pivot_table(values='count', index='bin_id', columns='pfam_clan_name').iloc[range_slider_bin_pfam_heatmap.value:range_slider_bin_pfam_heatmap.value+max_value_phm].T.replace(np.nan, 0)
    pfam_counts["sum"] = pfam_counts.sum(axis=1, numeric_only=True)
    sns.set(rc={"figure.figsize":(10, 10)})
    # Normalization here is ** 0.5 which is the same as 2log().
    normalize_pfam_clan_amount = 1
    if normalize_pfam_clan_count_heatmap.value:
        normalize_pfam_clan_amount = 0.5
    # Create heatmap.
    sns.heatmap(pfam_counts.sort_values(["sum"], ascending=False).head(pfam_clan_amount_heatmap.value).drop(columns="sum") ** normalize_pfam_clan_amount, cmap="inferno", square=True, vmin=0).set(xlabel="Bin Names", ylabel="Pfam Clan Names")
    plt.show()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""## KO <a id="ko"></a>""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ### KO Jaccard Similarity Index for Bins <a id="ko_jaccard_bins"></a>
    #### Prints a heatmap displaying how similar the molecular functionality is between the bins.
    """
    )
    return


@app.cell(hide_code=True)
def _(get_feature_dataframe, get_feature_ko_dataframe, mo):
    # Select which bins you want to print:
    all_bin_ko_option = get_feature_ko_dataframe().merge(get_feature_dataframe(), on="feature_id")["bin_id"].unique()
    bin_ko_x_multiselect = mo.ui.multiselect(options=all_bin_ko_option, label="Choose bin IDs for x axis:", value=all_bin_ko_option[:5])
    bin_ko_y_multiselect = mo.ui.multiselect(options=all_bin_ko_option, label="Choose bin IDs for y axis:", value=all_bin_ko_option[:5])
    mo.vstack([mo.hstack([bin_ko_x_multiselect]), mo.hstack([bin_ko_y_multiselect])])
    return bin_ko_x_multiselect, bin_ko_y_multiselect


@app.cell(hide_code=True)
def _(
    bin_ko_x_multiselect,
    bin_ko_y_multiselect,
    get_feature_dataframe,
    get_feature_ko_dataframe,
    np,
    pd,
    plt,
    sns,
):
    # Warning, this calculation will take a long time depending on the amount of samples present in your dataset.

    # Uses Jaccard Similarity Index to check how similar the bins are from eachother using found KO IDs.
    # The code checks to see if the KEGG orthologies count is the same for each instance for each individual bin.
    # One means complete similarity and zero means no similarity (when not counting normalization).
    ko_data = get_feature_ko_dataframe().merge(get_feature_dataframe(), on="feature_id").groupby(["bin_id", "ko_id"])["ko_id"].agg(count="count").reset_index()
    ko_data_bin = ko_data.pivot_table(values='count', index='ko_id', columns='bin_id')

    jaccards_values = []
    jaccards_values_bin = []
    for ko_data_bin_x in bin_ko_y_multiselect.value:
        for ko_data_bin_y in bin_ko_x_multiselect.value:
            ko_counts = pd.concat([ko_data_bin[ko_data_bin_x], ko_data_bin[ko_data_bin_y]], axis=1).T.max()-abs(ko_data_bin.replace(np.nan, 0)[ko_data_bin_x]-ko_data_bin.replace(np.nan, 0)[ko_data_bin_y])
            jaccards_values_bin.append(ko_counts.dropna().sum()/pd.concat([ko_data_bin[ko_data_bin_x], ko_data_bin[ko_data_bin_y]], axis=1).T.max().dropna().sum())
        jaccards_values.append(jaccards_values_bin)
        jaccards_values_bin = []
    jaccards_index = pd.DataFrame(jaccards_values, columns=bin_ko_x_multiselect.value, index=bin_ko_y_multiselect.value)
    sns.set(rc={"figure.figsize":(10, 10)})
    # Normalization here is ** 0.5 which is the same as 2log().
    sns.heatmap(jaccards_index ** 0.5, cmap="inferno", square=True, vmin=0, vmax=1).set(xlabel="Bins", ylabel="Bins", title="Jaccard similarity index KO")
    plt.show()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ### Pathway Modules completeness for Bins <a id="pathway_module_cov_bins"></a>
    #### Prints a heatmap displaying how complete each module is in a pathway for all bins.
    """
    )
    return


@app.cell(hide_code=True)
def _(get_bin_dataframe, get_pathway_dataframe, mo):
    # Here you can pick your pathway ID.
    pathway_module_dropdown = mo.ui.dropdown(options=get_pathway_dataframe()["pathway_id"], label="Choose pathway ID:", value="map00010")
    # Input the pathway ID instead of choosing it with the dropdown.
    input_pathway_module = mo.ui.text(placeholder="00000 (Numbers Only)", label="Type Pathway Instead:")
    max_value_khm = 30
    bin_count_khm = int(len(get_bin_dataframe())/max_value_khm)*max_value_khm
    if bin_count_khm < max_value_khm:
        max_value_khm = len(get_bin_dataframe())
    range_slider_bin_kegg_heatmap = mo.ui.slider(start=0, stop=bin_count_khm, value=0, step=max_value_khm, label="Display bin range:")
    mo.vstack([mo.hstack([pathway_module_dropdown]), mo.hstack([input_pathway_module]), mo.hstack([range_slider_bin_kegg_heatmap])])
    return (
        input_pathway_module,
        max_value_khm,
        pathway_module_dropdown,
        range_slider_bin_kegg_heatmap,
    )


@app.cell(hide_code=True)
def _(
    get_feature_dataframe,
    get_feature_ko_dataframe,
    get_ko_dataframe,
    get_module_dataframe,
    get_pathway_dataframe,
    get_pathway_module_dataframe,
    input_pathway_module,
    max_value_khm,
    pathway_module_dropdown,
    pd,
    plt,
    range_slider_bin_kegg_heatmap,
    sns,
):
    # Check if the text box has an input.
    if input_pathway_module.value:
        pathway_var = "map"+ input_pathway_module.value
    else:
        pathway_var = pathway_module_dropdown.value
    # Fetches the map name and modules within it.
    pathway_modules = list(get_pathway_dataframe().query(f"pathway_id == '{pathway_var}'").merge(get_pathway_module_dataframe(), on='pathway_id')['module_id'])
    pathway_name = '-'.join(get_pathway_dataframe().query(f"pathway_id == '{pathway_var}'").values[0])

    # Module class to help calculate module completeness for each bin.
    class Module:
        def __init__(self, module, ko_list_ids, entry_len):
            self.module = module
            self.ko_list_ids = ko_list_ids
            self.entry_len = entry_len
            self.entry_check_list = [x + 1 for x in range(entry_len)]
            self.bin_check_list = {}

        def check_bin_ko_list(self, bin_id, ko_list):
            self.bin_check_list[bin_id] = self.entry_check_list.copy()
            for ko in ko_list:
                if ko in self.ko_list_ids.keys():
                    if self.ko_list_ids[ko] in self.bin_check_list[bin_id]:
                        self.bin_check_list[bin_id].remove(self.ko_list_ids[ko])

        def get_module_completion(self, bin_id):
            return (1 - len(self.bin_check_list[bin_id]) / self.entry_len) * 100

    parentheses_layer = 1
    construct_KO = 'K'
    modules_list = []
    # Constructs a list of Module classes that will store for each contig how complete each module is. This is not entirely accurate but gives an idea to how complete the module is.
    for module, module_id in zip(get_module_dataframe()['module_order'], get_module_dataframe()['module_id']):
        entry_layer = 1
        module_ko_dict = {}
        for char in module:
            if char == '(':
                parentheses_layer = parentheses_layer + 1
            elif char == ')':
                parentheses_layer = parentheses_layer - 1
            elif char == ' ' and parentheses_layer == 1:
                entry_layer = entry_layer + 1
            elif char == ' ':
                entry_layer = entry_layer + 1
            elif char == ',':
                pass
            elif char == 'K':
                construct_KO = 'K'
            else:
                construct_KO = construct_KO + char
                if len(construct_KO) == 6:
                    module_ko_dict[construct_KO] = entry_layer
                    construct_KO = 'K'
        modules_list.append(Module(module_id, module_ko_dict, entry_layer))
    percentage_module_list = []
    for ko_list in get_feature_dataframe().merge(get_feature_ko_dataframe(), on='feature_id').merge(get_ko_dataframe(), on='ko_id').groupby(['bin_id'])['ko_id']:
        for module_num in modules_list:
            module_num.check_bin_ko_list(ko_list[0], ko_list[1])
            percentage_module_list.append([ko_list[0], module_num.module, module_num.get_module_completion(ko_list[0])])
    percentage_module_df = pd.DataFrame(data=percentage_module_list, columns=['bin_id', 'module_id', 'percentage']).merge(get_module_dataframe(), on='module_id').pivot_table(values='percentage', index='bin_id', columns=['module_id', 'module_desc']).iloc[range_slider_bin_kegg_heatmap.value:range_slider_bin_kegg_heatmap.value+max_value_khm].T
    if not percentage_module_df.query('module_id == @pathway_modules').empty:
        sns.set(rc={'figure.figsize': (10, 10)})
        sns.heatmap(percentage_module_df.query('module_id == @pathway_modules'), cmap='inferno', square=True, vmin=0, vmax=100).set(xlabel='Bins', ylabel='Modules', title=pathway_name)
        plt.show()
    else:
        print('Warning, the dataframe is empty!\nThis means no modules are connected to this pathway!')
    #percentage_module_df.query('module_id == @pathway_modules')
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        """
    ### Pathway map KO coverage for Bins <a id="pathway_ko_map_bins"></a>
    #### Prints the selected KEGG pathway map and highlights in red _the KO/list of KOs_ that are present in the genes of the selected bin ID. For bigger maps [iPath3](https://pathways.embl.de/) might be more useful.
    """
    )
    return


@app.cell(hide_code=True)
def _(get_bin_dataframe, get_pathway_dataframe, mo):
    # Here you can pick your bin ID.
    bin_id_pathway_dropdown = mo.ui.dropdown(options=get_bin_dataframe()["bin_id"], label="Choose bin ID:", value=get_bin_dataframe()["bin_id"].iloc[0])
    # Here you can pick your pathway ID.
    bin_pathway_module_dropdown = mo.ui.dropdown(options="ko" + get_pathway_dataframe()["pathway_id"].str[3:], label="Choose pathway ID:", value="ko00010")
    # Input the pathway ID instead of choosing it with the dropdown.
    input_pathway_module_bin = mo.ui.text(placeholder="00000 (Numbers Only)", label="Type Pathway Instead:")
    # If we should print the pathway next to the original pathway.
    should_print_original_pathway = mo.ui.checkbox(label="Print original pathway next to ortholog pathway:", value=False)
    # If we should print the list with KOs of the bin.
    should_print_list_of_kos = mo.ui.checkbox(label="Print list of KOs from the bin:", value=False)
    mo.vstack([mo.hstack([bin_id_pathway_dropdown]), mo.hstack([bin_pathway_module_dropdown]), mo.hstack([input_pathway_module_bin]), mo.hstack([should_print_original_pathway]), mo.hstack([should_print_list_of_kos])])
    return (
        bin_id_pathway_dropdown,
        bin_pathway_module_dropdown,
        input_pathway_module_bin,
        should_print_list_of_kos,
        should_print_original_pathway,
    )


@app.cell(hide_code=True)
def _(
    bin_pathway_module_dropdown,
    get_feature_ko_dataframe,
    get_pathway_dataframe,
    input_pathway_module_bin,
    mo,
    should_print_list_of_kos,
):
    # Check if the text box has an input.
    if input_pathway_module_bin.value:
        pathway_module_var = "ko"+ input_pathway_module_bin.value
    else:
        pathway_module_var = bin_pathway_module_dropdown.value
    if should_print_list_of_kos.value:
        print('\n'.join(get_feature_ko_dataframe()["ko_id"].to_list()))
    # Extract pathway description and fetch the proper link to the pathway.
    pathway_desc_name = get_pathway_dataframe().query(f"pathway_id == '{pathway_module_var.replace("ko", "map")}'")["pathway_desc"].values[0]
    mo.md(f"Link to pathway: __[{pathway_desc_name}](https://www.kegg.jp/entry/{pathway_module_var.replace("ko", "map")})__.")
    return (pathway_module_var,)


@app.cell(hide_code=True)
def _(
    BytesIO,
    KGMLCanvas,
    KGML_parser,
    PIL_Image,
    REST,
    bin_id_pathway_dropdown,
    convert_from_path,
    get_feature_dataframe,
    get_feature_ko_dataframe,
    get_ko_dataframe,
    get_ko_ec_dataframe,
    mo,
    os,
    pathway_module_var,
    should_print_original_pathway,
):
    # Get all KOs for the chosen bin.
    ko_bin_check_list = get_feature_dataframe().query(f"bin_id == '{bin_id_pathway_dropdown.value}'").merge(get_feature_ko_dataframe(), on='feature_id').merge(get_ko_dataframe(), on='ko_id')["ko_id"].values

    # Use KGML to fetch pathway data.
    kgml_pathway = KGML_parser.read(REST.kegg_get(pathway_module_var, "kgml"))
    ko_ec_table = get_ko_ec_dataframe()
    # Change the colour of the KOs in our pathway list.
    for kgml_ko in kgml_pathway.orthologs:
        ko_all_names = kgml_ko.name.replace("ko:", "").split(" ")
        for kgml_single_entry in ko_all_names:
            if kgml_single_entry in ko_bin_check_list:
                for kmgl_graphic in kgml_ko.graphics:
                    kmgl_graphic.bgcolor = "#FFBFBF"
                    kmgl_graphic.fgcolor = "#FF0000"
    # Put it in canvas to safe it to a PDF (couldn't print it to PNG immediately).
    kgml_canvas = KGMLCanvas(kgml_pathway, import_imagemap=True)
    kgml_canvas.draw("temp_pathway_map.pdf")
    # Convert first PDF page into PNG file.
    convert_from_path("temp_pathway_map.pdf")[0].save("temp_pathway_map.png")
    # Remove our PDF file.
    os.remove("temp_pathway_map.pdf")
    # Print the PNG file.
    temp_pathway_image = PIL_Image.open("temp_pathway_map.png")

    # If we should print the original pathway next to the other pathway.
    if should_print_original_pathway.value:
    # Open map in PIL using bytes.
        original_pathway_image = PIL_Image.open(BytesIO(REST.kegg_get(pathway_module_var.replace("ko", "map"), "image").read())).convert("RGBA")
        original_pathway_image = original_pathway_image.resize((temp_pathway_image.width, temp_pathway_image.height))
        combined_pathway_image = PIL_Image.new('RGBA', (temp_pathway_image.width*2, temp_pathway_image.height))
        combined_pathway_image.paste(temp_pathway_image, (0, 0))
        combined_pathway_image.paste(original_pathway_image, (temp_pathway_image.width, 0))
        combined_pathway_image.save("temp_pathway_map.png")
    # Print the PNG file.
    temp_pathway_image = PIL_Image.open("temp_pathway_map.png")
    mo.hstack([temp_pathway_image])
    # Remove our PNG file (cannot print if removed).
    #os.remove("temp_pathway_map.png")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        """
    ### Benchmark Run Information <a id=benchmark_table></a>
    #### Prints a table displaying the different benchmark features for each snakemake rule.
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    # If we merge benchmarks are for the same rules but have different bins.
    should_merge_same_rules = mo.ui.checkbox(label="Merge identical rules:", value=True)
    mo.vstack([mo.hstack([should_merge_same_rules])])
    return (should_merge_same_rules,)


@app.cell(hide_code=True)
def _(StringIO, get_run_dataframe, os, pd, should_merge_same_rules):
    # Get the benchmark map location.
    benchmark_table_output_map = get_run_dataframe()["od"].values[0] + "/benchmark"
    benchmark_table_dir_list = os.listdir(benchmark_table_output_map)

    # Make a Pandas table with all the benchmarks.
    all_benchmarks_string_table = "rule\ts\th:m:s\tmax_rss\tmax_vms\tmax_uss\tmax_pss\tio_in\tio_out\tmean_load\tcpu_time"
    for benchmark_table_file in benchmark_table_dir_list:
        with open(f"{benchmark_table_output_map}/{benchmark_table_file}") as input_handle_benchmark_table:
            name_tb = benchmark_table_file.replace('.benchmark', '')
            info_tb = input_handle_benchmark_table.readlines()[-1].strip("\n").replace('NA', '0')
            all_benchmarks_string_table += f"\n{name_tb}\t{info_tb}"
    benchmark_table_df = pd.read_csv(StringIO(all_benchmarks_string_table), sep ="\t")
    if should_merge_same_rules.value:
        # Grab the first part of the rule name.
        benchmark_table_df["rule"] = benchmark_table_df["rule"].str.split('.').str[0]
        # Calculate time times the mean load.
        benchmark_table_df["s_mean_load"] = benchmark_table_df["s"] * benchmark_table_df["mean_load"]
        # Merge benchmark entries with the same rules. 
        benchmark_table_df = benchmark_table_df.groupby(["rule"]).aggregate({"s": "sum", "max_rss": "max", "max_vms": "max", "max_uss": "max", "max_pss": "max", "io_in": "sum", "io_out": "sum", "s_mean_load": "sum", "cpu_time": "sum"}).reset_index()
        # Input correct Hours:Minutes:Seconds.
        benchmark_table_df["h:m:s"] = (benchmark_table_df["s"] / 3600).astype(int).astype(str) + pd.to_datetime(benchmark_table_df["s"], unit='s').dt.strftime(':%M:%S')
        # Recalculate correct mean_load (total load / total seconds).
        benchmark_table_df["mean_load"] = benchmark_table_df["s_mean_load"] / benchmark_table_df["s"]
    # Print table.
    benchmark_table_df[["rule", "s", "h:m:s", "max_rss", "max_vms", "max_uss", "max_pss", "io_in", "io_out", "mean_load", "cpu_time"]].sort_values("s", ascending=False)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        """
    ### Benchmark Pie Chart <a id=benchmark_pie_chart></a>
    #### Prints a pie chart displaying the selected benchmark features for each snakemake rule.
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    # All benchmark display options.
    benchmark_options = {"Time": ["s", "sum", "@rule: @hms"], 
                         "Maximum Resident Set Size": ["max_rss", "max", "@rule: @max_rss MB"], 
                         "Maximum Virtual Memory Size": ["max_vms", "max", "@rule: @max_vms MB"], 
                         "Unique Set Size": ["max_uss", "max", "@rule: @max_uss MB"], 
                         "Proportional Set Size": ["max_pss", "max", "@rule: @max_pss MB"], 
                         "MB Read": ["io_in", "sum", "@rule: @io_in MB"], 
                         "MB Written": ["io_out", "sum", "@rule: @io_out MB"], 
                         "CPU Usage Mean": ["mean_load", "sum", "@rule: @mean_load"], 
                         "CPU Time (Divided by 100)": ["cpu_time", "sum", "@rule: @cpu_hms"]}
    # What type of benchmark to display in the pie chart.
    benchmark_type_dropdown = mo.ui.dropdown(options=benchmark_options.keys(), label="Choose benchmark type:", value="Time")
    # How many unique benchmark entries to print in the pie chart.
    display_benchmarks_amount = mo.ui.number(start=3, stop=20, label="Amount of unique benchmark entries to display:", value=7)
    # Display UI.
    mo.vstack([mo.hstack([benchmark_type_dropdown]), mo.hstack([display_benchmarks_amount])])
    return (
        benchmark_options,
        benchmark_type_dropdown,
        display_benchmarks_amount,
    )


@app.cell(hide_code=True)
def _(
    Spectral11,
    StringIO,
    benchmark_options,
    benchmark_type_dropdown,
    cumsum,
    display_benchmarks_amount,
    figure,
    get_run_dataframe,
    interp_palette,
    math,
    mo,
    os,
    pd,
):
    # Get the benchmark map location.
    benchmark_output_map = get_run_dataframe()["od"].values[0] + "/benchmark"
    benchmark_dir_list = os.listdir(benchmark_output_map)

    # Make a Pandas table with all the benchmarks.
    all_benchmarks_string = "rule\ts\th:m:s\tmax_rss\tmax_vms\tmax_uss\tmax_pss\tio_in\tio_out\tmean_load\tcpu_time"
    for benchmark_file in benchmark_dir_list:
        with open(f"{benchmark_output_map}/{benchmark_file}") as input_handle_benchmark:
            name = benchmark_file.replace('.benchmark', '')
            info = input_handle_benchmark.readlines()[-1].strip("\n").replace('NA', '0')
            all_benchmarks_string += f"\n{name}\t{info}"
    benchmark_df = pd.read_csv(StringIO(all_benchmarks_string), sep ="\t")

    # Grab the rules with the highest run time and display the rest in the other rules row.
    benchmark_df["rule"] = benchmark_df["rule"].str.split('.').str[0]
    benchmark_df_merged = benchmark_df.groupby(["rule"]).aggregate({benchmark_options[benchmark_type_dropdown.value][0]: benchmark_options[benchmark_type_dropdown.value][1]})
    # Only display the benchmarks with the most time.
    top_benchmarks_listed = benchmark_df_merged.sort_values([benchmark_options[benchmark_type_dropdown.value][0]], ascending=False).index[:display_benchmarks_amount.value]
    # Display the rest of the benchmarks in the other rules category.
    benchmark_df.loc[~benchmark_df["rule"].isin(top_benchmarks_listed), "rule"] = "other rules"
    # Calculate time times the mean load.
    benchmark_df["s_mean_load"] = benchmark_df["s"] * benchmark_df["mean_load"]
    # Merge benchmark entries with the same rules. 
    benchmark_df = benchmark_df.groupby(["rule"]).aggregate({"s": "sum", "max_rss": "max", "max_vms": "max", "max_uss": "max", "max_pss": "max", "io_in": "sum", "io_out": "sum", "s_mean_load": "sum", "cpu_time": "sum"}).reset_index()
    # Input correct Hours:Minutes:Seconds.
    benchmark_df["hms"] = (benchmark_df["s"] / 3600).astype(int).astype(str) + pd.to_datetime(benchmark_df["s"], unit='s').dt.strftime(':%M:%S')
    # Input correct Hours:Minutes:Seconds.
    benchmark_df["cpu_hms"] = (benchmark_df["cpu_time"] / 360000).astype(int).astype(str) + pd.to_datetime(benchmark_df["cpu_time"] / 100, unit='s').dt.strftime(':%M:%S')
    # Recalculate correct mean_load (total load / total seconds).
    benchmark_df["mean_load"] = benchmark_df["s_mean_load"] / benchmark_df["s"]
    benchmark_df["angle"] = benchmark_df[benchmark_options[benchmark_type_dropdown.value][0]]/benchmark_df[benchmark_options[benchmark_type_dropdown.value][0]].sum() * 2 * math.pi
    benchmark_df["color"] = interp_palette(Spectral11, len(benchmark_df))
    if benchmark_type_dropdown.value == "Time":
        benchmark_df["legend"] = benchmark_df["rule"] + "\n" + benchmark_df["hms"]
    elif benchmark_type_dropdown.value == "CPU Time (Divided by 100)":
        benchmark_df["legend"] = benchmark_df["rule"] + "\n" + benchmark_df["cpu_hms"]
    elif benchmark_options[benchmark_type_dropdown.value][2][-2:] == "MB":
        benchmark_df["legend"] = benchmark_df["rule"] + "\n" + (round(benchmark_df[benchmark_options[benchmark_type_dropdown.value][0]] / 1024, 2)).astype(str) + " GB"
    else:
        benchmark_df["legend"] = benchmark_df["rule"] + "\n" + round(benchmark_df[benchmark_options[benchmark_type_dropdown.value][0]], 2).astype(str)

    # Make a pie chart from the benchmakrs.
    benchmark_pie_chart = figure(width=600, height=400, title=f"Benchmark Pie Chart - - - {benchmark_type_dropdown.value}", tooltips=benchmark_options[benchmark_type_dropdown.value][2], x_range=(-2, 1))
    benchmark_pie_chart.wedge(x=0, y=1, radius=0.8, start_angle=cumsum('angle', include_zero=True), end_angle=cumsum('angle'), line_color="white", color='color', alpha=0.7, legend_field='legend', source=benchmark_df)
    # Hide the grid from the pie chart.
    benchmark_pie_chart.axis.axis_label = None
    benchmark_pie_chart.axis.visible = False
    benchmark_pie_chart.grid.grid_line_color = None
    # Clean up legend if we selected a taxonomy level.
    benchmark_pie_chart.legend.location = "top_left"
    benchmark_pie_chart.legend.label_text_font_size = "8px"
    benchmark_pie_chart.legend.label_text_font_style = "bold"
    benchmark_pie_chart.legend.border_line_alpha = 0
    benchmark_pie_chart.legend.background_fill_alpha = 0
    mo.vstack([benchmark_pie_chart])
    return


@app.cell
def _():
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
