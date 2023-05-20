import pandas as pd
import os
from py.seq_search import random_name_gen


# Define a function to generate HTML links
def generate_link(row, config):
    link = f'<a href="{row[config["linked_column"]]}.html">{row[config["linked_column"]]}</a>'
    return link


def order_df_columns(df):
    # Get the column names and pop the last column
    columns = df.columns.tolist()
    last_column = columns.pop()
    # Reorder the column names and insert the last column at the beginning
    new_columns = [last_column] + columns
    # Reindex the DataFrame with the new column order
    df_out = df.reindex(columns=new_columns)

    return df_out


def dynamic_blastout_html(blastout_df_html, html_template_path):
    # Read the template file
    with open(html_template_path, 'r') as f:
        template = f.read()

    # Replace the placeholder with the table HTML
    template = template.replace('<div id="table"></div>', blastout_df_html)
    template = template.replace('<table border="1" class="dataframe">', '<table border="1" id="table">')
    template = template.replace('<thead>', '<thead id="thead">')
    template = template.replace('<tbody>', '<tbody id="tbody">')
    return template


# DEBUG INPUTS
# import yaml
# with open("config/render_result.yaml", "r") as f:
#     config = yaml.load(f, Loader=yaml.FullLoader)
def run(blastout_dict, config, random_file_prefix):
    # Import blastout parsed dictionary into pandas dataframe
    df = pd.DataFrame.from_dict(blastout_dict[list(blastout_dict.keys())[0]], orient='index')

    # Create a column for hit IDs
    df[config["linked_column"]] = df.index
    # Reorder the dataframe so the hit IDs are placed first
    df = order_df_columns(df)
    # Apply the function to links to the hit ID columns
    df[config["linked_column"]] = df.apply(lambda row: generate_link(row, config), axis=1)

    # Convert the DataFrame to HTML table
    df_blastout_html = df.to_html(escape=False, index=False)

    # Create and save a modified template to a new file
    blastout_html_template = dynamic_blastout_html(df_blastout_html, config["search_template_path"])

    return blastout_html_template
