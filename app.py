# Native modules
import os
import re
# External modules
from flask import Flask, request, render_template, jsonify
import yaml
# Project modules
import py.seq_search as seq_search
from py.render_search_result import run as render_blastout_table
from py.render_word_search import run as render_word_search
from py.render_wiki import run as render_wiki_data
from py.render_word_search import sql_table_to_df
from py.db_loadNupdate import psql_connect, get_absolute_path

# Flask setup
app = Flask(__name__)

# Import config files
#   -> Sequence search config
with open("config/seq_search.yaml", "r") as f:
    sq_search_config = yaml.load(f, Loader=yaml.FullLoader)
#   -> Table rendering config
with open("config/render_result.yaml", "r") as f:
    tbl_render_config: object = yaml.safe_load(f)
#   -> PostgreSQL DB interaction config
with open("config/db_interaction.yaml", "r") as f:
    psql_config = yaml.safe_load(f)


class DynamicHtmlTemplate:
    def __init__(self, render_config):
        # Generate a random file identifier for the job
        self.dynamic_file_prefix = seq_search.random_name_gen()
        # Page templates path
        self.root_path = f"{render_config['dir_root_template_path']}{os.sep}"
        # Create search output page internal  path/name
        self.internal_path = f'{render_config["temp_html_dir_path"]}{os.sep}{self.dynamic_file_prefix}.html'

    def export_html_template(self, table_content):
        with open(f"{self.root_path}{self.internal_path}", 'w') as html_template:
            html_template.write(table_content)
        # Apply Flask method to render the search output page
        return self.internal_path


@app.route('/buffet.html')
def buffet():
    """Define route to buffet.html"""
    return render_template('buffet.html')


@app.route('/cas_buffet', methods=['POST'])
def process_choices():
    data = request.form
    choices = {}
    for key, value in data.items():
        choices[key] = value
    # Process the choices as needed
    return f'Choices processed successfully {choices}'


@app.route('/index.html')
def index():
    """Define route to the front page"""
    return render_template('index.html')


@app.route('/wiki/<page>')
def wiki_page(page):
    """Define route to individual protein wiki pages"""
    wiki_entry = render_wiki_data(page, psql_config)
    return render_template(f'wiki/{page}')
    # TODO: Figure out how to best handle dynamic wiki HTML generation
    # return render_template(f'wiki/{page}',
    #                        properties=wiki_entry.properties,
    #                        resources=wiki_entry.resources,
    #                        text_summaries=wiki_entry.text_summaries,
    #                        gene_editing_human=wiki_entry.gene_editing_human,
    #                        gene_editing=wiki_entry.gene_editing,
    #                        tools=wiki_entry.tools,
    #                        variants=wiki_entry.variants,
    #                        exp_details=wiki_entry.variants
    #                        )


@app.route('/', methods=["POST", "GET"])
def gfg():
    if request.method == "POST":
        # HTML handling object - Handles temp files and page rendering
        dynamic_html_toolbox = DynamicHtmlTemplate(tbl_render_config)

        # GET INPUT from the search box
        user_raw_input = request.form["search-box"]
        user_clean_input = re.sub(r'\s+', '', user_raw_input).strip()
        # User input temp file handling
        temp_input_path = f"{sq_search_config['temp_fasta_dir']}{os.sep}"
        user_input_path = f"{temp_input_path}{dynamic_html_toolbox.dynamic_file_prefix}.in"
        # Create temp file to hold the user input
        with open(user_input_path, 'w') as temp_user_input:
            temp_user_input.write(user_raw_input)

        # INPUT ASSESSMENT: Test sequence format and proceed accordingly
        blastout_report_dict, format_check = seq_search.run(user_input_path,
                                                            sq_search_config,
                                                            dynamic_html_toolbox.dynamic_file_prefix)

        # SCENARIO 1: User input IS a sequence:
        #   Deliver sequence search output result
        if blastout_report_dict:
            print("Format blastout")
            # Clean up temp sequence
            os.remove(user_input_path)
            # Create sequence search output page
            html_blastout_tbl = render_blastout_table(
                blastout_report_dict,
                tbl_render_config
            )
            # Export blastP search output page
            html_template_path = dynamic_html_toolbox.export_html_template(html_blastout_tbl)
            # Render page
            print("Render blastout HTML result")
            return render_template(html_template_path)

        # SCENARIO 2: User input IS NOT a sequence:
        #   Deliver word search output result
        if not blastout_report_dict:
            print("No results found in BlastP")
            # Check if the FASTA format was validated
            #   Validated FASTAs with no report means: The sequence yielded no results
            if format_check:
                print("BlastP yielded no results")
                return render_template("search_out_page.html")

            print("Format word search")
            # Create word search output page
            html_word_search_tbl = render_word_search(
                str(user_clean_input),
                tbl_render_config, psql_config
            )
            # Export and render word search output page
            html_template_path = dynamic_html_toolbox.export_html_template(html_word_search_tbl)
            # Render page
            return render_template(html_template_path)

    else:
        return render_template("index.html")


if __name__ == '__main__':
    app.run()
