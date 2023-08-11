# Native modules
import os
import re
import shutil
import pickle
from pathlib import Path
# External modules
from flask import Flask, request, render_template, send_from_directory
import yaml
from urllib.parse import quote as encode
from urllib.parse import unquote as decode
# Project modules
import py.seq_search as seq_search
from py.render_search_result import run as render_blastout_table
from py.render_word_search import run as render_word_search
from py.tool_finder_search import run as tool_finder
from py.render_wiki import run as load_wiki_data
from py.render_tool_finder import run as load_tool_finder
# from py.render_word_search import sql_table_to_df
# from py.db_loadNupdate import psql_connect, get_absolute_path

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
        self.wiki_template_path = f'{render_config["wiki_template_path"]}'
        # Create search output page internal  path/name
        self.internal_path = f'{render_config["temp_html_dir_path"]}{os.sep}{self.dynamic_file_prefix}.html'
        self.wiki_dir_path = f'{render_config["temp_wiki_dir_path"]}{os.sep}'

    def export_html_template(self, table_content):
        with open(f"{self.root_path}{self.internal_path}", 'w') as html_template:
            html_template.write(table_content)
        # Apply Flask method to render the search output page
        return self.internal_path


# Register a custom filter to treat None as empty string
@app.template_filter('empty_string')
def empty_string(value):
    if value is None:
        return ''
    return value


@app.route('/index.html')
def index():
    """Define route to the front page"""
    return render_template('index.html')


@app.route('/coming_soon.html')
def coming():
    """Define route to the coming soon page"""
    return render_template('coming_soon.html')


@app.route('/phylogeny_viewer.html')
def phylo_viewer():
    """Define route to the philogeny viewer page"""
    return render_template('phylogeny_viewer.html')


@app.route('/fusion_proteins.html')
def fusion_prot():
    """Define route to the Fusion proteins page"""
    return render_template('fusion_proteins.html')


@app.route('/faq.html')
def faq_page():
    """Define route to the FAQ page"""
    return render_template('faq.html')


@app.route('/contact_us.html')
def contactus():
    return render_template('/contact_us.html')


@app.route('/example.html')
def example():
    return render_template('/SpyCas9a_backup.html')


@app.route('/tool_finder.html')
def buffet():
    """Define route to tool_finder.html"""

    # Load tool finder options from PSQL Database
    tool_finder_menu = load_tool_finder(psql_config)
    # Render the Wiki page and apply the empty_string filter to avoid ugly 'None's in the page
    return render_template(f'tool_finder.html',
                           targe_type=tool_finder_menu.target_type or empty_string,
                           trans_activity=tool_finder_menu.trans_activity or empty_string,
                           targeting_requirement=tool_finder_menu.targeting_requirement or empty_string,
                           multiplex=tool_finder_menu.multiplex or empty_string
                           )


@app.route('/tool_finder', methods=['POST'])
def process_choices():
    dynamic_html_toolbox = DynamicHtmlTemplate(tbl_render_config)
    data = request.form
    tool_features_input = {}
    for key, value in data.items():
        tool_features_input[str(key)] = str(value)

    print("Format Tool Finder search")
    # Create word search output page
    html_word_search_tbl = tool_finder(tool_features_input, tbl_render_config, psql_config)
    # Export and render word search output page
    html_template_path = dynamic_html_toolbox.export_html_template(html_word_search_tbl)
    # Render page
    return render_template(html_template_path)
    # Process the choices as needed


@app.route('/wiki/<page>')
def wiki_page(page):
    # Infer the entry identifier from the HTML page route
    entry_id = re.sub(fr"\b.html\b", '', page)

    # # Get the page name associated with a given Wiki entry and set a filepath
    new_wiki_path = f'{tbl_render_config["dir_root_template_path"]}{os.sep}' \
                    f'{tbl_render_config["temp_wiki_dir_path"]}{os.sep}{page}'
    #
    # Check whether that wiki entry exists or not. And create a new page if necessary
    # wiki_path_obj = Path(new_wiki_path)
    shutil.copy(tbl_render_config["wiki_template_path"], new_wiki_path)

    # Get the wiki pickles path associated with a given Wiki entry and set a filepath
    pickles_path = f'{psql_config["pickles_path"]}{os.sep}{entry_id}.pkl'
    pickle_path_obj = Path(pickles_path)

    # Check whether that wiki pickle entry exists or not. And create a new page if necessary
    # The pickles are created from the information stored in the PSQL Database
    # This is leveraged to create a DynamicWiki object
    #   Within this object all the information is formatted to be incorporated into the wiki HTML template
    if not pickle_path_obj.is_file():
        wiki_entry = load_wiki_data(page, psql_config)
    if pickle_path_obj.is_file():
        with open(pickles_path, "rb") as wiki_pickle:
            wiki_entry = pickle.load(wiki_pickle)
    # In case there's no content for a given entry render an error page
    if not wiki_entry.content_check:
        return render_template('wiki/error.html')

    # Render the Wiki page and apply the empty_string filter to avoid ugly 'None's in the page
    # encoded_page = encode.quote(page)
    return render_template(f'wiki/{page}',
                           page_name=wiki_entry.entry_id or empty_string,
                           classification=wiki_entry.classification or empty_string,
                           classification_sprites=wiki_entry.classification_sprites or empty_string,
                           properties=wiki_entry.properties or empty_string,
                           resources=wiki_entry.resources or empty_string,
                           sequences=wiki_entry.sequences or empty_string,
                           text_summaries=wiki_entry.text_summaries or empty_string,
                           # gene_editing_human=wiki_entry.gene_editing_human or empty_string,
                           gene_editing=wiki_entry.gene_editing or empty_string,
                           tools=wiki_entry.tools or empty_string,
                           variants=wiki_entry.variants or empty_string,
                           exp_details=wiki_entry.exp_details or empty_string,
                           sequence_browser=wiki_entry.sequence_browser or empty_string,
                           pfam=wiki_entry.pfam or empty_string,
                           domains=wiki_entry.domains or empty_string,
                           structure=wiki_entry.structure or empty_string,
                           references=wiki_entry.formatted_references or empty_string
                           )


@app.route('/static/fasta/<filename>')
def render_file(filename):
    # Get the absolute path of the file in the /fasta/ folder
    file_path = os.path.join('static/fasta', filename)

    # Check if the file exists
    if os.path.exists(file_path):
        # Send the file as download
        # return send_file(file_path, as_attachment=True)
        # Pop file on the browser
        return send_from_directory('static/fasta', filename, mimetype='text/plain')
    else:
        # File not found, return an error response
        return f"File not found: {filename}", 404


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
