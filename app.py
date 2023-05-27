# Native modules
import os
# External modules
from flask import Flask, request, render_template
# Project modules
import py.seq_search as seq_search
from py.render_search_result import run as render_blastout_table
from py.render_word_search import run as render_word_search
import yaml
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


@app.route('/', methods=["POST", "GET"])
def gfg():
    if request.method == "POST":
        # HTML handling object - Handles temp files and page rendering
        dynamic_html_toolbox = DynamicHtmlTemplate(tbl_render_config)

        # GET INPUT from the search box
        user_raw_input = request.form["search-box"]
        # FASTA temp file handling
        temp_fasta_path = f"{sq_search_config['temp_fasta_dir']}{os.sep}"
        user_input_path = f"{temp_fasta_path}{dynamic_html_toolbox.dynamic_file_prefix}.fasta"
        # Create temp file to hold the user input
        with open(user_input_path, 'w') as temp_user_input:
            temp_user_input.write(user_raw_input)

        # INPUT ASSESSMENT: Test sequence format and proceed accordingly
        blastout_report_dict = seq_search.run(user_input_path,
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
            return render_template(html_template_path)

        # SCENARIO 2: User input IS NOT a sequence:
        #   Deliver word search output result
        if not blastout_report_dict:
            print("Format word search")
            # Create word search output page
            html_word_search_tbl = render_word_search(
                str(user_raw_input.strip("\s+").rstrip("\s+")),
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
