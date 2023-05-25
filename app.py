# Native modules
import os
# External modules
from flask import Flask, request, render_template
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
    tbl_render_config = yaml.safe_load(f)
#   -> PostgreSQL Database config
with open("config/table_update.yaml", "r") as f:
    psql_config = yaml.safe_load(f)


@app.route('/', methods=["POST", "GET"])
def gfg():
    if request.method == "POST":
        # Generate a random file identifier for the job
        file_dynamic_prefix: str = seq_search.random_name_gen()
        # Page templates path
        root_template_dir_path = tbl_render_config['root_template_dir_path']

        # Get input from search box
        user_raw_input = request.form["search-box"]
        user_input_path = f"{sq_search_config['temp_fasta_dir']}{os.sep}{file_dynamic_prefix}.fasta"
        with open(user_input_path, 'w') as temp_fasta:
            temp_fasta.write(user_raw_input)

        # Test sequence format and perform sequence search if applicable
        blastout_report_dict = seq_search.run(user_input_path, sq_search_config, file_dynamic_prefix)

        # Deliver sequence search output result
        if blastout_report_dict:
            print("Format blastout")
            # Clean up temp sequence
            os.remove(user_input_path)
            # Create sequence search output page
            html_blastout_tbl = render_blastout_table(blastout_report_dict, tbl_render_config)
            # Create search output page internal  path/name
            html_tbl_template_path = f'{tbl_render_config["temp_html_dir_path"]}{os.sep}{file_dynamic_prefix}.html'

            # Export search output page  to temporary HTML file
            with open(f"{root_template_dir_path}{os.sep}{html_tbl_template_path}", 'w') as blastout_template:
                blastout_template.write(html_blastout_tbl
                                        )
            return render_template(html_tbl_template_path)

        # Deliver word search output result
        if not blastout_report_dict:
            print("Format word search")

            #TODO: Continuar com integracao do word search
            # Create word search output page
            html_word_search_tbl = render_word_search(user_input_path, psql_config)
    else:
        return render_template("index.html")


if __name__ == '__main__':
    app.run()
