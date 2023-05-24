# Native modules
import os
# External modules
from flask import Flask, request, render_template
import py.seq_search as seq_search
from py.render_search_result import run as render_blastout_table
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


@app.route('/', methods=["POST", "GET"])
def gfg():
    if request.method == "POST":
        # Generate a random file identifier for the job
        file_dynamic_prefix: str = seq_search.random_name_gen()

        # Get input from search box
        fasta_sequence = request.form["search-box"]
        seq_path = f"{sq_search_config['temp_fasta_dir']}{os.sep}{file_dynamic_prefix}.fasta"
        with open(seq_path, 'w') as temp_fasta:
            temp_fasta.write(fasta_sequence)

        # Test sequence format and perform sequence search if applicable
        blast_result = seq_search.run(seq_path, sq_search_config, file_dynamic_prefix)

        # Deliver search output result
        if blast_result:
            print("Format blastout")
            # Clean up temp sequence
            os.remove(seq_path)
            # Create search output page
            html_blastout_tbl = render_blastout_table(blast_result, tbl_render_config)
            # Create search output page path/name
            html_tbl_template_path = f'{tbl_render_config["root_template_dir_path"]}{os.sep}{file_dynamic_prefix}.html'
            # Export search output page
            with open(html_tbl_template_path, 'w') as blastout_template:
                blastout_template.write(html_blastout_tbl)

            return render_template(f"{file_dynamic_prefix}.html")

    else:
        return render_template("index.html")

    # Clean up temp html page
    os.remove(f"{file_dynamic_prefix}.html")


if __name__ == '__main__':
    app.run()
