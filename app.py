# Native modules
import os
# External modules
from flask import Flask, request, jsonify, render_template
from py.seq_search import run as seq_search
from py.seq_search import random_name_gen
app = Flask(__name__)


@app.route('/', methods=["POST", "GET"])
def gfg():
    if request.method == "POST":

        # Get input from search box
        fasta_sequence = request.form["search-box"]
        seq_path = f"jobs/seqIN_{random_name_gen()}.fasta"
        with open(seq_path, 'w') as f:
            f.write(fasta_sequence)

        # Test sequence format and perform sequence search if applicable
        blast_result = seq_search(seq_path)

        # Deliver search output result
        if blast_result:
            os.remove(seq_path)
            return render_template(blast_result)

    else:
        return render_template("index.html")


if __name__ == '__main__':
    app.run()
