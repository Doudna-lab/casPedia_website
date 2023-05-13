import os

from flask import Flask, request, jsonify, render_template
from py.seq_search import run as seq_search
import random
import string
import datetime
app = Flask(__name__)


def random_name_gen():
    # Generate Unique filename to store blastouts
    n = 20
    prefix = ''.join(random.choices(string.ascii_letters +
                                    string.digits, k=n))
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
    return f"{prefix}_{timestamp}"


@app.route('/', methods=["POST", "GET"])
def gfg():
    if request.method == "POST":

        # getting input from search box
        sequence = request.form["search-box"]
        seq_path = f"jobs/seqIN_{random_name_gen()}.fasta"
        with open(seq_path, 'w') as f:
            f.write(sequence)
        #
        blast_result = seq_search(seq_path)
        os.remove(seq_path)
        return render_template(blast_result)
    else:
        return render_template("index.html")


if __name__ =='__main__':
    app.run()
