from flask import Flask, request, jsonify, render_template
from py.seq_search import run as seq_search
app = Flask(__name__)

# ph_seq = [
#     {
#         "record": "P0DPB8",
#         "itens": [
#             {
#                 "id": "cas13",
#                 "seq": "MWISIKTLIHHLGVLFFCDYMYNRREKKIIEVKTMRITKVEVDRKKVLISRDKNGGKLVY",
#                 "species": "Listeria seeligeri"
#             }
#         ]
#     }
# ]


# @app.get("/seq") # Endpoint: http://127.0.0.1:5000/seq
# def get_seq():
#     return {"seq": ph_seq}
#

@app.route('/', methods=["POST", "GET"])
def gfg():
    if request.method == "POST":
        # getting input with name = fname in HTML form
        sequence = request.form["search-box"]
        # return "Your sequence is " + sequence
        return "Your sequence " + sequence
    return render_template("index.html")

#
# @app.post("/seq")
# def create_seq():
#     request_data = request.get_json()
#     # Creates an internal dictionary that is structured following the /seq format
#     new_seq = {"record": request_data["record"], "itens": []}
#     ph_seq.append(new_seq)
#     # Accepts the submission and OK with the 201 status code
#     return new_seq, 201
#
# #
# @app.post("/seq/<string:id>/item")
# def create_item(id):
#     request_data = request.get_json()
#     for s in ph_seq:
#         if s["record"] == id:
#             new_item = {"id": request_data["id"], "seq": request_data['seq']}
#             s["itens"].append(new_item)
#             return new_item, 201
#     return {"message": "Entry not found"}, 404

# @app.route('/')
# def home():
#     return render_template('index.html')
#
#
# @app.route('/form/submit', methods=['GET'])
# def process_form_data():
#     selected_mode = request.args.get('mode')
#     if selected_mode == 'mode_1':
#         # process user input for mode 1
#         pass
#     elif selected_mode == 'mode_2':
#         # process user input for mode 2
#         pass
#     elif selected_mode == 'mode_3':
#         # process user input for mode 3
#         pass
#     else:
#         # handle error for invalid input
#         pass
#
#     # return
#
#


if __name__=='__main__':
    app.run()
