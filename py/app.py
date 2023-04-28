from flask import Flask, request, render_template

app = Flask(__name__)

@app.route('/')
def home():
    return render_template('index.html')


@app.route('/form/submit', methods=['GET'])
def process_form_data():
    selected_mode = request.args.get('mode')
    if selected_mode == 'mode_1':
        # process user input for mode 1
        pass
    elif selected_mode == 'mode_2':
        # process user input for mode 2
        pass
    elif selected_mode == 'mode_3':
        # process user input for mode 3
        pass
    else:
        # handle error for invalid input
        pass

    # return


if __name__ == '__main__':
    app.run(debug=True)
