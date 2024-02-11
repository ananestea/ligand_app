from flask import Flask, render_template, request
from ligand.bsa_calc import Bsa
import uuid
import os

# from flask_cors import CORS

model = Bsa()


app = Flask(__name__)
# CORS(app, resources={r"/*": {"origins": "*"}})

@app.route('/', methods=['GET', 'POST'])
def index():
    return render_template('index.html')


# @app.route('/calc-ligand', methods=['POST'])
@app.route('/calc-ligand')
def calc_ligand():
    # Сгенерим папку с уникальным id для расчета
    dir_calc_model = f'./ligand/{uuid.uuid4()}'
    try:
        path = os.path.join(dir_calc_model)
        os.mkdir(path)
    except:
        print("бля вернкть 500 ошибка")
    # Переименнуем pdb файл в ...  и сохраним его в уникальную папку
    
    # Запустим расчет 
    model.main(20, 1.75, dir_calc_model)

    # print(request)

    return

@app.route('/find-your-ligand', methods=['POST'])
def find_your_ligand():

    return
    
    
if __name__ == "__main__":
    app.run(debug=True)