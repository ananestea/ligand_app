from flask import Flask, render_template, request, send_file
from ligand.bsa_calc import Bsa
import uuid
import os
import zipfile
import shutil

# from flask_cors import CORS

model = Bsa()

app = Flask(__name__)
# CORS(app, resources={r"/*": {"origins": "*"}})

@app.route('/', methods=['GET', 'POST'])
def index():
    return render_template('index.html')

def save_zip(arch, folder_list, mode):
    # Счетчики
    num = 0
    num_ignore = 0
    # Создание нового архива
    z = zipfile.ZipFile(arch, mode, zipfile.ZIP_DEFLATED, True)
    # Получаем папки из списка папок.
    for add_folder in folder_list:
        # Список всех файлов и папок в директории add_folder
        for root, dirs, files in os.walk(add_folder):
            for file in files:
                # Создание относительных путей и запись файлов в архив
                path = os.path.join(root, file)
                z.write(path)
                print(num, path)
                num += 1
    z.close()
    print("------------------------------")
    print("Добавлено: ", num)
    print("Проигнорировано: ", num_ignore)


def calc_ligand_func(model_num, distance, dir_calc_model, id):
    # Запустим расчет 
    model.main(model_num, distance, dir_calc_model)
    # Сохраним результаты в zip файле в отдельной папке
    save_zip(f"./ligand_zip/ligand_{id}.zip", [f'./ligand/current/{id}'], "w")
    # Удалим папку
    shutil.rmtree(f'./ligand/current/{id}')



@app.route('/calc-ligand', methods=['POST'])
def calc_ligand():
    # Сгенерим папку с уникальным id для расчета
    id = uuid.uuid4()
    dir_calc_model = f'./ligand/current/{id}'
    try:
        path = os.path.join(dir_calc_model)
        os.mkdir(path)
    except:
        print("бля вернкть 500 ошибка")
    # Переименнуем pdb файл и сохраним его в уникальную папку
    model_num = int(request.form['model-num'])
    distance = float(request.form['distance'])
    pdb_file = request.files['file']
    pdb_file.save(f'{dir_calc_model}/Cyclodextrine_ligand_names_good.pdb')
    # Запустим расчет
    calc_ligand_func(model_num, distance, dir_calc_model, id)
    #Зарендерить объект или страницу с информацией что типа ваш расчет готовится
    return f"Обязательно запомните ваш уникальный id {id}"


@app.route('/find-your-ligand', methods=['POST'])
def find_your_ligand():
    id = request.form['ligand-id']
    try:
        # Выполним поиск папки
       file = f'./ligand_zip/ligand_{id}.zip'
       return send_file(file, as_attachment=True)
    except: 
        return "Не найден такой id"
    
    
if __name__ == "__main__":
    app.run(debug=True)