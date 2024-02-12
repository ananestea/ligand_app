import zipfile
import os

backup_folders = ['./ligand']  # Список папок для архивации
arch_name = "ligand_" + "egwetryewrgewu" + ".zip"  # имя архива!
ignore_file = []  # если надо исключить файлы

def mybackup(arch, folder_list, mode):
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
                if file in ignore_file:  # Исключаем лишние файлы
                    print("Исключен! ", str(file))
                    num_ignore += 1
                    continue
                # Создание относительных путей и запись файлов в архив
                path = os.path.join(root, file)
                z.write(path)
                print(num, path)
                num += 1
    z.close()
    print("------------------------------")
    print("Добавлено: ", num)
    print("Проигнорировано: ", num_ignore)


mybackup(arch_name, backup_folders, "w")