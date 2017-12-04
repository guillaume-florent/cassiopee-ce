# Genere la liste des fichiers HDF5 a compiler
# dans le fichier files
import os

file = open('files', 'w')
listFiles = os.listdir('.')
for i in listFiles:
    if (i.find('.c') != -1):
        file.write('           \'Converter/HDF/'+i+'\',\n')
file.close()

