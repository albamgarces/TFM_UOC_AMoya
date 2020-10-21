'''
Created on 21 oct. 2020

@author: alba
'''

if __name__ == '__main__':
    pass

def extension (file_name):
    import os.path
    compt = {".gbf", ".gff", ".fna", ".faa"}
    name, extension = os.path.splitext(file_name)
    if extension in compt:
        print ("The file extension is {}". format(extension))
    else:
        print("This file is not compatible")
extension("example_annot.gbf")
extension("example_annot.gbp")
