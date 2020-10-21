'''
Created on 21 oct. 2020

@author: alba
'''

if __name__ == '__main__':
    pass

# file_name = input("Input file name:")
# name_parts = file_name.split(".")
# print(name_parts)

#extract extension function

def extension (file_name):
    import os.path
    name, extension = os.path.splitext(file_name)
    print ("The file extension is {}". format(extension))
extension("example_annot.gbf")