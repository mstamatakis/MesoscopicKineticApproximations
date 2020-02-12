# Performs comparison of two files and reports if they are not the same (w.r.p. to pre-defined tolerance)  
from os.path import join, exists

def loadtext(path):
  with open(path, 'r') as file:
    file.readline() 
    return [ [float(u) for u in line.split()] for line in file ] 

def run(input_file_name, directory_results, directory_reference):
  res = 0
  num_tests = len(directory_reference)

  input_file_name_full = open(join(directory_reference, input_file_name), 'r')
  approx_name = input_file_name_full.readline()   # read the approximation name

  # print('Reading file for: ', approx_name)
  str_tmp = ''
  file_name = str_tmp.join((approx_name, '_Fortran_Theta_vs_Mu.txt'))

  calculated_data = loadtext(join(directory_results, file_name))
  reference_data  = loadtext(join(directory_reference, file_name))

  if (len(calculated_data) != len(reference_data)):
    res = 1
  else:
    res = 0
  if (all(abs(x-y) < 5e-6 for l0, l1 in zip(calculated_data, reference_data)
                            for x, y in zip(l0, l1) )):
    res = 0
  else:
    res = 1

  return res, approx_name

if __name__ == '__main__':
  directory_results = '001'
  input_file_name = 'calculation_input.dat'
  directory_reference = ['../tests/data/0/', '../tests/data/1/', '../tests/data/2/']

  print('\n')

  for n in range(0, len(directory_reference)):
    res, approx_name = run(input_file_name, directory_results, directory_reference[n])
    if (res == 0):
      print approx_name, ':', 'Succeeded!'
    else:
      print approx_name, ':', 'Failed!'
  
  print('\n')