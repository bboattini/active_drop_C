import os
import numpy as np

# Find the parent directory of the current path
PATH = str(os.path.abspath(__file__))
PATH = PATH.replace('/'+ PATH.split("/")[-2] +'/' + PATH.split("/")[-1], "")
items = os.listdir(PATH)
COLORS = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w', 'black', '#33FFF5', '#FF8C33', '#8C33FF', '#33FF8C', '#FF3333', '#33A1FF']
MARKERS = ['o','d', 's', 'p', '*', 'h', '+', 'x']

def a_from_file (File):
  fs = File.split("/")[-1].split("_")
  fl = len(fs)
  pad = int(fs[8].replace("", ""))
  return pad # return a value

def h_from_file (File):
  fs = File.split("/")[-1].split("_")
  fl = len(fs)
  pad = int(fs[10].replace("", ""))
  return pad # return h value

def state_from_file (File):
  fs = File.split("/")[-1].split("_")
  fl = len(fs)
  pad = fs[0].replace("/", "")
  return pad # return state value

def fo_from_file (File):
  fs = File.split("/")[-2].split("_")
  fl = len(fs)
  pad = fs[-1]
  return pad # return state value

def file_crawler(path = PATH):
    items = os.listdir(path)

    configuration_files = []
    last_configuration_files = []
    measures_files = []
    general_data_files = []
    for item in items:
        full_item_path = os.path.join(path, item)
        if os.path.isdir(full_item_path):
            if 'dados_' in item:
                # List the files inside this directory
                files = os.listdir(full_item_path)
                #print(files)
                for File in files:
                    if 'conf.dsf' in File:
                        configuration_files.append(os.path.join(full_item_path, File))
                    elif 'LAST.dsf' in File:
                        last_configuration_files.append(os.path.join(full_item_path, File))
                    elif '.out' in File:
                        general_data_files.append(os.path.join(full_item_path, File))
                    elif '.dsf' in File:
                        measures_files.append(os.path.join(full_item_path, File))
                    

    return {'config': configuration_files ,'last': last_configuration_files, 'measures': measures_files, 'general': general_data_files, }

def time_series_to_frequency_array(x_data, y_data = [], bin = 1):
    if len(y_data) == 0:
        l = len(x_data)
        y_data = x_data[l//2:]
        x_data = x_data[:l//2]
    # Get the maximum and minimum values of x_data or y_data
    max_value = max(max(np.append(x_data,0)), max(np.append(y_data,0)))
    min_value = min(min(np.append(x_data,1000)), min(np.append(y_data,1000)))
    # Values range
    values_range = np.arange(int(min_value), int(max_value+1), 1/bin)
    # Create a dictionary to store the frequency of each value
    frequency_dict = {}
    # Iterate over the unique values
    for value in values_range:
        # Count the number of occurrences of the value in the x]_data and y_data array
        #frequency = np.sum(x_data.astype(int) == value) + np.sum(y_data.astype(int) == value)
        value = np.round(value * bin) / bin
        frequency = np.sum(np.round(x_data*bin)/bin == value) + np.sum(np.round(y_data*bin)/bin == value)

        # Store the frequency in the dictionary
        frequency_dict[value] = frequency
    # Create an array with the unique values and their frequencies
    frequency_array = np.array([[key, frequency_dict[key]] for key in frequency_dict])
    return frequency_array

