# -*- coding: utf-8 -*-
"""
Created on Tue Jan 12 00:02:36 2021

@author: khaik
"""
import os, json, requests, pickle
import numpy as np

def hex_to_RGB(hex):
  ''' "#FFFFFF" -> [255,255,255] '''
  return [int(hex[i:i+2], 16) for i in range(1,6,2)]

def RGB_to_hex(RGB):
  ''' [255,255,255] -> "#FFFFFF" '''
  RGB = [int(x) for x in RGB]
  return "#"+"".join(["0{0:x}".format(v) if v < 16 else
            "{0:x}".format(v) for v in RGB])

def color_dict(gradient):
  ''' Takes in a list of RGB sub-lists and returns dictionary of
    colors in RGB and hex form for use in a graphing function
    defined later on '''
  return {"hex":[RGB_to_hex(RGB) for RGB in gradient],
      "r":[RGB[0] for RGB in gradient],
      "g":[RGB[1] for RGB in gradient],
      "b":[RGB[2] for RGB in gradient]}


def linear_gradient(start_hex, finish_hex="#FFFFFF", n=10):
  ''' returns a gradient list of (n) colors between
    two hex colors. start_hex and finish_hex
    should be the full six-digit color string,
    inlcuding the number sign ("#FFFFFF") '''
  # Starting and ending colors in RGB form
  s = hex_to_RGB(start_hex)
  f = hex_to_RGB(finish_hex)
  # Initilize a list of the output colors with the starting color
  RGB_list = [s]
  # Calcuate a color at each evenly spaced value of t from 1 to n
  for t in range(1, n):
    # Interpolate RGB vector for color at the current value of t
    curr_vector = [
      int(s[j] + (float(t)/(n-1))*(f[j]-s[j]))
      for j in range(3)
    ]
    # Add it to our list of output colors
    RGB_list.append(curr_vector)

  return color_dict(RGB_list)

def getAPI(url):
    response = requests.get(url)
    if response.status_code != 200:
        raise Exception('oops')
    return response.json()

def jprint(data):
    print(json.dumps(data, sort_keys=True, indent=4))
    

def savecsv(save_path, filenames, files, ftype='array', 
            save_index=False, label=None):
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    
    if ftype=='array':
        [np.savetxt(os.path.join(save_path, filenames[i]), files[i], delimiter=",") for i in range(len(filenames))];

    elif ftype=='dataframe': 
        [files[i].to_csv(os.path.join(save_path, filenames[i]), index=save_index, index_label=label) for i in range(len(filenames))];


def savepickle(save_path, filenames, files):
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    
    for i in range(len(filenames)):
        filename=os.path.join(save_path,filenames[i])
        with open(filename, 'wb') as f:
            pickle.dump(files[i], f)