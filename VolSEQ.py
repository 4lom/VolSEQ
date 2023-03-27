## MIT License
##
## Copyright (c) 2023 Zeissloff Louis
##
## Permission is hereby granted, free of charge, to any person obtaining a copy
## of this software and associated documentation files (the "Software"), to deal
## in the Software without restriction, including without limitation the rights
## to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
## copies of the Software, and to permit persons to whom the Software is
## furnished to do so, subject to the following conditions:
##
## The above copyright notice and this permission notice shall be included in all
## copies or substantial portions of the Software.
##
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
## IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
## FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
## AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
## LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
## OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
## SOFTWARE.

# Import Definition:
import base64
import PySimpleGUI as sg

from pathlib import Path
from tempfile import TemporaryDirectory
import warnings

import io
from PIL import Image, ImageTk

from pydicom import dcmread
from pydicom.fileset import FileSet
from pydicom.uid import generate_uid
from pydicom.pixel_data_handlers.util import apply_modality_lut
from pydicom.errors import InvalidDicomError
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg



import numpy as np

import csv

# CONST Definition:
_VARS = {'WINDOW': False}
_VARS = {'APP_NAME': "VolSEQ"}

xStartingColumnCSV = 19
yStartingColumnCSV = 20
deltaColumnCSV = 5
invertROIImageCount = 1

strTooltipDICOMDIR = "Choose DICOMDIR file. The file MUST be grouped with corresponding DICOM files (same folder)"
strTooltipROI = "Upload ROI file as a csv file"
strTooltipXColumnROI = "Input column number for 1st pxX value. . These have the values of the points of the ROI."
strTooltipYColumnROI = "Input column number for 1st pxY value. . These have the values of the points of the ROI."
strTooltipdeltaColumnROI = "Calculate number of columns separating each pxX/pxY value"
strTooltipINMIN = " "
strTooltipINMAX = " "
strTooltipINVERTVALUE = "Invert the order in which the dicom Image are processed. In case the ROI and the DICOM does not have the same order."
strTooltipProcess = "Process Data"
strTooltipROISUM = "Volume of segmented tissue of choice"
strTooltipROICLEANSUM = "Total ROI volume"
strTooltipSLICESSUM = "Sum of pixels for the current slice"
strTooltipSERIE = "If there are multiple CT scans associated to that DICOM file, these can be chosen. Note: if the series number is changed, the analysis has to be reprocessed."
strTooltipSLIDER = "Slide to choose a specific slice in the series (mostly for graphic purposes)"
strTooltipIMAGE = "Image of CT scan with ROI overlap (in red) and pixels of only segmented tissue (in white)"
strTooltipSAVEFILE = "The chosen image can be saved here"

# Funtion Definition:

def isStrEmpty(string):
    return ("".__eq__(string))

def draw_figure(canvas, figure):
    figure_canvas_agg = FigureCanvasTkAgg(figure, canvas)
    figure_canvas_agg.draw()
    figure_canvas_agg.get_tk_widget().pack(side='top', fill='both', expand=1)
    return figure_canvas_agg

def displayArray(arr):
    plt.imshow(arr)
    plt.show()

def displayPolygon(poly_path):
    fig, ax = plt.subplots()
    patch = patches.PathPatch(poly_path, facecolor='orange', lw=2)
    ax.add_patch(patch)
    ax.set_xlim(200, 400)
    ax.set_ylim(200, 400)
    plt.show()

def saveFigOfSliceWithROI(arr, poly_path):
    im = plt.imshow(arr, cmap='gray')
    patch = patches.PathPatch(poly_path, facecolor=None, edgecolor='red', fill=False, lw=1)
    ax = plt.gca()
    ax.add_patch(patch)

    #Save plot into a base64 encoding
    io_buf = io.BytesIO()
    plt.savefig(io_buf, format='png')
    img = Image.open(io_buf)
    img.save(io_buf, format="png")
    imageBuffer = base64.b64encode(io_buf.getvalue())
    io_buf.close()

    ax.patches.pop()

    return imageBuffer

#Retrieve all row of a CSV file excluding the header.
def processCSVFile(fileName):
    imageNumbs = []
    xsRows = []
    ysRows = []
    with open(fileName, 'r') as CSVFile:
        CSVReader = csv.reader(CSVFile, delimiter=',')
        header = next(CSVReader)
        for row in CSVReader:
            #Only save the Image Name and the pxX and pyY colons
            if(len(row) > 0):
                rowArray = np.array(row)
                imageNumbs.append(int(rowArray[0]))
                xsRow = np.array(rowArray[xStartingColumnCSV::deltaColumnCSV], float)
                ysRow = np.array(rowArray[yStartingColumnCSV::deltaColumnCSV], float)
                xsRows.append(xsRow)
                ysRows.append(ysRow)
    CSVFile.close()
    return imageNumbs, xsRows, ysRows

def createROIPolygon(currxs, currys):
    #Set all values under zero to zero
    currxs[currxs<0]
    currys[currys<0]

    polygon = np.array([(currxs), (currys)])
    polygon = polygon.T

    return Path(polygon)

def createHUMask(minHU, maxHU, HU):
    mask = (HU > minHU) & (HU < maxHU)
    return mask

def computeSliceArray(ds, array):
    # [row, column] [mm, mm]
    pixelSpacingArray = ds.PixelSpacing
    pixelSpacing = pixelSpacingArray[0]*pixelSpacingArray[1]
    # In mm
    sliceThickness = ds.SliceThickness
    sliceTotal = np.sum(array)*sliceThickness*np.prod(pixelSpacing)
    return sliceTotal

def applyROIFilter(HUMask, poly_path):
    HUMaskFiltered = np.array(HUMask)

    #Get each coordinate where useful HU are.
    filterHUIndexList = np.where(HUMaskFiltered == True)

    #Convert to numpy array of 2xY
    filterHUIndexList = np.array(list(zip(filterHUIndexList[0], filterHUIndexList[1])))

    #Get each useful HU in ROI polygon
    isInPolyList = poly_path.contains_points(filterHUIndexList)
    
    #Get the indexes of filterHUIndexList not in ROI polygon
    notROIIndexes = np.where(isInPolyList == False)

    #Set filterHUIndexList values not in ROI polygon to 0.  
    for notROIIndex in notROIIndexes[0]:
        index = filterHUIndexList[notROIIndex]
        HUMaskFiltered[index[0], index[1]] = 0

    return HUMaskFiltered

def applyROIFilterWithPicture(HUMask, poly_path):
    HUMaskFiltered = np.array(HUMask)

    #Get each coordinate where useful HU are.
    filterHUIndexList = np.where(HUMaskFiltered == True)

    #Convert to numpy array of 2xY
    filterHUIndexList = np.array(list(zip(filterHUIndexList[0], filterHUIndexList[1])))

    #Get each useful HU in ROI polygon
    isInPolyList = poly_path.contains_points(filterHUIndexList)
    
    #Get the indexes of filterHUIndexList not in ROI polygon
    notROIIndexes = np.where(isInPolyList == False)

    #Set filterHUIndexList values not in ROI polygon to 0.  
    for notROIIndex in notROIIndexes[0]:
        index = filterHUIndexList[notROIIndex]
        HUMaskFiltered[index[0], index[1]] = 0

    return HUMaskFiltered


def readDICOMDIRData(DICOMDIRpath):
    # fetch the path to the test data


    fs = FileSet(DICOMDIRpath)
    print(fs)
    print()

    # We can search the File-set
    patient_ids = fs.find_values("PatientID")
    
    #SOP Instance of first Patient ID
    result = fs.find(PatientID=patient_ids[0])

    #SOP Instance of first Patient studie
    study_uids = fs.find_values("StudyInstanceUID", instances=result)
    result = fs.find(PatientID=patient_ids[0], StudyInstanceUID=study_uids[0])

    # Search available series
    series_uids = fs.find_values("SeriesInstanceUID", instances=result)
    
    for series_uid in series_uids:
        result = fs.find(
                    PatientID=patient_ids,
                    StudyInstanceUID=study_uids,
                    SeriesInstanceUID=series_uids)

                
    return series_uids

#Engine
def processInputData(serie, DICOMDIRpath, ROIpath, minHU, maxHU):
    # A File-set can be loaded from the path to its DICOMDIR dataset or the
    #   dataset itself
    fs = FileSet(DICOMDIRpath)
    # A summary of the File-set's contents can be seen when printing
    print(fs)
    print()

    #Find all images of the serie
    images = fs.find(SeriesInstanceUID=serie)
    
    maxImageNumber = len(images)
    imageNumb, xs, ys = processCSVFile(ROIpath)

    fullArray = np.ones((500,500))
    figListOfFigOfSliceWithROI = []
    figListOfFilteredFigOfSliceWithROI = []
    ROITotal = 0
    ROICleanTotal = 0
    slicesTotal = 0
    # Iterating over the FileSet yields FileInstance objects
    for instance in images:
        # Load the corresponding SOP Instance dataset
        ds = instance.load()

        #ROI slice number start at the end of the dicom Instance number.
        if( invertROIImageCount < 1 ):
            imageIndex = ds.InstanceNumber
        else:
            imageIndex = maxImageNumber - ds.InstanceNumber +1
        
        #Skip image not in the ROI list
        if(imageNumb.count(imageIndex) > 0):
            imageIndex = imageNumb.index(imageIndex)

            #Apply a modality lookup table or rescale operation to arr.
            #eq to int32(image * info.RescaleSlope) + info.RescaleIntercept.
            HU = apply_modality_lut(ds.pixel_array, ds)

            #Get the HU array and create a boolean array mask to have all value between minHU and maxHU
            HUMask = createHUMask(minHU, maxHU, HU)

            #Create a polygon with the ROI coordinate
            poly_path = createROIPolygon(xs[imageIndex], ys[imageIndex])

            newFig = saveFigOfSliceWithROI(HU, poly_path)
            figListOfFigOfSliceWithROI.append(newFig)

            newFilteredFig = saveFigOfSliceWithROI(HUMask, poly_path)
            figListOfFilteredFigOfSliceWithROI.append(newFilteredFig)
            
            filteredHUMask = applyROIFilterWithPicture(HUMask, poly_path)
            
            
            #Compute HU between minHU and maxHU in the ROI area
            ROISliceTotal = computeSliceArray(ds, filteredHUMask)

            #Compute all pixel in the ROI area
            filteredfullArray = applyROIFilter(fullArray, poly_path)
            ROISliceTotalClean = computeSliceArray(ds, filteredfullArray)

            ROITotal = ROITotal + ROISliceTotal

            ROICleanTotal = ROICleanTotal + ROISliceTotalClean

            #Compute all HU between minHU and maxHU
            sliceTotal = computeSliceArray(ds, HUMask)
            slicesTotal = slicesTotal + sliceTotal
            
    return slicesTotal, ROICleanTotal, ROITotal, figListOfFigOfSliceWithROI, figListOfFilteredFigOfSliceWithROI

def verifyInputData(DICOMDIRpath, ROIpath, minHU, maxHU):
    if(isStrEmpty(DICOMDIRpath) and isStrEmpty(ROIpath) and isStrEmpty(minHU) and isStrEmpty(maxHU)):
        sg.Popup('Error:', 'Some input are missing. Please enter all input.')
        return False
    #elif(minHU > maxHU):
    #    sg.Popup('Error:', 'Minimum HU value is greater than maximum HU value.')
    #    return False
    try:
        f1 = open(DICOMDIRpath, "r")
        f2 = open(ROIpath, "r")
        float(values['-INMIN-'])
        float(values['-INMAX-'])
        f1.close()
        f2.close()
    except IOError:
        sg.Popup('Error:', 'File does not exist.')
        return False
    except ValueError:
        sg.Popup('Error:', 'Input HU values are not correct numbers.')
        return False
    return True

def saveImageInPNGFile(imageData, filePath):
    #Decod image
    #img = Image.open(io.BytesIO(base64.b64decode(imageData)))
    
    images = [Image.open(io.BytesIO(base64.b64decode(x))) for x in imageData ]
    widths, heights = zip(*(i.size for i in images))
    total_width = sum(widths)
    max_height = max(heights)

    img = Image.new('RGB', (total_width, max_height))

    x_offset = 0
    for im in images:
      img.paste(im, (x_offset,0))
      x_offset += im.size[0]
    
    img.save(filePath) 

def displaySliceWithROI(imageData):
    #Display image
    imageBuffer = ImageTk.PhotoImage(image=Image.open(io.BytesIO(base64.b64decode(imageData))))
    _VARS['WINDOW']['-IMAGE-'].update(data=imageBuffer)

# Main Definition:
sg.theme('DarkTeal2')   # Theme
#Window layout

tissueComboChoice = [ 'Bone',
                      'Liver',
                      'White mater',
                      'Grey mater',
                      'Blood',
                      'Muscle',
                      'Kidney',
                      'CSP',
                      'Water',
                      'Fat',
                      'Air',
                      'Custom']



layout = [  [sg.T("")], [sg.Text("Choose a DICOMDIR file: "), sg.Input('', key='-DICOMDIR-', enable_events=True, tooltip = strTooltipDICOMDIR), sg.FileBrowse()],
            [sg.T("")], [sg.Text("Choose a ROI file: "), sg.Input('', key='-ROI-', enable_events=True, tooltip = strTooltipROI), sg.FileBrowse()],
            [sg.Text('Starting column for pxX in ROI file :'), sg.Input(xStartingColumnCSV, key='-XColumnROI-', enable_events=True, tooltip = strTooltipXColumnROI)],
            [sg.Text('Starting column for pyY in ROI file :'), sg.Input(yStartingColumnCSV, key='-YColumnROI-', enable_events=True, tooltip = strTooltipYColumnROI)],
            [sg.Text('Number of columns between pxX/pxY values in ROI file :'), sg.Input(deltaColumnCSV, key='-deltaColumnROI-', enable_events=True, tooltip = strTooltipdeltaColumnROI)],
            [sg.Text('Enter HU minimum '), sg.Input('0', key='-INMIN-', enable_events=True, tooltip = strTooltipINMIN)],
            [sg.Text('Enter HU maximum '), sg.Input('0', key='-INMAX-', enable_events=True, tooltip = strTooltipINMAX)],
            #[sg.Text('Decrement DICOM image number'), sg.Slider(range=(0,1), key='-INVERTVALUE-', enable_events = True, orientation='h', size=(34, 20), default_value=1, tooltip = strTooltipINVERTVALUE)],
            [sg.Checkbox('Decrement DICOM image number', key='-INVERTVALUE-', default=True, tooltip = strTooltipINVERTVALUE)],
            [sg.Button('Process', tooltip = strTooltipProcess)],
            [sg.Text('Segmented tissue value:'), sg.Text('0.0', key='-ROISUM-', tooltip = strTooltipROISUM),sg.Text('mm^3', tooltip = strTooltipROISUM)],
            [sg.Text('ROI total volume:'), sg.Text('0.0', key='-ROICLEANSUM-', tooltip = strTooltipROICLEANSUM),sg.Text('mm^3', tooltip = strTooltipROICLEANSUM)],
            #[sg.Text('Slices value:'), sg.Text('0.0', key='-SLICESSUM-', tooltip = strTooltipSLICESSUM),sg.Text('mm^3', tooltip = strTooltipSLICESSUM)],
            [sg.Text('Series:'), sg.Spin([i for i in range(1,1)], initial_value=1, key='-SERIE-', tooltip = strTooltipSERIE)],
            [sg.Text("Choose a Slice with ROI area:"), sg.Slider(range=(1, 1), key='-SLIDER-', disabled = True, enable_events = True, orientation='h', size=(34, 20), default_value=1, tooltip = strTooltipSLIDER)],
            [sg.Image(key='-IMAGE-', tooltip = strTooltipIMAGE)],
            [sg.In(key='-SAVEFILEPATH-', enable_events = True ,visible = False), sg.SaveAs(key='-SAVEFILE-', file_types=((('PNG', '*.png'),)), default_extension='.png', disabled = True, tooltip = strTooltipSAVEFILE)]]

# Create the Window
_VARS['WINDOW'] = sg.Window(_VARS['APP_NAME'],
                            layout,
                            finalize=True,
                            resizable=True,
                            element_justification="left")

figListOfFigOfSliceWithROI = []
figListOfFilteredFigOfSliceWithROI = []
series = []
# Event Loop to process "events" and get the "values" of the inputs
while True:
    event, values = _VARS['WINDOW'].read()
    if event == sg.WIN_CLOSED or event == 'Cancel': # if user closes window or clicks cancel
        break
    if event == '-DICOMDIR-':
        try:
            series = readDICOMDIRData(values['-DICOMDIR-'])
            valid_value = [i for i in range(1,len(series)+1)]
            _VARS['WINDOW']['-SERIE-'].update(values=valid_value, value=1)
        except FileNotFoundError:
            sg.Popup('Error:', 'DICOMDIR file not found.')
        except IndexError:
                sg.Popup('Error:', 'DICOMDIR file can not be read.')
    if event == '-XColumnROI-' and values['-XColumnROI-'] and values['-XColumnROI-'][-1] not in ('-0123456789.'):
        _VARS['WINDOW']['-XColumnROI-'].update(values['-XColumnROI-'][:-1])
    if event == '-YColumnROI-' and values['-YColumnROI-'] and values['-YColumnROI-'][-1] not in ('-0123456789.'):
        _VARS['WINDOW']['-YColumnROI-'].update(values['-YColumnROI-'][:-1])
    if event == '-deltaColumnROI-' and values['-deltaColumnROI-'] and values['-deltaColumnROI-'][-1] not in ('-0123456789.'):
        _VARS['WINDOW']['-deltaColumnROI-'].update(values['-deltaColumnROI-'][:-1])
    if event == '-INMIN-' and values['-INMIN-'] and values['-INMIN-'][-1] not in ('-0123456789.'):
        _VARS['WINDOW']['-INMIN-'].update(values['-INMIN-'][:-1])
    if event == '-INMAX-' and values['-INMAX-'] and values['-INMAX-'][-1] not in ('-0123456789.'):
        _VARS['WINDOW']['-INMAX-'].update(values['-INMAX-'][:-1])

    if event == 'Process':
        inputOk = verifyInputData(values['-DICOMDIR-'], values['-ROI-'], values['-INMIN-'], values['-INMAX-'])
        if(inputOk):
            try:
                xStartingColumnCSV = int(values['-XColumnROI-']) -1
                yStartingColumnCSV = int(values['-YColumnROI-']) -1
                deltaColumnCSV = int(values['-deltaColumnROI-'])
                _VARS['WINDOW']['-SLIDER-'].update(disabled = True)
                _VARS['WINDOW']['-SAVEFILE-'].update(disabled = True)
                slicesSum, ROICleanSum, ROISum, figListOfFigOfSliceWithROI, figListOfFilteredFigOfSliceWithROI = processInputData(series[int(values['-SERIE-'])-1], values['-DICOMDIR-'], values['-ROI-'], float(values['-INMIN-']), float(values['-INMAX-']))
                _VARS['WINDOW']['-ROISUM-'].update(ROISum)
                _VARS['WINDOW']['-ROICLEANSUM-'].update(ROICleanSum)
                #_VARS['WINDOW']['-SLICESSUM-'].update(slicesSum)
                _VARS['WINDOW']['-SAVEFILE-'].update(disabled = False)
                _VARS['WINDOW']['-SLIDER-'].update(value = 1, range=(1, len(figListOfFigOfSliceWithROI)), disabled = False)
                if(len(figListOfFigOfSliceWithROI)> 0):
                    displaySliceWithROI(figListOfFilteredFigOfSliceWithROI[0])
            except InvalidDicomError:
                sg.Popup('Error:', 'DICOMDIR input is not a DICOMDIR file.')
            except IndexError:
                sg.Popup('Error:', 'DICOMDIR file can not be read.')
            except ValueError:
                sg.Popup('Error:', 'ROI .csv file data can not be used.')
    if event == '-INVERTVALUE-':
        invertROIImageCount = int(values['-INVERTVALUE-'])
    if event == '-SLIDER-':
        sliderValue = int(values['-SLIDER-'])
        if(len(figListOfFigOfSliceWithROI)> 0):
            displaySliceWithROI(figListOfFilteredFigOfSliceWithROI[sliderValue-1])
    if event == '-SAVEFILEPATH-':
        file_path = values['-SAVEFILEPATH-']
        if file_path:
            sliderValue = int(values['-SLIDER-'])
            saveImageInPNGFile( [figListOfFigOfSliceWithROI[sliderValue-1], figListOfFilteredFigOfSliceWithROI[sliderValue-1]], file_path)
        
_VARS['WINDOW'].close()
