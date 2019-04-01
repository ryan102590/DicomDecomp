import pydicom
import os
import numpy
import numpy as np
import natsort
from natsort import natsorted
from os import system, name
try:
    import Tkinter as tk # this is for python2
except:
    import tkinter as tk # this is for python3
import glob2
import tkinter.filedialog
import pydicom.uid 
from pydicom.uid import ExplicitVRBigEndian
from pydicom.uid import ExplicitVRLittleEndian 
from pydicom.uid import ImplicitVRLittleEndian
from pydicom.uid import DeflatedExplicitVRLittleEndian 
import math
import logging
logger = logging.getLogger(__name__)

print('called all protocols')
def get_INPUT_directory():
    root = tk.Tk()
    root.withdraw()
    INPUT_DIR = tkinter.filedialog.askdirectory(title = 'Please Select CT Dicom directory Reference')
    return INPUT_DIR
	
def get_INPUT2_directory():
    root = tk.Tk()
    root.withdraw()
    INPUT_DIR = tkinter.filedialog.askdirectory(title = 'Please Select MR Dicom directory for conversion')
    return INPUT_DIR

def _slice_positions(slice_datasets):
    image_orientation = slice_datasets[0].ImageOrientationPatient
    row_cosine, column_cosine, slice_cosine = _extract_cosines(image_orientation)
    return [np.dot(slice_cosine, d.ImagePositionPatient) for d in slice_datasets]

def _extract_cosines(image_orientation):
    row_cosine = np.array(image_orientation[:3])
    column_cosine = np.array(image_orientation[3:])
    slice_cosine = np.cross(row_cosine, column_cosine)
    return row_cosine, column_cosine, slice_cosine


def check_for_missing_slices(slice_positions):
    slice_positions_diffs = np.diff(sorted(slice_positions))
    print('there does not appear to be any missing slices')
    if not np.allclose(slice_positions_diffs, slice_positions_diffs[0], atol=0, rtol=1e-5):
        # TODO: figure out how we should handle non-even slice spacing
        msg = "The slice spacing is non-uniform but there doesnt appear to be any missing. Slice spacings:\n{}"
        logger.warning(msg.format(slice_positions_diffs))
    if not np.allclose(slice_positions_diffs, slice_positions_diffs[0], atol=0, rtol=1e-1):
        raise DicomImportException('It appears there are missing slices')
        


def _slice_spacing(slice_datasets):
    try:
        slice_thickness = np.abs(slices[0].ImagePositionPatient[2] - slices[1].ImagePositionPatient[2])
    except:
        slice_thickness = np.abs(slices[0].SliceLocation - slices[1].SliceLocation)
    for s in slices:
        s.SliceThickness = slice_thickness
    if len(slice_datasets) > 1:
        slice_positions = _slice_positions(slice_datasets) 
        slice_positions_diffs = np.diff(sorted(slice_positions))
        return np.mean(slice_positions_diffs)
    else:
        return 0.0

# define our clear function
def clear():
        
    # for windows
    if name == 'nt':
        _ = system('cls')
    else:
        _ = system('clear')
        
clear()
print('created definitions and cleared temporary data')
class DicomImportException(Exception):
    pass
#%%
# GET dicom directory and print number of files found and sort
lstFilesDCM = [] # create an empty list
INPUT_DIR = get_INPUT_directory()
for dirName, subdirList, fileList in os.walk(INPUT_DIR):
    for filename in fileList:
        if ".dcm" or ".img" or "" in filename.lower():  # check whether the file's DICOM
            lstFilesDCM.append(os.path.join(dirName,filename))

# GET dicom directory and print number of files found and sort
lstFilesDCM2 = [] # create an empty list
INPUT_DIR2 = get_INPUT2_directory()
for dirName, subdirList, fileList in os.walk(INPUT_DIR2):
    for filename in fileList:
        if ".dcm" or ".img" or ""  in filename.lower():  # check whether the file's DICOM
            lstFilesDCM2.append(os.path.join(dirName,filename))
slices = [pydicom.read_file(filenameDCM, force = True) for filenameDCM in lstFilesDCM2]
print('list of files created')
try:
    slices.sort(key=lambda x: int(x.InstanceNumber))
except:
    slices = natsorted(slices)
#Find actual Slice Thickness
try:
    slice_thickness = np.abs(slices[0].ImagePositionPatient[2] - slices[1].ImagePositionPatient[2])
except:
    if "SliceLocation" in slices:
        slice_thickness = np.abs(slices[0].SliceLocation - slices[1].SliceLocation)
    else:
        slice_thickness = input(" can't determine slice thickness, please find and input manually now, only integer")
for s in slices:
    s.SliceThickness = slice_thickness  
print("found {"+format(len(slices))+"} Dicom Files")
print('Creating Pixel Array')
# Get ref file
RefDs = pydicom.dcmread(lstFilesDCM2[0], force = True)

#check for issues in the Dicom
try:
    slice_positions = _slice_positions(slices)
except:
    pass
try:
    check_for_missing_slices(slice_positions)
except:
    pass
try:
    spacing = _slice_spacing(slices)
except:
    pass
print('number of slices is ', len(slices))
try:
    print("Orientation: ", slices[0].ImageOrientationPatient)
except:
    print('no orientation provided')
try:
    print('spacing between slices written', RefDs.SpacingBetweenSlices)
except:
    print('no spacing obtained so calculating')
try:
    Actualspacing = spacing - slice_thickness
    print('Actual spacing between slices is:', Actualspacing)
except:
    pass
# In[103]:
try:
    print('manufacturer is', slices[0].Manufacturer)
except:
    print('no manufacturer stated')
try:
    print('kernel is ', slices[0].ConvolutionKernel)
except:
    print('No kernel Stated')
try:
    cos_value = (slices[0].ImageOrientationPatient[0])
    cos_degree = round(math.degrees(math.acos(cos_value)),2)
    print("cosign value_", cos_value)
    print("cosign degree_", cos_degree)
except:
    pass

print("Columns: ", RefDs.Columns)
print("Rows: ", RefDs.Rows)
# Load dimensions based on the number of rows, columns, and slices (along the Z axis)
ConstPixelDims = (int(RefDs.Rows), int(RefDs.Columns), len(slices))


# Load spacing values (in mm)
print('written pixel spacing', RefDs.PixelSpacing)
print('calculated pixel dimensions', slices[0].PixelSpacing)
pixel_spacing = slices[0].PixelSpacing
print('Written slice thickness', RefDs.SliceThickness)
print('calculated slice thickness')
pixel_spacing.append(slices[0].SliceThickness)
ActualThickness = slices[0].SliceThickness
print(ActualThickness)

ConstPixelSpacing = (float(RefDs.PixelSpacing[0]), float(RefDs.PixelSpacing[1]), float(ActualThickness))

x = np.arange(0.0, (ConstPixelDims[0]+1)*ConstPixelSpacing[0], ConstPixelSpacing[0])
y = np.arange(0.0, (ConstPixelDims[1]+1)*ConstPixelSpacing[1], ConstPixelSpacing[1])
z = np.arange(0.0, (ConstPixelDims[2]+1)*ConstPixelSpacing[2], ConstPixelSpacing[2])

# The array is sized based on 'ConstPixelDims'
ArrayDicom = np.zeros(ConstPixelDims, dtype=RefDs.pixel_array.dtype)

OUT_DIR = INPUT_DIR2 + "_Converted/"
if not os.path.exists(OUT_DIR):
    os.mkdir(OUT_DIR)
# loop through all the DICOM files


print("Setting file meta information...")
for filenameDCM in lstFilesDCM: 
    ds = pydicom.dcmread(filenameDCM, force = True)
    ds.decompress()
    ds.decode()
print("Setting dataset values...")
print('MR dicom To CT....')
for filenameDCM in lstFilesDCM2:
    rs = pydicom.dcmread(filenameDCM, force = True)
    rs.decompress()
    rs.decode()
    ArrayDicom[:, :, lstFilesDCM2.index(filenameDCM)] = rs.pixel_array
    var1 = rs.is_little_endian 
    var2 = rs.is_implicit_VR 
    ds.PixelData = ArrayDicom[:, :, lstFilesDCM2.index(filenameDCM)].tobytes()
    if var1 == True:
        ds.is_little_endian = True
        if var2 == True:
            ds.is_implicit_VR = True
            ds.file_meta.TransferSyntaxUID = pydicom.uid.ImplicitVRLittleEndian
        else:
            ds.is_implicit_VR = False
            ds.file_meta.TransferSyntaxUID = pydicom.uid.ExplicitVRLittleEndian
    else: 
        ds.is_little_endian = False
        if var2 == True:
            ds.is_implicit_VR = True
            ds.file_meta.TransferSyntaxUID = pydicom.uid.DeflatedExplicitVRLittleEndian
        else:
            ds.is_implicit_VR = False
            ds.file_meta.TransferSyntaxUID = pydicom.uid.ExplicitVRBigEndian
    try:
        NEWDCM = OUT_DIR+str(ds.InstanceNumber)+".dcm"
    except:
        NEWDCM = OUT_DIR+str(filenameDCM)
    ds.Modality = "CT" 
    try:
        ds.PhotometricInterpretation = rs.PhotometricInterpretation
        ds.PixelRepresentation = rs.PixelRepresentation
    except: 
        pass
    try:
        ds.WindowCenter = rs.WindowCenter
        ds.WindowWidth = rs.WindowWidth
    except:
        pass
    try:
        ds[0x28, 0x1055] = rs[0x28,0x1055]
    except:
        pass
    try: 
        ds.ConvolutionKernel = ds.ConvolutionKernel
    except:
        ds.ConvolutionKernel = 'BonePlus'
    try: 
        ds.Manufacturer = ds.Manufacturer
    except:
        ds.Manufacturer = 'GE MEDICAL SYSTEMS'
    ds.SpacingBetweenSlices = Actualspacing
    ds.SliceThickness = ActualThickness
    ds.Rows = rs.Rows
    ds.Columns = rs.Columns
    ds.PixelSpacing = rs.PixelSpacing
    ds.PatientName = rs.PatientName
    ds.PatientID = rs.PatientID
    ds.StudyID = rs.StudyID
    try:
        ds.FrameofReferenceUID = rs.FrameofReferenceUID
    except:
        pass
    ds.SeriesDescription = 'EssentialTremorPETMR2CT/MR'
    # Determine Slice Location and calculate if not present
    try:
        ds.SliceLocation = rs.SliceLocation
    except:
        if ("ImagePositionPatient" in rs)      and \
            rs.ImagePositionPatient            and \
            ("ImageOrientationPatient" in rs)  and \
            (len(rs.ImageOrientationPatient) >= 6):
            o = rs.ImageOrientationPatient
            slice_normals = [ (o[1] * o[5]) - (o[2] * o[4]),
                              (o[2] * o[3]) - (o[0] * o[5]),
                              (o[0] * o[4]) - (o[1] * o[3]),
                            ]
            ds.SliceLocation = sum([a * b for a, b in \
                            zip(slice_normals, rs.ImagePositionPatient)])
        else:
            print('cant determine slice location')
    try:
            ds.ImagePositionPatient = rs.ImagePositionPatient
            ds.ImageOrientationPatient = rs.ImageOrientationPatient
    except:
        pass
    ds.file_meta.MediaStorageSOPInstanceUID = rs.file_meta.MediaStorageSOPInstanceUID 
    try:
        ds.InstanceNumber = rs.InstanceNumber
    except:
        pass
    try:
        ds.SpecificCharacterSet = rs.SpecificCharacterSet
    except:
        pass
    try:
        ds.SOPInstanceUID = rs.SOPInstanceUID
    except:
        pass
    pydicom.filewriter.dcmwrite(NEWDCM,ds,write_like_original=False)
input("Press enter to close program")
clear()