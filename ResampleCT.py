import pydicom
import pydicom.pixel_data_handlers.gdcm_handler as gdcm_handler
pydicom.config.image_handlers = [None, gdcm_handler]
from joblib import dump
import sys
import os
import numpy
import numpy as np
import natsort
from natsort import natsorted
from os import system, name
from tkinter import *
from tkinter import filedialog
import pydicom.uid 
from pydicom.uid import ExplicitVRBigEndian
from pydicom.uid import ExplicitVRLittleEndian 
from pydicom.uid import ImplicitVRLittleEndian
from pydicom.uid import DeflatedExplicitVRLittleEndian 
from pydicom.dicomdir import DicomDir
from pathlib import Path
import math
import scipy.ndimage
import matplotlib.pyplot as plt
import logging
from skimage import measure, morphology, exposure
from progressbar import *               # just a simple progress bar
import warnings
from tkinter import messagebox

warnings.filterwarnings('ignore', '.*output shape of zoom.*')
widgets = ['Test: ', Percentage(), ' ', Bar(marker='0',left='[',right=']'),
           ' ', ETA(), ' ', FileTransferSpeed()] #see docs for other options
def get_INPUT_directory():
    root = Tk()
    root.withdraw()
    INPUT_DIR = filedialog.askdirectory(title = 'Please Select CT Dicom directory Reference')
    return INPUT_DIR
# GET dicom directory and print number of files found and sort
lstFilesDCM = [] # create an empty list
INPUT_DIR = get_INPUT_directory()
""" Main Directory"""
OUT_DIR1 = INPUT_DIR + "NEW"
if not os.path.exists(OUT_DIR1):
    os.mkdir(OUT_DIR1)
class Logger(object):
    def __init__(self):
        DicomCheck = Path(OUT_DIR1 + '/DicomCheck.txt')
        DicomCheck.touch(exist_ok=True)  # will create file, if it exists will do nothing
        self.terminal = sys.stdout
        self.log = open(DicomCheck, 'w+')

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)  

    def flush(self):
        #this flush method is needed for python 3 compatibility.
        #this handles the flush command by doing nothing.
        #you might want to specify some extra behavior here.
        pass    

sys.stdout = Logger()
print('called all protocols')
def get_INPUT_directory():
    root = Tk()
    root.withdraw()
    INPUT_DIR = filedialog.askdirectory(title = 'Please Select CT Dicom directory Reference')
    return INPUT_DIR
def normalise_zero_one(image):
    """Image normalisation. Normalises image to fit [0, 1] range."""
    image = image.astype(np.float32)
    minimum = np.min(image)
    maximum = np.max(image)
    if maximum > minimum:
        ret = (image - minimum) / (maximum - minimum)
    else:
        ret = image * 0.
    return ret

def normalise_one_one(image):
    """Image normalisation. Normalises image to fit [-1, 1] range."""
    ret = normalise_zero_one(image)
    ret *= 2.
    ret -= 1.
    return ret
def get_INPUT2_directory():
    root = Tk()
    root.withdraw()
    INPUT_DIR = filedialog.askdirectory(title = 'Please Select MR Dicom directory for conversion')
    return INPUT_DIR

def vector(x1,y1,z1,x2,y2,z2):
    uv = [0,0,0]
    v = [x2-x1,y2-y1,z2-z1]
    len = math.sqrt(math.pow(x2-x1, 2) + math.pow(y2-y1, 2) +math.pow(z2-z1,2)*1.0)
    uv[0] = v[0]/len
    uv[1] =v[1]/len
    uv[2] = v[2]/len
    return uv
	
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
    if not np.allclose(slice_positions_diffs, slice_positions_diffs[0], atol=0, rtol=1e-5):
        # TODO: figure out how we should handle non-even slice spacing
        msg = "The slice spacing is non-uniform but there doesnt appear to be any missing. Slice spacings:\n{}"
        logger.warning(msg.format(slice_positions_diffs))
        s=sum((slice_positions_diffs))
        print('BSpline correction will be applied')
        return True, s
    if not np.allclose(slice_positions_diffs, slice_positions_diffs[0], atol=0, rtol=1e-1):
        print('It appears there are missing slices, Bspline applied')
        s=sum((slice_positions_diffs))
        return True, s
    
def is_compressed(self):
    """Return True if a compressed transfer syntax UID."""
    if self.is_transfer_syntax:
        # Explicit VR Little Endian
        # Implicit VR Little Endian
        # Explicit VR Big Endian
        # Deflated Explicit VR Little Endian
        if self in ['1.2.840.10008.1.2', '1.2.840.10008.1.2.1',
                    '1.2.840.10008.1.2.2', '1.2.840.10008.1.2.1.99']:
            return False
        print('Dicom is Compressed')
        # All encapsulated transfer syntaxes
        return True
    raise ValueError('UID is not a transfer syntax.')

def _slice_spacing(slice_datasets):
    try:
        slice_thickness = np.abs(slices[0].ImagePositionPatient[2] - slices[1].ImagePositionPatient[2])
    except:
        slice_thickness = np.abs(slices[0].SliceLocation - slices[1].SliceLocation)
    for s in slices:
        s.SliceThickness = slice_thickness
    if len(slice_datasets) > 1:
        slice_positions = _slice_positions(slice_datasets) 
        Min_S = min(_slice_positions(slice_datasets))
        Max_S = max(_slice_positions(slice_datasets))
        R = Max_S - Min_S
        slice_positions_diffs = np.diff(sorted(slice_positions))
        return R, np.mean(slice_positions_diffs)
    else:
        R = len(slice_datasets)*slice_datasets.SliceThickness
        return R, 0.0

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
for dirName, subdirList, fileList in os.walk(INPUT_DIR):
    for filename in fileList:
        if ".dcm" or ".img" or "" in filename.lower():  # check whether the file's DICOM
            if filename == 'DICOMDIR':
                continue
            lstFilesDCM.append(os.path.join(dirName,filename))
slices = [pydicom.dcmread(filenameDCM, force = True) for filenameDCM in lstFilesDCM]
S = 0
try:
    slices.sort(key=lambda x: int(x.InstanceNumber))
    S = 1
    print('sorting by instancenumber', S)
except:
    slices.sort(key = lambda x: float(x.ImagePositionPatient[2]))
    S = 1
    print('sorting by Image Position Patient', S)
if S == 1:
    pass
else:
    print('could not find an easy way of sorting the slices but should still work check order in radiant')

# Get ref file
RefDs = pydicom.dcmread(lstFilesDCM[0], force = True)
var = is_compressed(RefDs.file_meta.TransferSyntaxUID)
x1=slices[0].ImagePositionPatient[0]
y1=slices[0].ImagePositionPatient[1]
z1=slices[0].ImagePositionPatient[2]
x2=slices[1].ImagePositionPatient[0]
y2=slices[1].ImagePositionPatient[1]
z2=slices[1].ImagePositionPatient[2]
vector = vector(x1,y1,z1,x2,y2,z2)
print('the unit vector of the slices is :', vector)
""" Corrected Header Output"""
OUT_DIR = OUT_DIR1 + "/_HeaderCorrected/"
if not os.path.exists(OUT_DIR):
    os.mkdir(OUT_DIR)
""" Isomorphic Output"""
OUT_DIR2 = OUT_DIR1 +"/_Isomorphic/"
if not os.path.exists(OUT_DIR2):
    os.mkdir(OUT_DIR2)

if (var == True):
    OUT_DIR0 = OUT_DIR1 +"/_Decompressed/"
    if not os.path.exists(OUT_DIR0):
        os.mkdir(OUT_DIR0)
print(INPUT_DIR)
print('list of files created')

print("found {"+format(len(slices))+"} Dicom Files")
#Find actual Slice Thickness
try:
    slice_thickness = np.abs(slices[0].ImagePositionPatient[2] - slices[1].ImagePositionPatient[2])
except:
    if "SliceLocation" in slices:
        slice_thickness = np.abs(slices[0].SliceLocation - slices[1].SliceLocation)
for s in slices:
    s.SliceThickness = slice_thickness  
print('Creating Pixel Array')

#check for issues in the Dicom
try:
    slice_positions = _slice_positions(slices)
except:
    pass
try:
    va, s = check_for_missing_slices(slice_positions)
    print('Does the Dicom need to be resampled', va)
    print('number of slices need after correction', s)
except:
    pass
try:
    R, spacing = _slice_spacing(slices)
except:
    pass
print('number of slices right now is ', len(slices))
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
    Actualspacing = 0
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
ArrayDicomCor = np.zeros(ConstPixelDims, dtype=RefDs.pixel_array.dtype)
ArrayDicomISO = np.zeros(ConstPixelDims, dtype=RefDs.pixel_array.dtype)
print('A new folder containing the CT with corrected header will be made')
# loop through all the DICOM files
print("Correcting file meta information and exporting")
pbar = ProgressBar(widgets=widgets, maxval=len(lstFilesDCM))
pbar.start()
count = 0
for filenameDCM in lstFilesDCM: 
    count += 1
    ds = pydicom.dcmread(filenameDCM, force = True)
    if (var == True):
        ds.decompress()
        try:
            NEWDCM = OUT_DIR0+str(ds.InstanceNumber)+".dcm"
            pydicom.filewriter.dcmwrite(NEWDCM,ds,write_like_original=False)
        except:
            NEWDCM = OUT_DIR0+str(filenameDCM)
            pydicom.filewriter.dcmwrite(NEWDCM,ds,write_like_original=False)
    ds.decode()
    intercept = ds.RescaleIntercept
    slope = ds.RescaleSlope 
    ds.InstanceNumber = count
    ArrayDicom[:, :, lstFilesDCM.index(filenameDCM)] = ds.pixel_array
    var1 = ds.is_little_endian 
    var2 = ds.is_implicit_VR 
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
    ds.SpacingBetweenSlices = Actualspacing
    ds.SliceThickness = ActualThickness
    ds.SeriesDescription = 'HeaderCorrected'
    # Determine Slice Location and calculate if not present
    try:
        ds.SliceLocation = ds.SliceLocation
    except:
        if ("ImagePositionPatient" in ds)      and \
            ds.ImagePositionPatient            and \
            ("ImageOrientationPatient" in ds)  and \
            (len(ds.ImageOrientationPatient) >= 6):
            o = ds.ImageOrientationPatient
            slice_normals = [ (o[1] * o[5]) - (o[2] * o[4]),
                              (o[2] * o[3]) - (o[0] * o[5]),
                              (o[0] * o[4]) - (o[1] * o[3]),
                            ]
            ds.SliceLocation = sum([a * b for a, b in \
                            zip(slice_normals, ds.ImagePositionPatient)])
        else:
            print('cant determine slice location')
    #ArrayDicom[:, :, lstFilesDCM.index(filenameDCM)] = normalise_one_one(ArrayDicom[:, :, lstFilesDCM.index(filenameDCM)])
    ds.PixelData = ArrayDicom[:, :, lstFilesDCM.index(filenameDCM)].tobytes()
    intercept = slices[0].RescaleIntercept
    slope = slices[0].RescaleSlope
    if slope != 1:
        ArrayDicom[:, :, lstFilesDCM.index(filenameDCM)] = slope * ArrayDicom[:, :, lstFilesDCM.index(filenameDCM)].astype(np.float64)
        ArrayDicom[:, :, lstFilesDCM.index(filenameDCM)] = ArrayDicom[:, :, lstFilesDCM.index(filenameDCM)]
        ArrayDicom[:, :, lstFilesDCM.index(filenameDCM)] += np.int16(intercept)
    pydicom.filewriter.dcmwrite(NEWDCM,ds,write_like_original=False)
    pbar.update(count)
pbar.finish()
if (ActualThickness <= 0.999):
    Resampleme = 1
    print('A new folder containing a 1 mm reformat will be made since it is less that 1mm')
    OUT_DIR4 = OUT_DIR1 +"/_1mm_Reformat/"
    if not os.path.exists(OUT_DIR4):
        os.mkdir(OUT_DIR4)
    lstFilesDCM2 =[]
    for dirName, subdirList, fileList in os.walk(OUT_DIR):
        for filename in fileList:
            if ".dcm" or ".img" or "" in filename.lower():  # check whether the file's DICOM
                lstFilesDCM2.append(os.path.join(dirName,filename))
    ps = [pydicom.read_file(filenameDCM, force = True) for filenameDCM in lstFilesDCM2]
    ps.sort(key = lambda x: float(x.ImagePositionPatient[2]))
    print("Shape before resampling\t", ArrayDicom.shape)
    pbar = ProgressBar(widgets=widgets, maxval= 7)
    pbar.start()
    count = 0
    spacing = np.array(list(ps[0].PixelSpacing) + [ps[0].SliceThickness], dtype=np.float32)
    count = 1
    pbar.update(count)
    new_spacing = [ps[0].PixelSpacing[0],ps[0].PixelSpacing[1],1]
    count = 2
    pbar.update(count)
    resize_factor = spacing / new_spacing
    count = 3
    pbar.update(count)
    new_real_shape = ArrayDicom.shape * resize_factor
    new_shape = np.round(new_real_shape)
    count = 4
    pbar.update(count)
    real_resize_factor = new_shape / ArrayDicom.shape
    new_spacing = spacing / real_resize_factor
    count = 5
    pbar.update(count)
    ArrayDicomCor = scipy.ndimage.interpolation.zoom(ArrayDicom, real_resize_factor, order=5, mode='nearest')
    count = 6
    pbar.update(count)
    lstFilesDCM2 = lstFilesDCM2[0:ArrayDicomCor.shape[2]]
    count = 7
    pbar.update(count)
    pbar.finish()
    print("Shape after resampling\t", ArrayDicomCor.shape)
    
    pbar = ProgressBar(widgets=widgets, maxval=len(lstFilesDCM2))
    pbar.start()
    count = 0
    BaseSL = RefDs.SliceLocation
    BaseIPP = RefDs.ImagePositionPatient
    for filenameDCM in lstFilesDCM2: 
        rs = pydicom.dcmread(filenameDCM, force = True)
        rs.decompress()
        rs.decode()
        voxelH = (ps[0].SpacingBetweenSlices+ps[0].SliceThickness)*(ArrayDicom.shape[2]/ArrayDicomCor.shape[2])
        rs.SliceThickness =ps[0].SliceThickness*(ArrayDicom.shape[2]/ArrayDicomCor.shape[2]) 
        rs.SpacingBetweenSlices = ps[0].SpacingBetweenSlices*(ArrayDicom.shape[2]/ArrayDicomCor.shape[2])
        rs.PixelSpacing[0] = ps[0].PixelSpacing[0]*(ps[0].Columns/rs.Columns)
        rs.PixelSpacing[1] = ps[0].PixelSpacing[1]*(ps[0].Rows/rs.Rows)
        rs.SliceLocation = BaseSL
        rs.ImagePositionPatient = BaseIPP
        BaseSL += (voxelH)*vector[2]
        BaseIPP[0] += (voxelH)*vector[0]
        BaseIPP[1] += (voxelH)*vector[1]
        BaseIPP[2] += (voxelH)*vector[2]
        var1 = rs.is_little_endian 
        var2 = rs.is_implicit_VR 
        if var1 == True:
            rs.is_little_endian = True
            if var2 == True:
                rs.is_implicit_VR = True
                rs.file_meta.TransferSyntaxUID = pydicom.uid.ImplicitVRLittleEndian
            else:
                rs.is_implicit_VR = False
                rs.file_meta.TransferSyntaxUID = pydicom.uid.ExplicitVRLittleEndian
        else: 
            rs.is_little_endian = False
            if var2 == True:
                rs.is_implicit_VR = True
                rs.file_meta.TransferSyntaxUID = pydicom.uid.DeflatedExplicitVRLittleEndian
            else:
                rs.is_implicit_VR = False
                rs.file_meta.TransferSyntaxUID = pydicom.uid.ExplicitVRBigEndian
        count += 1
        rs.InstanceNumber = count
        
        try:
            NEWDCM4 = OUT_DIR4+str(rs.InstanceNumber)+".dcm"
        except:
            NEWDCM4 = OUT_DIR4+str(filenameDCM)
        rs.PixelData = ArrayDicomCor[:, :, lstFilesDCM2.index(filenameDCM)].tobytes()
        rs.SeriesDescription = '1mm Reformat'
        pydicom.filewriter.dcmwrite(NEWDCM4,rs,write_like_original=False)
        intercept = slices[0].RescaleIntercept
        slope = slices[0].RescaleSlope
        if slope != 1:
            ArrayDicomCor[:, :, lstFilesDCM2.index(filenameDCM)] = slope * ArrayDicomCor[:, :, lstFilesDCM2.index(filenameDCM)].astype(np.float64)
            ArrayDicomCor[:, :, lstFilesDCM2.index(filenameDCM)] = ArrayDicomCor[:, :, lstFilesDCM2.index(filenameDCM)]
            ArrayDicomCor[:, :, lstFilesDCM2.index(filenameDCM)] += np.int16(intercept)
        pbar.update(count)
    pbar.finish()
else:
    Resampleme = 0
Isomorphic = messagebox.askyesno("Create an Isomorphic  Image Set", "This will try to compensate for differences in CT Spatial Resolution")
if Isomorphic == True:
	print('re-read the dicom files to confirm writing and resample to isomorphic')
	lstFilesDCM3 =[]
	for dirName, subdirList, fileList in os.walk(OUT_DIR):
		for filename in fileList:
			if ".dcm" or ".img" or "" in filename.lower():  # check whether the file's DICOM
				lstFilesDCM3.append(os.path.join(dirName,filename))
	ds = [pydicom.read_file(filenameDCM, force = True) for filenameDCM in lstFilesDCM3]
	ds.sort(key = lambda x: float(x.ImagePositionPatient[2]))
	print("Shape before resampling\t", ArrayDicom.shape)
	spacing = np.array(list(ds[0].PixelSpacing) + [ds[0].SliceThickness], dtype=np.float32)
	new_spacing = [1,1,1]
	resize_factor = spacing / new_spacing
	new_real_shape = ArrayDicom.shape * resize_factor
	new_shape = np.round(new_real_shape)
	real_resize_factor = new_shape / ArrayDicom.shape
	new_spacing = spacing / real_resize_factor
	ArrayDicomISO = scipy.ndimage.interpolation.zoom(ArrayDicom, real_resize_factor, order=5, mode='nearest')
	lstFilesDCM3 = lstFilesDCM3[0:ArrayDicomISO.shape[2]]
	print("Shape after resampling\t", ArrayDicomISO.shape)
	pbar = ProgressBar(widgets=widgets, maxval=len(lstFilesDCM3))
	pbar.start()
	count = 0
	BaseSL = RefDs.SliceLocation
	BaseIPP = RefDs.ImagePositionPatient
	for filenameDCM in lstFilesDCM3: 
		iso = pydicom.dcmread(filenameDCM, force = True)
		iso.decompress()
		iso.decode()
		iso.Columns = ArrayDicomISO.shape[0]
		iso.Rows = ArrayDicomISO.shape[1]
		voxelH = (ds[0].SpacingBetweenSlices+ds[0].SliceThickness)*(ArrayDicom.shape[2]/ArrayDicomISO.shape[2])
		iso.SliceThickness =ds[0].SliceThickness*(ArrayDicom.shape[2]/ArrayDicomISO.shape[2]) 
		iso.SpacingBetweenSlices = ds[0].SpacingBetweenSlices*(ArrayDicom.shape[2]/ArrayDicomISO.shape[2])
		iso.PixelSpacing[0] = ds[0].PixelSpacing[0]*(ds[0].Columns/iso.Columns)
		iso.PixelSpacing[1] = ds[0].PixelSpacing[1]*(ds[0].Rows/iso.Rows)
		iso.SliceLocation = BaseSL
		iso.ImagePositionPatient = BaseIPP
		BaseSL += (voxelH)*vector[2]
		BaseIPP[0] += (voxelH)*vector[0]
		BaseIPP[1] += (voxelH)*vector[1]
		BaseIPP[2] += (voxelH)*vector[2]
		count += 1
		iso.InstanceNumber = count
		var1 = iso.is_little_endian 
		var2 = iso.is_implicit_VR 
		if var1 == True:
			iso.is_little_endian = True
			if var2 == True:
				iso.is_implicit_VR = True
				iso.file_meta.TransferSyntaxUID = pydicom.uid.ImplicitVRLittleEndian
			else:
			   iso.is_implicit_VR = False
			   iso.file_meta.TransferSyntaxUID = pydicom.uid.ExplicitVRLittleEndian
		else: 
			iso.is_little_endian = False
			if var2 == True:
				iso.is_implicit_VR = True
				iso.file_meta.TransferSyntaxUID = pydicom.uid.DeflatedExplicitVRLittleEndian
			else:
				iso.is_implicit_VR = False
				iso.file_meta.TransferSyntaxUID = pydicom.uid.ExplicitVRBigEndian
		try:
			NEWDCM2 = OUT_DIR2+str(iso.InstanceNumber)+".dcm"
		except:
			NEWDCM2 = OUT_DIR2+str(filenameDCM)
		iso.PixelData = ArrayDicomISO[:, :, lstFilesDCM3.index(filenameDCM)].tobytes()
		iso.SeriesDescription = 'Isomorphic Resampling'
		pydicom.filewriter.dcmwrite(NEWDCM2,iso,write_like_original=False)
		intercept = iso.RescaleIntercept
		slope = iso.RescaleSlope
		if slope != 1:
			ArrayDicomISO[:, :, lstFilesDCM3.index(filenameDCM)] = slope * ArrayDicomISO[:, :, lstFilesDCM3.index(filenameDCM)].astype(np.float64)
			ArrayDicomISO[:, :, lstFilesDCM3.index(filenameDCM)] = ArrayDicomISO[:, :, lstFilesDCM3.index(filenameDCM)]
			ArrayDicomISO[:, :, lstFilesDCM3.index(filenameDCM)] += np.int16(intercept)
		pbar.update(count)
	pbar.finish()
	print('checking for comparison need')
result = messagebox.askyesno("Comparison To another CT","Do you want to compare CT values to another CT scan (this section will not correct the new ct)")
if result == True:
	print('you selected yes')
	lstFilesDCM4 = [] # create an empty list
	INPUT_DIR2 = get_INPUT2_directory()
	for dirName, subdirList, fileList in os.walk(INPUT_DIR2):
		for filename in fileList:
			if ".dcm" or ".img" or ""  in filename.lower():  # check whether the file's DICOM
				lstFilesDCM4.append(os.path.join(dirName,filename))
	slices = [pydicom.read_file(filenameDCM, force = True) for filenameDCM in lstFilesDCM4]
	try:
		slice_thickness = np.abs(slices[1].ImagePositionPatient[2] - slices[2].ImagePositionPatient[2])
	except:
		if "SliceLocation" in slices:
			slice_thickness = np.abs(slices[1].SliceLocation - slices[2].SliceLocation)
	for s in slices:
		s.SliceThickness = slice_thickness  
	print('Creating Pixel Array')
	fDs = pydicom.dcmread(lstFilesDCM4[0], force = True)
	ConstPixelDims = (int(fDs.Rows), int(fDs.Columns), len(slices))

	ConstPixelSpacing = (float(slices[1].PixelSpacing[0]), float(slices[1].PixelSpacing[1]), float(slice_thickness))

	x = np.arange(0.0, (ConstPixelDims[0]+1)*ConstPixelSpacing[0], ConstPixelSpacing[0])
	y = np.arange(0.0, (ConstPixelDims[1]+1)*ConstPixelSpacing[1], ConstPixelSpacing[1])
	z = np.arange(0.0, (ConstPixelDims[2]+1)*ConstPixelSpacing[2], ConstPixelSpacing[2])
	# The array is sized based on 'ConstPixelDims'
	ArrayDicom4 = np.zeros(ConstPixelDims, dtype=fDs.pixel_array.dtype)
	for filenameDCM in lstFilesDCM4: 
		com = pydicom.dcmread(filenameDCM, force = True)
		com.decompress()
		com.decode()
		ArrayDicom4[:, :, lstFilesDCM4.index(filenameDCM)] = com.pixel_array
		intercept = com.RescaleIntercept
		slope = com.RescaleSlope
		if slope != 1:
			ArrayDicom4[:, :, lstFilesDCM4.index(filenameDCM)] = slope * ArrayDicom4[:, :, lstFilesDCM4.index(filenameDCM)].astype(np.float64)
			ArrayDicom4[:, :, lstFilesDCM4.index(filenameDCM)] = ArrayDicom4[:, :, lstFilesDCM4.index(filenameDCM)]
			ArrayDicom4[:, :, lstFilesDCM4.index(filenameDCM)] += np.int16(intercept)
TEMP = messagebox.askyesno("Template","Create a Template")
if TEMP == True:
	OUT_DIR5 = OUT_DIR1 +"/_BlankCTTEMPLATE/"
	if not os.path.exists(OUT_DIR5):
		os.mkdir(OUT_DIR5)
	print('you selected yes')
	lstFilesDCM5 = [] # create an empty list
	for dirName, subdirList, fileList in os.walk(OUT_DIR):
		for filename in fileList:
			if ".dcm" or ".img" or ""  in filename.lower():  # check whether the file's DICOM
				lstFilesDCM5.append(os.path.join(dirName,filename))
	slices = [pydicom.read_file(filenameDCM, force = True) for filenameDCM in lstFilesDCM5]
	try:
		slice_thickness = np.abs(slices[1].ImagePositionPatient[2] - slices[2].ImagePositionPatient[2])
	except:
		if "SliceLocation" in slices:
			slice_thickness = np.abs(slices[1].SliceLocation - slices[2].SliceLocation)
	for s in slices:
		s.SliceThickness = slice_thickness  
	print('Creating Pixel Array')
	fDs = pydicom.dcmread(lstFilesDCM5[0], force = True)
	ConstPixelDims = (int(fDs.Rows), int(fDs.Columns), len(slices))

	ConstPixelSpacing = (float(slices[1].PixelSpacing[0]), float(slices[1].PixelSpacing[1]), float(slice_thickness))

	x = np.arange(0.0, (ConstPixelDims[0]+1)*ConstPixelSpacing[0], ConstPixelSpacing[0])
	y = np.arange(0.0, (ConstPixelDims[1]+1)*ConstPixelSpacing[1], ConstPixelSpacing[1])
	z = np.arange(0.0, (ConstPixelDims[2]+1)*ConstPixelSpacing[2], ConstPixelSpacing[2])
	# The array is sized based on 'ConstPixelDims'
	ArrayDicom5 = np.zeros(ConstPixelDims, dtype=fDs.pixel_array.dtype)
	for filenameDCM in lstFilesDCM5: 
		com = pydicom.dcmread(filenameDCM, force = True)
		com.decompress()
		com.decode()
		com.PixelData=ArrayDicom5[:, :, lstFilesDCM5.index(filenameDCM)].tobytes
		intercept = com.RescaleIntercept
		slope = com.RescaleSlope
		NEWDCM5 = str(com.InstanceNumber)+'DS.pkl'
		dump(com,NEWDCM5)

print('Creating a CT value histogram for data analysis')
fig = plt.hist(ArrayDicom.flatten(), bins=1000, range = [400 , 2500], color='c', label='Original')
if Resampleme == 1:
	fig = plt.hist(ArrayDicomCor.flatten(), bins=1000, range = [400 , 2500], color='g', label='1mm_Resample')
else:
	pass
if Isomorphic == True:
    fig = plt.hist(ArrayDicomISO.flatten(), bins=1000, range = [400 , 2500], color='b', label='Isomorphic Resample')
if result == True:
    fig = plt.hist(ArrayDicom4.flatten(), bins=1000, range = [400 , 2500], color='b', label='ComparableCT')
plt.legend(loc='upper left')
plt.title('CT Value Histogram')
plt.xlabel("Hounsfield Units (HU)")
plt.ylabel("Frequency")
print('The output Directory is', OUT_DIR1)
print("Press enter to close program")
plt.savefig(OUT_DIR1 + '/HistogramFig.png')
input()
plt.close("all")
clear()
