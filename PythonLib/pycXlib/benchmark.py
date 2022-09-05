""" cellXpress library: Benchmark functions
    ============================================================================
    Written by Lit-Hsin Loo (c) 2022, A*STAR
"""
### Load the libraries
import scipy.optimize
import h5py    
import numpy as np 
import roifile
import cv2
import math
import tifffile
import cmapy
import pandas as pd
import pprint
import sys
import matplotlib.pyplot as plt

### Define common functions
def LoadcXSeg(seg_opt: dict) -> np.array:
    """
        Load cellXpress segmentation
    """
    ### Load the cellXpress segmentation
    with h5py.File(seg_opt['filename'], 'r') as seg_file:
        cX_segmask = seg_file['seg_out'][
            f'frame{seg_opt["frame_idx"]}-region0'
        ][...]

    return cX_segmask

def LoadcXMergedSeg(seg_opt: dict):
    """
        Load cellXpress Merged segmentation
    """

def LoadTIFFSeg(seg_opt: dict) -> np.array:
    """
        Load cell mask from TIFF
    """
    segmask = tifffile.imread(seg_opt['filename'])
    if segmask.dtype == np.float32:
        if np.max(segmask) > 65535:
            raise ValueError('segmask has more than 65535 cells')
        
        segmask = segmask.astype(np.uint16)
    return segmask

def GetCellMaskFromSegMask(seg_mask: np.array, cell_id: int):
    """
        Get the cell mask of the given cell from seg_mask
    """
    mkrset_num = seg_mask.shape[0]
    for mkrset_idx in range(mkrset_num):
        cell_mask = seg_mask[mkrset_idx, :, :] == cell_id
        if (np.count_nonzero(cell_mask)) > 0:
            break

    return cell_mask

def GetBBox(img: np.array) -> tuple[int, int, int, int]:
    """
        Find the bounding box
    """
    rows = np.any(img, axis=1)
    cols = np.any(img, axis=0)
    rmin, rmax = np.where(rows)[0][[0, -1]]
    cmin, cmax = np.where(cols)[0][[0, -1]]

    return rmin, rmax, cmin, cmax

def GenerateROIMask(roi: roifile.roifile.ImagejRoi, img_shape: tuple):
    """
    Generate ROI Mask
    """
    roi_mask = np.zeros(img_shape, np.uint8)
    contour = np.rint(roi.coordinates()).astype(int)
    cv2.drawContours(roi_mask, [contour], 0, 255, -1)
    return roi_mask

def DrawSegPerf(rois, segmask, roi_def):
    """
    Draw seg mask and color it with the performance
    ----------------------------------------------------------------------------
    """
    #img_out = cv2.normalize(segmask, dst=None, alpha=0, beta=255, norm_type=cv2.NORM_MINMAX)
    #img_out = cv2.convertScaleAbs(img_out)
    #img_out = cv2.cvtColor(img_out, cv2.COLOR_GRAY2BGR)
    
    #img_out = cv2.normalize(segmask, dst=None, alpha=0, beta=65525, norm_type=cv2.NORM_MINMAX)
    #cmap = plt.get_cmap("viridis", 2**16)
    #img_out = (cmap(img_out)* 2**8).astype(np.uint8)[:,:,:3]
    #img_out = cv2.cvtColor(img_out, cv2.COLOR_RGB2BGR)

    #cmap = plt.get_cmap("viridis", 256)
    scale_num = 256
    cmap = plt.get_cmap("viridis", scale_num)
    img_out = np.zeros((segmask.shape[0], segmask.shape[1], 4), np.uint8)

    ### Determine line width
    if (roi_def.right - roi_def.left + 1) > 1000:
        line_width = 2
    else:
        line_width = 1

    ### Draw the ground truth
    np.random.seed(0)
    for roi in rois:
        contour = np.rint(roi.coordinates()).astype(int)
        color = cmap(np.random.randint(0, scale_num))
        #color = (color[0], color[1], color[2], 0.5)
        color = tuple(np.floor(255*x) for x in color)
        cv2.drawContours(img_out, [contour], 0, color, -1)
    
    ### Draw the segmask
    cell_ids = FindCellIDs(segmask, roi_def)
    for cell_id in cell_ids:
        cell_img = (segmask == cell_id).astype(np.uint8)
        contours, hierarchy = cv2.findContours(
            cell_img, cv2.RETR_LIST, cv2.CHAIN_APPROX_SIMPLE)
        cv2.drawContours(img_out, contours, -1, (255, 0, 0, 255), line_width)

    img_out = cv2.cvtColor(img_out, cv2.COLOR_RGB2BGR)
    return img_out

def DrawROI(rois, segmask, roi_def):
    """
    Draw ROI over seg mask
    ----------------------------------------------------------------------------
    """
    #img_out = cv2.normalize(segmask, dst=None, alpha=0, beta=255, norm_type=cv2.NORM_MINMAX)
    #img_out = cv2.convertScaleAbs(img_out)
    #img_out = cv2.cvtColor(img_out, cv2.COLOR_GRAY2BGR)
    
    #img_out = cv2.normalize(segmask, dst=None, alpha=0, beta=65525, norm_type=cv2.NORM_MINMAX)
    #cmap = plt.get_cmap("viridis", 2**16)
    #img_out = (cmap(img_out)* 2**8).astype(np.uint8)[:,:,:3]
    #img_out = cv2.cvtColor(img_out, cv2.COLOR_RGB2BGR)

    #cmap = plt.get_cmap("viridis", 256)
    scale_num = 256
    cmap = plt.get_cmap("viridis", scale_num)
    img_out = np.zeros((segmask.shape[1], segmask.shape[2], 4), np.uint8)

    ### Determine line width
    if (roi_def.right - roi_def.left + 1) > 1000:
        line_width = 2
    else:
        line_width = 1

    ### Draw the ground truth
    alpha = 0.25
    np.random.seed(0)    
    #for roi in rois:
    for roi_idx, roi in enumerate(rois):
        contour = np.rint(roi.coordinates()).astype(int)
        color = cmap(np.random.randint(0, scale_num))
        color = (color[0], color[1], color[2], 1)
        color = tuple(np.floor(255*x) for x in color)
        
        cell_out = np.zeros((segmask.shape[1], segmask.shape[2], 4), np.uint8)
        cv2.drawContours(cell_out, [contour], 0, color, -1)
        #if roi_idx == 100:
        #    cv2.imwrite('figures3/img_out.png', img_out)
        #    cv2.imwrite('figures3/cell_out.png', cell_out)

        img_out = cv2.add(img_out, cell_out)
        #img_out = cv2.bitwise_xor(img_out, cell_out)
        #img_out = cv2.addWeighted(img_out, 1-alpha, cell_out, alpha, 0)
        #if roi_idx == 100:
        #    cv2.imwrite('figures3/img_out2.png', img_out)
        #    sys.exit('My error message')
        
    
    ### Draw the segmask
    cell_ids = FindCellIDs(segmask, roi_def)
    for cell_id in cell_ids:
        cell_mask = GetCellMaskFromSegMask(segmask, cell_id)
        cell_img = cell_mask.astype(np.uint8)
        contours, hierarchy = cv2.findContours(
            cell_img, cv2.RETR_LIST, cv2.CHAIN_APPROX_SIMPLE)
        cv2.drawContours(img_out, contours, -1, (255, 0, 0, 255), line_width)

    img_out = cv2.cvtColor(img_out, cv2.COLOR_RGB2BGR)
    return img_out

class BBoxClass:
    """
    Bounding box class
    """
    ### Variable
    left:float = math.inf
    top:float  = math.inf
    right:float = 0.0
    bottom:float = 0.0

    ### Properties
    #def __init__(self):
    
    ### Functions
    def __str__(self):
        return f"BBox: width = {self.right - self.left + 1}, height = {self.bottom - self.top + 1}"

def IsOverlap(roi1, roi2) -> bool:
    """
    Check if two rois are overlapped
    ---------------------------------------------------------------------------
    """
    return not (
        roi1.right < roi2.left or
        roi1.left > roi2.right or
        roi1.top < roi2.bottom or
        roi1.bottom > roi2.top
    #    self.top_right.x < other.bottom_left.x or 
    #    self.bottom_left.x > other.top_right.x or 
    #    self.top_right.y < other.bottom_left.y or 
    #    self.bottom_left.y > other.top_right.y
    )
    

def LoadIJROI(root_path: str, roi_name: str):
    """
    Load and process ImageJ ROI files
    ----------------------------------------------------------------------------
    """
    ### Construct the path name
    roi_path = f'{root_path}/{roi_name}_Cell.zip'

    ### Load the rois
    rois = roifile.ImagejRoi.fromfile(roi_path)

    ### The bouding box for all ROIs
    roi_all_bbox = BBoxClass()
    
    roi_def_idx = -1
    for roi_idx, roi in enumerate(rois):
        if (
            roi.name.startswith('ROI') and
            roi.roitype == roifile.ROI_TYPE.RECT and
            roi_def_idx < 0
        ):
            roi_def_idx = roi_idx
        elif roi.roitype != roifile.ROI_TYPE.RECT:
            ### Add to the bounding box
            roi_all_bbox.left   = min(roi.left, roi_all_bbox.left)
            roi_all_bbox.top    = min(roi.top, roi_all_bbox.top)
            roi_all_bbox.right  = max(roi.right, roi_all_bbox.right)
            roi_all_bbox.bottom = max(roi.bottom, roi_all_bbox.bottom)

    ### Can we find the ROI definition?
    roi_def = None

    if roi_def_idx < 0:
        print(f'Cannot find ROI specification in {roi_name}')
    else:
        ### roi_def is defined
        roi_def = rois[roi_def_idx]
        del rois[roi_def_idx]

        ### We need to remove rois not bounded by roi_def
        rois[:] = [roi for roi in rois if not IsOverlap(roi, roi_def)]
    
    return rois, roi_all_bbox, roi_def

def _safe_divide(x, y, eps=1e-10):
    """
        Computes a safe divide which returns 0 if y is zero
    """
    if np.isscalar(x) and np.isscalar(y):
        return x/y if np.abs(y)>eps else 0.0
    else:
        res = np.zeros(np.broadcast(x,y).shape, np.float32)
        np.divide(x,y, out=res, where=np.abs(y)>eps)
        return res

def precision(tp:float, fp:float, fn:float):
    return tp/(tp+fp) if tp > 0 else 0

def recall(tp:float, fp:float, fn:float):
    return tp/(tp+fn) if tp > 0 else 0

def accuracy(tp:float, fp:float, fn:float):    
    # -> https://www.kaggle.com/c/data-science-bowl-2018#evaluation
    return tp/(tp+fp+fn) if tp > 0 else 0

def f1(tp:float, fp:float, fn:float):
    return (2*tp)/(2*tp+fp+fn) if tp > 0 else 0

def FindCellIDs(segmask, roi_def):
    """
    Find all the cell IDs within roi_def
    """
    if not roi_def:
        ### Use the whole image
        cell_ids = np.unique(segmask)
    else:
        ### Only use roi def
        cell_ids = np.unique(
            segmask[
                :,
                roi_def.top:(roi_def.bottom+1),
                roi_def.left:(roi_def.right+1)
            ]
        )
    
    cell_ids = cell_ids[cell_ids > 0]

    return cell_ids

def FilterROI(mask, roi_def):
    """
    Filter a mask based on a roi definition
    ---------------------------------------------------------------------------
    """
    mask_roi = mask.copy()

    if len(mask.shape) == 3:        
        mkrset_num = mask.shape[0]
        for mkrset_idx in range(mkrset_num):
            mask_roi[mkrset_idx, :, :roi_def.left]      = 0  # the last index is not inclusive!!
            mask_roi[mkrset_idx, :, (roi_def.right+1):] = 0  # the first index is inclusive
            mask_roi[mkrset_idx, :roi_def.top, :]       = 0
            mask_roi[mkrset_idx, (roi_def.bottom+1):,:] = 0
    else:
        mask_roi[:, :roi_def.left]      = 0  # the last index is not inclusive!!
        mask_roi[:, (roi_def.right+1):] = 0  # the first index is inclusive
        mask_roi[:roi_def.top, :]       = 0
        mask_roi[(roi_def.bottom+1):,:] = 0

    return mask_roi

def ComputeSegPerf(rois, segmask, roi_def):
    """
    Compute segmentation performance
    * segmask is the segmentation results, which is an numpy array that may
    have 3 dimensions for merged cells (chn_num, width, height)
    * rois are the annotated cell boundaries
    * roi_def is the region where comparison should be made
    """
    ### Only focus on ROI
    segmask_roi = FilterROI(segmask, roi_def)

    ### Find all the cell_ids within the roi_def
    cell_ids = FindCellIDs(segmask_roi, None)
    #print(cell_ids)

    ### Define the output assignment matrix.  row = rois, col = cell_ids
    iou_mat = np.zeros((len(rois), len(cell_ids)), np.float64)

    ### Find the areas of all segmask
    cell_areas = []
    for cell_id in cell_ids:
        cell_areas.append(np.count_nonzero(segmask_roi == cell_id))
    
    ### Loop through all the annotated ROIs ("true cells")
    roi_areas = []
    for roi_idx, roi in enumerate(rois):

        ### Generate a mask for the ROI
        roi_mask = GenerateROIMask(roi, segmask_roi.shape[1:3])
        roi_mask = FilterROI(roi_mask, roi_def)

        ### The roi_area may be zero because outside of the ROI def
        roi_area = np.count_nonzero(roi_mask)
        roi_areas.append(roi_area)

        ### Only determine IOU for non-zero roi
        if roi_area > 0:
            ## Find all the predicted cells overlapped with the roi
            overlap_areas = dict()

            for mkrset_idx in range(segmask_roi.shape[0]):
                overlap_ids_cur, counts_cur = np.unique(
                    segmask_roi[mkrset_idx, roi_mask>0],
                    return_counts=True
                )
                overlap_areas_cur = dict(zip(overlap_ids_cur, counts_cur))
                if 0 in overlap_areas_cur.keys():
                    del overlap_areas_cur[0]

                overlap_areas |= overlap_areas_cur

                #print(f"mkrset_idx={mkrset_idx}")
                #print(overlap_areas)

            #roi_area = sum(overlap_areas.values())
                
            ### Loop through overlaps
            for overlap_id in overlap_areas.keys():
                ### We only check non-zero cell ids
                if overlap_id != 0:

                    ## The overlap cells may not be within roi_def
                    cell_idx = np.where(cell_ids == overlap_id)[0]

                    if cell_idx.size > 0:
                        cell_area    = cell_areas[cell_idx[0]]
                        overlap_area = overlap_areas[overlap_id]
                        union_area   = cell_area + roi_area - overlap_area
                        iou = overlap_area/ union_area
                        #print(f"roi_area     = {roi_area}")
                        #print(f"cell_area    = {cell_area}")
                        #print(f"overlap_area = {overlap_area}")
                        #print(f"union_area   = {union_area}")
                        #print(f"iou          = {iou}")
                        #if iou >= 0.5:
                        iou_mat[roi_idx, cell_idx] = iou

    #pd.DataFrame(iou_mat).to_csv("temp.csv")

    ### Remove ROIs with empty area
    print(f"Empty ROIs = {roi_areas.count(0)}")
    if roi_areas.count(0) > 0:
        nonempty_idx = [i for i, x in enumerate(roi_areas) if x > 0]

        rois = [rois[i] for i in nonempty_idx]
        roi_areas = [roi_areas[i] for i in nonempty_idx]
        iou_mat = iou_mat[nonempty_idx, :]

    ### Find the optimum matching based on iou
    matched_roi_idxs, matched_cell_idxs = scipy.optimize.linear_sum_assignment(
        iou_mat, True)
    #print(iou_mat.shape)
    #print(f"length of matched_roi_idxs  = {len(matched_roi_idxs)}")
    #print(f"length of matched_cell_idxs = {len(matched_cell_idxs)}")
    #roi_matches = dict(zip(matched_roi_idxs, matched_cell_idxs))

    ### Calculate sum_matched score
    sum_matched_score = iou_mat[matched_roi_idxs, matched_cell_idxs].sum()

    ### Create output image
    img_out = np.zeros((segmask.shape[1], segmask.shape[2], 4), np.uint8)

    scale_num = 256
    cmap = plt.get_cmap("viridis", scale_num)

    ### Draw the segmask (predicted cell region)
    area_ratio_log2_max = 2

    for cell_idx, cell_id in enumerate(cell_ids):
        ### Find the predicted area
        predicted_area  = cell_areas[cell_idx]

        ### Is there any match?
        match_idx = np.where(cell_idx == matched_cell_idxs)[0]

        if match_idx.size > 0:
            ### There is a match
            true_area = roi_areas[matched_roi_idxs[match_idx[0]]]
            area_ratio = np.log2(predicted_area/true_area)            
        else:
            ### There is no match, or predicted_area too large!!
            area_ratio = area_ratio_log2_max
        
        if area_ratio > area_ratio_log2_max:
            area_ratio = area_ratio_log2_max

        if area_ratio < (-area_ratio_log2_max):
            area_ratio = (-area_ratio_log2_max)

        color = cmap((area_ratio/area_ratio_log2_max)*(scale_num/2))
        color = tuple(np.floor(255*x) for x in color)

        ### Draw the predicted segmask        
        cell_img = GetCellMaskFromSegMask(segmask, cell_id).astype(np.uint8)
        contours, hierarchy = cv2.findContours(
            cell_img, cv2.RETR_LIST, cv2.CHAIN_APPROX_SIMPLE)
        cv2.drawContours(img_out, contours, 0, color, -1)

    img_out = cv2.cvtColor(img_out, cv2.COLOR_RGB2BGR)

    ### Calculate the accuracy
    nt = len(cell_ids)
    tp = len(matched_roi_idxs)
    fp = len(cell_ids) - tp
    fn = len(rois) - tp
    #print(f"tp = {tp}")

    # Compute performance metrics
    perf_scores = {
        'precision' : precision(tp, fp, fn),
        'recall'    : recall(tp, fp, fn),
        'accuracy'  : accuracy(tp, fp, fn),
        'f1'        : f1(tp, fp, fn),
        # the score average over all matched cells
        'mean_matched_score' : _safe_divide(sum_matched_score, tp),
        # the score average over all true cells
        'mean_true_score'    : _safe_divide(sum_matched_score, nt),
        'panoptic_quality'   : _safe_divide(sum_matched_score, tp+fp/2+fn/2)
    }
    
    return perf_scores, img_out