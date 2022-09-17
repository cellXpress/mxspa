import sys
from pathlib import Path

### Setup the paths
repo_root = str(Path(__file__).absolute().parents[2])

python_lib_root = str(Path(__file__).absolute().parents[2] / 'PythonLib')
sys.path.append(python_lib_root)

### Load the libraries
from pycXlib.benchmark import *

### Script root
bench_root = str(Path(__file__).absolute().parent)

### Manual annotation data folder
annot_dir = f'{repo_root}/Reference_data/Annotations'

### cellXpress project root folder
### The four reference projects should be installed under this folder 
cX_proj_dir = str(Path(__file__).absolute().parents[3])

### Define the benchmark options
benchmarks = [
    ### 0. BEAS2B Lung Cells B05_p2
    {
        'name': 'BEAS2B_Lung_Cells',
        'annot_root_path' : f'{annot_dir}/2022_BEAS2B_Lung_Cells',
        'roi_names' : [
            'JY15135-B05_p2'
        ],
        ### Segmentations
        'seg_opts' :[
            {
                'type': 'cellXpress1',
                'filename' : f'{cX_proj_dir}/2022_07_cX2_BEAS2B_Lung_Cells/220726_BEAS2B_Lung_Cells_Project/markerset0/segmentation/0/seg28.cXd',
                'frame_idx': 1
            },
            {
                'type': 'cellXpress2',
                'filename' : f'{cX_proj_dir}/2022_07_cX2_BEAS2B_Lung_Cells/220726_BEAS2B_Lung_Cells_Project/markerset1/segmentation/0/seg28.cXd',
                'frame_idx': 1
            },
            {
                'type': 'qupath',
                'filename' : f'{bench_root}/Software_Projects/QuPath/1_BEAS2B_Lung_Cells/rois_final/JY15135-B05_p2t0.ome.tif - Image0-labels.ome.tif'
            },
            {
                'type': 'deepcell',
                'filename' : f'{bench_root}/Software_Projects/DeepCell/results_final/BEAS2B_B05_p2.tif'
            },
            {
                'type': 'cellprofiler',
                'filename' : f'{bench_root}/Software_Projects/CellProfiler/results_final/JY15135-B05_p2t0.ome_CellMask.tiff'
            }
        ]
    },
    ### 1. BEAS2B Lung Cells B05_p3
    {
        'name': 'BEAS2B_Lung_Cells',
        'annot_root_path' : f'{annot_dir}/2022_BEAS2B_Lung_Cells',
        'roi_names' : [
            'JY15135-B05_p3'
        ],
        ### Segmentations
        'seg_opts' :[
            {
                'type': 'cellXpress1',
                'filename' : f'{cX_proj_dir}/2022_07_cX2_BEAS2B_Lung_Cells/220726_BEAS2B_Lung_Cells_Project/markerset0/segmentation/0/seg28.cXd',
                'frame_idx': 2
            },
            {
                'type': 'cellXpress2',
                'filename' : f'{cX_proj_dir}/2022_07_cX2_BEAS2B_Lung_Cells/220726_BEAS2B_Lung_Cells_Project/markerset1/segmentation/0/seg28.cXd',
                'frame_idx': 2
            },
            {
                'type': 'qupath',
                'filename' : f'{bench_root}/Software_Projects/QuPath/1_BEAS2B_Lung_Cells/rois_final/JY15135-B05_p3t0.ome.tif - Image0-labels.ome.tif'
            },
            {
                'type': 'deepcell',
                'filename' : f'{bench_root}/Software_Projects/DeepCell/results_final/BEAS2B_B05_p3.tif'
            },
            {
                'type': 'cellprofiler',
                'filename' : f'{bench_root}/Software_Projects/CellProfiler/results_final/JY15135-B05_p3t0.ome_CellMask.tiff'
            }
        ]
    },
    ### 2. BEAS2B Lung Cells C11_p2
    {
        'name': 'BEAS2B_Lung_Cells',
        'annot_root_path' : f'{annot_dir}/2022_BEAS2B_Lung_Cells',
        'roi_names' : [
            'JY15135-C11_p2'
        ],
        ### Segmentations
        'seg_opts' :[
            {
                'type': 'cellXpress1',
                'filename' : f'{cX_proj_dir}/2022_07_cX2_BEAS2B_Lung_Cells/220726_BEAS2B_Lung_Cells_Project/markerset0/segmentation/0/seg58.cXd',
                'frame_idx': 1
            },
            {
                'type': 'cellXpress2',
                'filename' : f'{cX_proj_dir}/2022_07_cX2_BEAS2B_Lung_Cells/220726_BEAS2B_Lung_Cells_Project/markerset1/segmentation/0/seg58.cXd',
                'frame_idx': 1
            },
            {
                'type': 'qupath',
                'filename' : f'{bench_root}/Software_Projects/QuPath/1_BEAS2B_Lung_Cells/rois_final/JY15135-C11_p2t0.ome.tif - Image0-labels.ome.tif'
            },
            {
                'type': 'deepcell',
                'filename' : f'{bench_root}/Software_Projects/DeepCell/results_final/BEAS2B_C11_p2.tif'
            },
            {
                'type': 'cellprofiler',
                'filename' : f'{bench_root}/Software_Projects/CellProfiler/results_final/JY15135-C11_p2t0.ome_CellMask.tiff'
            }
        ]
    },
    ### 3. BEAS2B Lung Cell C11_p3
    {
        'name': 'BEAS2B_Lung_Cells',
        'annot_root_path' : f'{annot_dir}/2022_BEAS2B_Lung_Cells',
        'roi_names' : [
            'JY15135-C11_p3'
        ],
        ### Segmentations
        'seg_opts' :[
            {
                'type': 'cellXpress1',
                'filename' : f'{cX_proj_dir}/2022_07_cX2_BEAS2B_Lung_Cells/220726_BEAS2B_Lung_Cells_Project/markerset0/segmentation/0/seg58.cXd',
                'frame_idx': 2
            },
            {
                'type': 'cellXpress2',
                'filename' : f'{cX_proj_dir}/2022_07_cX2_BEAS2B_Lung_Cells/220726_BEAS2B_Lung_Cells_Project/markerset1/segmentation/0/seg58.cXd',
                'frame_idx': 2
            },
            {
                'type': 'qupath',
                'filename' : f'{bench_root}/Software_Projects/QuPath/1_BEAS2B_Lung_Cells/rois_final/JY15135-C11_p3t0.ome.tif - Image0-labels.ome.tif'
            },
            {
                'type': 'deepcell',
                'filename' : f'{bench_root}/Software_Projects/DeepCell/results_final/BEAS2B_C11_p3.tif'
            },
            {
                'type': 'cellprofiler',
                'filename' : f'{bench_root}/Software_Projects/CellProfiler/results_final/JY15135-C11_p3t0.ome_CellMask.tiff'
            }
        ]
    },
    ### 4. Lung LUNG IO TMA P#1A_Core[1,4,B]
    {
        ### Annotations (ground truth)
        'name': 'Lung_Tumors',
        'annot_root_path' : f'{annot_dir}/2022_Lung_Tumors',
        'roi_names' : [
            'LUNG IO TMA P#1A_Core[1,4,B]_[47045,13477]_p0t0_ROI1'
        ],
        ### Segmentations
        'seg_opts' :[
            {
                'type': 'cellXpress1',
                'filename' : f'{cX_proj_dir}/2022_07_cX2_Lung_Tumors/220203_Lung_Tumors_Single_Project/markerset0/segmentation/0/seg0.cXd',
                'frame_idx': 0
            },
            {
                'type': 'cX2_merged',
                'filename' : f'{bench_root}/Software_Projects/cX2/results_final/Lung_tumors_[1,4,B]_merged_cell_masks.tif',
            },
            {
                'type': 'qupath',
                'filename' : f'{bench_root}/Software_Projects/QuPath/3_Lung_Tumors/rois_final/1A_Core[1,4,B]_[47045,13477]_p0t0-labels.ome.tif'
            },
            {
                'type': 'deepcell',
                'filename' : f'{bench_root}/Software_Projects/DeepCell/results_final/1A_Core[1,4,B]_mask.tif'
            },
            {
                'type': 'cellprofiler',
                'filename' : f'{bench_root}/Software_Projects/CellProfiler/results_final/1A_Core[1,4,B]_[47045,13477]_p0t0_CellMask.tiff'
            }
        ]
    },
    ### 5. Lung LUNG IO TMA P#1A_Core[1,4,F]
    {
        ### Annotations (ground truth)
        'name': 'Lung_Tumors',
        'annot_root_path' : f'{annot_dir}/2022_Lung_Tumors',
        'roi_names' : [
            'LUNG IO TMA P#1A_Core[1,4,F]_[58113,12778]_p0t0_ROI1'
        ],
        ### Segmentations
        'seg_opts' :[
            {
                'type': 'cellXpress1',
                'filename' : f'{cX_proj_dir}/2022_07_cX2_Lung_Tumors/220203_Lung_Tumors_Single_Project/markerset0/segmentation/0/seg2.cXd',
                'frame_idx': 0
            },
            {
                'type': 'cX2_merged',
                'filename' : f'{bench_root}/Software_Projects/cX2/results_final/Lung_tumors_[1,4,F]_merged_cell_masks.tif',
            },
            {
                'type': 'qupath',
                'filename' : f'{bench_root}/Software_Projects/QuPath/3_Lung_Tumors/rois_final/1A_Core[1,4,F]_[58113,12778]_p0t0-labels.ome.tif'
            },
            {
                'type': 'deepcell',
                'filename' : f'{bench_root}/Software_Projects/DeepCell/results_final/1A_Core[1,4,F]_mask.tif'
            },
            {
                'type': 'cellprofiler',
                'filename' : f'{bench_root}/Software_Projects/CellProfiler/results_final/1A_Core[1,4,F]_[58113,12778]_p0t0_CellMask.tiff'
            }
        ]
    },
    ### 6. Lung LUNG IO TMA P#1C_Core[1,6,C]
    {
        ### Annotations (ground truth)
        'name': 'Lung_Tumors',
        'annot_root_path' : f'{annot_dir}/2022_Lung_Tumors',
        'roi_names' : [
            'LUNG IO TMA P#1C_Core[1,6,C]_[56306,18521]_p0t0_ROI2'
        ],
        ### Segmentations
        'seg_opts' :[
            {
                'type': 'cellXpress1',
                'filename' : f'{cX_proj_dir}/2022_07_cX2_Lung_Tumors/220203_Lung_Tumors_Single_Project/markerset0/segmentation/0/seg9.cXd',
                'frame_idx': 0
            },
            {
                'type': 'cX2_merged',
                'filename' : f'{bench_root}/Software_Projects/cX2/results_final/Lung_tumors_[1,6,C]_merged_cell_masks.tif',
            },
            {
                'type': 'qupath',
                'filename' : f'{bench_root}/Software_Projects/QuPath/3_Lung_Tumors/rois_final/1C_Core[1,6,C]_[56306,18521]_p0t0-labels.ome.tif'
            },
            {
                'type': 'deepcell',
                'filename' : f'{bench_root}/Software_Projects/DeepCell/results_final/1C_Core[1,6,C]_mask.tif'
            },
            {
                'type': 'cellprofiler',
                'filename' : f'{bench_root}/Software_Projects/CellProfiler/results_final/1C_Core[1,6,C]_[56306,18521]_p0t0_CellMask.tiff'
            }
        ]
    },    
    ### 7. COVID Lung High C1q
    {
        'name': 'COVID_Lung_Tissues',
        'annot_root_path' : f'{annot_dir}/2022_COVID_Lung_Tissues',
        'roi_names' : [
            'High C1q_stacked_ROI1',
            'High C1q_stacked_ROI2'
        ],
        ### Segmentations
        'seg_opts' :[
            {
                'type': 'cellXpress1',
                'filename' : f'{cX_proj_dir}/2022_07_cX2_COVID_Lung_Tissues/220726_COVID19_Lung_Tissues_Project/markerset1/segmentation/0/seg0.cXd',
                'frame_idx': 0
            },
            {
                'type': 'cellXpress2',
                'filename' : f'{cX_proj_dir}/2022_07_cX2_COVID_Lung_Tissues/220726_COVID19_Lung_Tissues_Project/markerset0/segmentation/0/seg0.cXd',
                'frame_idx': 0
            },
            {
                'type': 'qupath',
                'filename' : f'{bench_root}/Software_Projects/QuPath/2_COVID_Lung_Tissues/rois_final/High_C1q-A01_p1t0-labels.ome.tif'
            },
            {
                'type': 'deepcell',
                'filename' : f'{bench_root}/Software_Projects/DeepCell/results_final/MARIO_C1high_mask.tif'
            },
            {
                'type': 'cellprofiler',
                'filename' : f'{bench_root}/Software_Projects/CellProfiler/results_final/High_C1q-A01_p1t0.ome_CellMask.tiff'
            }
        ]
    },
    ### 8. COVID Lung Low C1q
    {
        ### Annotations (ground truth)
        'name': 'COVID_Lung_Tissues',
        'annot_root_path' : f'{annot_dir}/2022_COVID_Lung_Tissues',
        'roi_names' : [
            'Low C1q_stacked_ROI1',
            'Low C1q_stacked_ROI2'
        ],
        ### Segmentations
        'seg_opts' :[
            {
                'type': 'cellXpress1',
                'filename' : f'{cX_proj_dir}/2022_07_cX2_COVID_Lung_Tissues/220726_COVID19_Lung_Tissues_Project/markerset1/segmentation/0/seg1.cXd',
                'frame_idx': 0
            },
            {
                'type': 'cellXpress2',
                'filename' : f'{cX_proj_dir}/2022_07_cX2_COVID_Lung_Tissues/220726_COVID19_Lung_Tissues_Project/markerset0/segmentation/0/seg1.cXd',
                'frame_idx': 0
            },
            {
                'type': 'qupath',
                'filename' : f'{bench_root}/Software_Projects/QuPath/2_COVID_Lung_Tissues/rois_final/LowC1q-A02_p1t0-labels.ome.tif',
            },
            {
                'type': 'deepcell',
                'filename' : f'{bench_root}/Software_Projects/DeepCell/results_final/MARIO_C1low_mask.tif'
            },
            {
                'type': 'cellprofiler',
                'filename' : f'{bench_root}/Software_Projects/CellProfiler/results_final/LowC1q-A02_p1t0.ome_CellMask.tiff'
            }
        ]
    },
    ### 9. Tonsil Region 1
    {
        ### Annotations (ground truth)
        'name': 'Tonsil_Tissues',
        'annot_root_path' : f'{annot_dir}/2022_Tonsil_Tissues',
        'roi_names' : [
            'Region1_Tonsil_CODEX_ROI'
        ],
        ### Segmentations
        'seg_opts' :[
            {
                'type': 'cellXpress1',
                'filename' : f'{cX_proj_dir}/2022_07_cX2_Tonsil_Tissues/220311_Tonsil_Tissues_Single_Project/markerset0/segmentation/0/seg0.cXd',
                'frame_idx': 0
            },
            {
                'type': 'cX2_merged',
                'filename' : f'{bench_root}/Software_Projects/cX2/results_final/Tonsil_tissues_A01_merged_cell_masks.tif',
                'frame_idx': 0
            },
            {
                'type': 'qupath',
                'filename' : f'{bench_root}/Software_Projects/QuPath/4_Tonsil_Tissues/rois_final/Region1_A01_p0t0-labels.ome.tif',
            },
            {
                'type': 'deepcell',
                'filename' : f'{bench_root}/Software_Projects/DeepCell/results_final/Tonsil_crop1_mask.tif'
            },
            {
                'type': 'cellprofiler',
                'filename' : f'{bench_root}/Software_Projects/CellProfiler/results_final/Region1_A01_p0t0_CellMask.tiff'
            }
        ]
    },
    ### 10. Tonsil Region 2
    {
        ### Annotations (ground truth)
        'name': 'Tonsil_Tissues',
        'annot_root_path' : f'{annot_dir}/2022_Tonsil_Tissues',
        'roi_names' : [
            'Region2_Tonsil_CODEX_ROI1',
            'Region2_Tonsil_CODEX_ROI2'
        ],
        ### Segmentations
        'seg_opts' :[
            {
                'type': 'cellXpress1',
                'filename' : f'{cX_proj_dir}/2022_07_cX2_Tonsil_Tissues/220311_Tonsil_Tissues_Single_Project/markerset0/segmentation/0/seg1.cXd',
                'frame_idx': 0
            },
            {
                'type': 'cX2_merged',
                'filename' : f'{bench_root}/Software_Projects/cX2/results_final/Tonsil_tissues_A02_merged_cell_masks.tif',
                'frame_idx': 0
            },
            {
                'type': 'qupath',
                'filename' : f'{bench_root}/Software_Projects/QuPath/4_Tonsil_Tissues/rois_final/Region2_A02_p0t0-labels.ome.tif',
            },
            {
                'type': 'deepcell',
                'filename' : f'{bench_root}/Software_Projects/DeepCell/results_final/Tonsil_crop2_mask.tif'
            },
            {
                'type': 'cellprofiler',
                'filename' : f'{bench_root}/Software_Projects/CellProfiler/results_final/Region2_A02_p0t0_CellMask.tiff'
            }
        ]
    }
]

### Init output data
dataset_names = []
benchmark_idxs = []
roi_names = []
seg_types = []
accuracies = []
f1_scores = []
mean_matched_scores = []
panoptic_qualities = []

### Loop through all the benchmarks
for benchmark_idx, benchmark in enumerate(benchmarks):
    
    print(f'Benchmark image #{benchmark_idx}...')

    ### Loop through all the ROIs
    perf_scores = []

    for seg_opt in benchmark['seg_opts']:
        ### Load the segmentation mask
        if seg_opt['type'].startswith('cellXpress'):
            segmask = LoadcXSeg(seg_opt)
        else:
            segmask = LoadTIFFSeg(seg_opt)

        ### Make sure segmask is three dimension
        if len(segmask.shape) == 2:
            segmask = segmask[np.newaxis, :, :]

        ### Loop through seg_opts
        perf_scores_tmp = []

        for roi_name in benchmark['roi_names']:
            ### Load the ROI
            rois, roi_all_bbox, roi_def = LoadIJROI(
                benchmark['annot_root_path'], roi_name
            )
            
            ### Compute the performance
            perf_scores_cur, segmask_colored = ComputeSegPerf(
                rois, segmask, roi_def
            )

            ### Update the results
            dataset_names.append(benchmark['name'])
            benchmark_idxs.append(benchmark_idx)
            roi_names.append(roi_name)
            seg_types.append(seg_opt['type'])
            accuracies.append(perf_scores_cur['accuracy'])
            f1_scores.append(perf_scores_cur['f1'])
            mean_matched_scores.append(perf_scores_cur['mean_matched_score'])
            panoptic_qualities.append(perf_scores_cur['panoptic_quality'])

            ### Determine the segmentation performance
            perf_scores_tmp.append(perf_scores_cur)

            ### Save the output
            img_outline = DrawROISegmask(rois, segmask, roi_def)
            
            cv2.imwrite(
                f"{bench_root}/figures/{benchmark['name']}_{roi_name}_{seg_opt['type']}.png",
                img_outline[roi_def.top:roi_def.bottom,roi_def.left:roi_def.right]
                #img_outline
            )

            ### Save the ROI only
            img_roi = DrawROI(segmask.shape[1], segmask.shape[2], rois)
            img_roi = cv2.cvtColor(img_roi, cv2.COLOR_RGB2BGR)

            cv2.imwrite(
                f"{bench_root}/figures/{benchmark['name']}_{roi_name}_.png",
                img_roi
            )

        ### Add the results
        perf_scores.append(
            {k: [dic[k] for dic in perf_scores_tmp] for k in perf_scores_tmp[0]}
        )

    #pprint.pprint(perf_scores)

### Save the data
perf_data = pd.DataFrame(
    {
        "Dataset_name": dataset_names,
        "ROI_name": roi_names,
        "Seg_method": seg_types,
        "Accuracy": accuracies,
        "F1_score": f1_scores,
        "MMS": mean_matched_scores,
        "PQ": panoptic_qualities
    }
)
perf_data.to_csv(f"{bench_root}/results/performance.csv", index=False)

