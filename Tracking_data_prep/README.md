# Noradrenaline_proj
##Tracking_data_prep

parse_into_sessions.R
- Script that uses pre-recorded onsets and offsets (frame#) for text session of behavioural expirement to subset data 

Convert_to_DLC_format.R
- Using processing scripts that require data to have column names formated as would be outputted from DLC
- Converting SLEAP format to DLC before processing 

SLEAP_data_cleaning.R 
- Filter predictions with low prediction accuracy and interpolated missing frames, compute movement metrics of speed, acceleration ect. 

SLEAP_combine_data.R 
- Combine individual tracking files into one large df with metadata 
