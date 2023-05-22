# TEA-BISCUIT Workshop Day 2: Bulk Retrotranscriptomic Analysis


## Schedule:
Tuesday May 23 – Practicals I: Downstream Bulk Retrotranscriptomics

| Time | Event | Speaker | 
| ------------- | ------------- | ------------- |
| 9:00 AM - 9:20 AM | Coffee + Log in to RStudio | Dr. Matthew Bendall | 
| 9:30 AM - 10:30 AM | Unsupervised Principal Component Analysis (PCA) | Bhavya |
| 10:30 AM - 11:00 AM | Coffee break | Coffee speaks for itself |
| 11:00 - 12:00 AM | Differential expression analysis using DESeq2 | Bhavya |
| 12:00 - 1:00 PM | Performing Gene-set enrichment analysis | Bhavya |
| 1:00 PM - 2:30 PM | Lunch | |
| 2:30 PM - 3:30 PM | Performing rTWAS | Dr. Rodrigo Duarte | 
| 3:30 PM - 4:00 PM | Indian Tea and Biscuits | Parle-G |
| 4:00 PM - 4:30 PM | Performing rTWAS | Dr. Rodrigo Duarte | 
| 8:00 PM - TBD | Meet at Sour Mouse (East Village) for drinks and games (optional) | |


## Setup:

Option 1:

Log into AWS as per Matthew Bendall's tutorial from [yesterday](https://github.com/nixonlab/teabiscuit), and symlink this entire directory into your home directory. All you *really* need are the the `results` and `fgsea` directories, but if you symlink the entire directory (as suggested), you'll be able to access intermediate files in case one of the commands fails (I hope not though!)

```
cd /fsx/users/[youruid]
ln -s /efs/projects/tea-biscuit-retro . 
```

Option 2:

Ask Bhavya for *secret link* to download (sorry, cannot post the link here).
