# TEA-BISCUIT Workshop Day 2: Bulk Retrotranscriptomic Analysis


## Schedule
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


## Setup

### Option 1

Log into AWS as per Matthew Bendall's tutorial from [yesterday](https://github.com/nixonlab/teabiscuit), and copy this entire directory into your home directory. All you *really* need are the the `results`, `refs`, and `fgsea` directories, but if you copy the entire directory (as suggested), you'll be able to access intermediate files in case one of the commands fails (I hope not though!)

```
cp -R /efs/projects/tea-biscuit-retro . 
cd tea-biscuit-retro
```

### Option 2:

Ask Bhavya for *secret link* to download (sorry, cannot post the link here).

## How to Interact

### Option 1

You can create your own RScript and copy-paste commands from the `tutorial.md` or `.html` file. Or, you can open the tutorial.Rmd and just run things with me. 

Ceate your own R script: `File -> New File -> R Script`. Follow along by copying and pasting from the `tutorial.md` file.

### Option 2

You can simply follow along and watch. There is no pressure! Play along however you see fit. 
