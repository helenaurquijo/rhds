# rhds

## Setup

To setup this project 

```
cp config-template.env config.env
```

Update the paths in `config.env`.


## Installation

Use `renv` to install relevant packages. In R:

```
install.packages("renv")
renv::restore()
```

## How to run 

### 1. Download the data

```
bash scripts/download-data.sh
Rscript scripts/download-pan-cancer-clinical.r
```


### 2. Extract / clean relevant data

```
Rscript scripts/extract-data.r
```

### 3. Clean clinical data

```
Rscript scripts/clean-clinical.r
```

### 4. Predict proteins

```
Rscript scripts/predict-proteins.r
```

### 5. Combine data

```
Rscript scripts/combine.r
```

### 6. Analyse data

```
Rscript scripts/analysis.r
```

